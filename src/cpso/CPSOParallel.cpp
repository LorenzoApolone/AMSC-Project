/**
 * @file CPSOParallel.cpp
 * @brief Implementation of the CPSOParallel class methods inheriting from CPSOBase
 */
#include "CPSOParallel.hpp"
#include "NumericValidation.hpp"
#include "SubSwarmOwnershipUtils.hpp"
#include "SwarmMetrics.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#if __has_include(<mpi.h>)
#include <mpi.h>
#define CPSO_HAVE_MPI 1
#else
#define CPSO_HAVE_MPI 0
#endif

namespace {

// Aggregates time spent in the MPI primitives used by the solver.
struct MpiTimingStats {
  double allgather_s = 0.0;
  double bcast_s = 0.0;
  double allreduce_s = 0.0;
  double barrier_s = 0.0;

  double total() const { return allgather_s + bcast_s + allreduce_s; }
  double wait_total() const { return barrier_s; }
};

struct SyncedFailureState {
  bool failed = false;
  int rank = -1;
  std::string message;
};

struct ParallelMergeResult {
  std::vector<double> context;
  double fitness = std::numeric_limits<double>::infinity();
  bool failed = false;
  std::string failure_message;
};

#if CPSO_HAVE_MPI
// Fn is a callable and we measure how the mpi calls is doing.
template <typename Fn>
void time_mpi_call(double &accumulator, Fn &&fn) {
  const double start = MPI_Wtime();
  fn();
  accumulator += MPI_Wtime() - start;
}

double allreduce_max_double(double local_value) {
  double global_value = local_value;
  MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  return global_value;
}

// Used to synchronize the fact that one MPI-process failed and notified to all other ranks
// SyncedFailureState is used to store the failure, the rank and the message (see above)
SyncedFailureState sync_failure_state(bool local_failed,
                                      const std::string &local_message,
                                      int mpi_rank,
                                      MpiTimingStats &mpi_timings) {
  SyncedFailureState state;

  // Checking if something failed globally
  int local_failed_int = local_failed ? 1 : 0;
  int global_failed_int = 0;
  time_mpi_call(mpi_timings.allreduce_s, [&]() {
    MPI_Allreduce(&local_failed_int, &global_failed_int, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
  });

  if (global_failed_int == 0) {
    return state;
  }

  // Since something failed, let's find the rank who failed
  state.failed = true;
  int local_failure_rank = local_failed ? mpi_rank : std::numeric_limits<int>::max();
  int failure_rank = local_failure_rank;
  time_mpi_call(mpi_timings.allreduce_s, [&]() {
    MPI_Allreduce(&local_failure_rank, &failure_rank, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
  });

  // Then the message error
  int message_size =
      (local_failed && mpi_rank == failure_rank)
          ? static_cast<int>(local_message.size())
          : 0;
  time_mpi_call(mpi_timings.bcast_s, [&]() {
    MPI_Bcast(&message_size, 1, MPI_INT, failure_rank, MPI_COMM_WORLD);
  });

  // Broadcast the message to all the mpi ranks
  std::string broadcast_message(static_cast<size_t>(std::max(message_size, 0)),
                                '\0');
  if (local_failed && mpi_rank == failure_rank && message_size > 0) {
    broadcast_message = local_message;
  }

  if (message_size > 0) {
    time_mpi_call(mpi_timings.bcast_s, [&]() {
      MPI_Bcast(broadcast_message.data(), message_size, MPI_CHAR, failure_rank,
                MPI_COMM_WORLD);
    });
  }

  state.rank = failure_rank;
  state.message = broadcast_message.empty()
                      ? "parallel CPSO numeric failure"
                      : broadcast_message;
  return state;
}

#else
template <typename Fn>
void time_mpi_call(double &, Fn &&fn) {
  fn();
}

double allreduce_max_double(double local_value) { return local_value; }

SyncedFailureState sync_failure_state(bool local_failed,
                                      const std::string &local_message,
                                      int,
                                      MpiTimingStats &) {
  SyncedFailureState state;
  if (local_failed) {
    state.failed = true;
    state.rank = 0;
    state.message = local_message;
  }
  return state;
}
#endif

// Apply changes to the global context vector
// Parallel merges exchange only the coordinates of the active dims of a sub-swarm instead of the full context.
void apply_sparse_delta(std::vector<double> &target,
                        const std::vector<int> &active_dims,
                        const std::vector<double> &packed_delta_values,
                        size_t packed_offset = 0) {
  if (packed_delta_values.size() < packed_offset + active_dims.size()) {
    throw std::invalid_argument(
        "packed delta buffer is too small for the active dimensions");
  }

  for (size_t d = 0; d < active_dims.size(); ++d) {
    const int dim = active_dims[d];
    if (dim < 0 || static_cast<size_t>(dim) >= target.size()) {
      throw std::out_of_range("packed delta refers to an invalid dimension");
    }

    const double delta = packed_delta_values[packed_offset + d];
    ensure_finite_value(delta, "packed subswarm delta");
    target[dim] += delta;
    ensure_finite_value(target[dim], "context coordinate after packed delta");
  }
}

// Function to calculate the actual diff of the evaluation of the particle in order to not send all the context vector
std::vector<double> pack_swarm_delta(const SubSwarm &swarm,
                                     const std::vector<double> &base_context,
                                     int padded_size) {
  if (padded_size < 0) {
    throw std::invalid_argument("packed delta size must be >= 0");
  }

  const auto &active_dims = swarm.get_active_dims();
  const auto &best_pos = swarm.get_gbest_pos();
  if (active_dims.size() != best_pos.size()) {
    throw std::runtime_error(
        "subswarm active dimensions and best position are inconsistent");
  }
  if (padded_size < static_cast<int>(active_dims.size())) {
    throw std::invalid_argument(
        "packed delta size must cover the active dimensions");
  }

  ensure_finite_vector(base_context, "base context before delta packing");
  if (!is_finite_value(swarm.get_gbest_val())) {
    return std::vector<double>(static_cast<size_t>(padded_size), 0.0);
  }
  ensure_finite_vector(best_pos, "subswarm best position before delta packing");

  std::vector<double> packed(static_cast<size_t>(padded_size), 0.0);
  for (size_t d = 0; d < active_dims.size(); ++d) {
    const int dim = active_dims[d];
    if (dim < 0 || static_cast<size_t>(dim) >= base_context.size()) {
      throw std::out_of_range("active dimension out of bounds while packing delta");
    }

    const double delta = best_pos[d] - base_context[dim];
    ensure_finite_value(delta, "packed subswarm delta");
    packed[d] = delta;
  }

  return packed;
}

// Finding a new global candidate for the global context vector
std::vector<double> build_batch_candidate_vector(
    const std::vector<double> &base_context, int batch_index,
    const std::vector<std::pair<int, int>> &swarm_ranges,
    const std::vector<SubSwarm> &swarms,
    const std::vector<double> &gathered_deltas, int packed_delta_size,
    const std::vector<int> *include_flags = nullptr) {
  std::vector<double> candidate = base_context;
  ensure_finite_vector(candidate, "base context before gathered merge");

  for (int rank = 0; rank < static_cast<int>(swarm_ranges.size()); ++rank) {
    const int swarm_idx = compute_batch_swarm_index(swarm_ranges, rank, batch_index);
    if (swarm_idx < 0) {
      continue;
    }
    if (include_flags != nullptr) {
      if (rank >= static_cast<int>(include_flags->size()) || (*include_flags)[rank] == 0) {
        continue;
      }
    }

    apply_sparse_delta(candidate, swarms[swarm_idx].get_active_dims(),
                       gathered_deltas,
                       static_cast<size_t>(rank * packed_delta_size));
  }

  return candidate;
}

double evaluate_candidate_fitness(const TestFunction &f,
                                  const std::vector<double> &candidate,
                                  const char *context_label) {
  ensure_finite_vector(candidate, context_label);
  return sanitize_fitness(f.value(candidate));
}

bool try_evaluate_candidate_fitness(const TestFunction &f,
                                    const std::vector<double> &candidate,
                                    const char *context_label,
                                    double &fitness_out,
                                    std::string &failure_message) {
  try {
    fitness_out = evaluate_candidate_fitness(f, candidate, context_label);
    return true;
  } catch (const std::exception &ex) {
    fitness_out = std::numeric_limits<double>::infinity();
    failure_message = ex.what();
    return false;
  }
}

// Calculating all the mpi timings
void fill_parallel_timing_artifacts(CpsoRunArtifacts &artifacts,
                                    const MpiTimingStats &mpi_timings,
                                    double optimization_wall_time_s) {
  const double local_comm_total = mpi_timings.total();
  const double local_wait_total = mpi_timings.wait_total();
  const double local_compute_total =
      std::max(0.0, optimization_wall_time_s - local_comm_total - local_wait_total);

  artifacts.comm_total_s = allreduce_max_double(local_comm_total);
  artifacts.comm_allgather_s = allreduce_max_double(mpi_timings.allgather_s);
  artifacts.comm_bcast_s = allreduce_max_double(mpi_timings.bcast_s);
  artifacts.comm_allreduce_s = allreduce_max_double(mpi_timings.allreduce_s);
  artifacts.comm_barrier_s = allreduce_max_double(mpi_timings.barrier_s);
  artifacts.wait_total_s = allreduce_max_double(local_wait_total);
  artifacts.compute_total_s = allreduce_max_double(local_compute_total);
}

CpsoRunArtifacts build_parallel_failure_artifacts(
    const ContextVector &context, int mpi_size, int iter,
    const std::vector<double> &fitness_history, const std::string &message,
    const MpiTimingStats &mpi_timings, double optimization_wall_time_s) {
  CpsoRunArtifacts artifacts;
  artifacts.best_position = context.get_full_vector();
  artifacts.best_fitness = std::numeric_limits<double>::infinity();
  artifacts.best_fitness_history = fitness_history;
  artifacts.cores = mpi_size;
  artifacts.iterations = iter;
  artifacts.stop_reason = CpsoStopReason::NUMERIC_FAILURE;
  artifacts.failure_message = message;
  fill_parallel_timing_artifacts(artifacts, mpi_timings,
                                 optimization_wall_time_s);
  return artifacts;
}

bool any_rank_has_candidate(const std::vector<int> &flags) {
  return std::any_of(flags.begin(), flags.end(),
                     [](int flag) { return flag != 0; });
}

bool env_flag_enabled(const char *name) {
  const char *raw_value = std::getenv(name);
  if (raw_value == nullptr || *raw_value == '\0') {
    return false;
  }

  std::string value(raw_value);
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });
  return value == "1" || value == "true" || value == "yes" ||
         value == "on";
}

bool skip_parallel_greedy_merge_fallback() {
  static const bool disabled =
      env_flag_enabled("CPSO_MPI_DISABLE_GREEDY_MERGE") ||
      env_flag_enabled("CPSO_PARALLEL_DISABLE_GREEDY_MERGE");
  return disabled;
}

// Incremental merge based on best-evaluated local deltas
ParallelMergeResult parallel_greedy_merge(
    const TestFunction &f, const std::vector<double> &base_context,
    double base_fitness, const std::vector<SubSwarm> &swarms,
    const std::vector<std::pair<int, int>> &swarm_ranges, int batch_index,
    std::vector<double> local_candidate_delta_values, bool has_local_candidate,
    int packed_delta_size, int mpi_rank, int mpi_size,
    MpiTimingStats &mpi_timings) {
  ParallelMergeResult result;
  std::vector<double> current_greedy = base_context;
  double current_greedy_fitness = base_fitness;
  const int local_swarm_idx =
      compute_batch_swarm_index(swarm_ranges, mpi_rank, batch_index);

#if CPSO_HAVE_MPI
  std::vector<double> gathered_fitnesses(
      mpi_size, std::numeric_limits<double>::infinity());
#endif

  // Accept only candidate deltas that improve the current batch context.
  while (true) {
    double local_candidate_fitness = std::numeric_limits<double>::infinity();
    bool local_candidate_failed = false;
    std::string local_failure_message;

    if (has_local_candidate && local_swarm_idx >= 0) {
      try {
        std::vector<double> test_greedy = current_greedy;
        apply_sparse_delta(test_greedy, swarms[local_swarm_idx].get_active_dims(),
                           local_candidate_delta_values);
        local_candidate_fitness =
            evaluate_candidate_fitness(f, test_greedy, "greedy merged context");
      } catch (const std::exception &ex) {
        local_candidate_failed = true;
        local_failure_message =
            std::string("parallel greedy merge batch ") +
            std::to_string(batch_index) + ": " + ex.what();
      }
    }

    const SyncedFailureState local_failure =
        sync_failure_state(local_candidate_failed, local_failure_message,
                           mpi_rank, mpi_timings);
    if (local_failure.failed) {
      result.failed = true;
      result.failure_message = local_failure.message;
      return result;
    }

#if CPSO_HAVE_MPI
    time_mpi_call(mpi_timings.allgather_s, [&]() {
      MPI_Allgather(&local_candidate_fitness, 1, MPI_DOUBLE,
                    gathered_fitnesses.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
    });

    int selected_rank = -1;
    double selected_fitness = current_greedy_fitness;
    for (int rank = 0; rank < mpi_size; ++rank) {
      if (gathered_fitnesses[rank] < selected_fitness) {
        selected_fitness = gathered_fitnesses[rank];
        selected_rank = rank;
      }
    }

    if (selected_rank < 0) {
      break;
    }

    const int selected_swarm_idx =
        compute_batch_swarm_index(swarm_ranges, selected_rank, batch_index);
    if (selected_swarm_idx < 0) {
      result.failed = true;
      result.failure_message =
          "selected MPI rank does not own a swarm in this batch";
      return result;
    }

    std::vector<double> selected_delta_values(
        static_cast<size_t>(packed_delta_size), 0.0);
    if (mpi_rank == selected_rank) {
      selected_delta_values = local_candidate_delta_values;
      has_local_candidate = false;
      std::fill(local_candidate_delta_values.begin(),
                local_candidate_delta_values.end(), 0.0);
    }

    time_mpi_call(mpi_timings.bcast_s, [&]() {
      MPI_Bcast(selected_delta_values.data(), packed_delta_size, MPI_DOUBLE,
                selected_rank, MPI_COMM_WORLD);
    });

    bool local_apply_failed = false;
    std::string local_apply_failure_message;
    try {
      apply_sparse_delta(current_greedy,
                         swarms[selected_swarm_idx].get_active_dims(),
                         selected_delta_values);
      current_greedy_fitness = selected_fitness;
    } catch (const std::exception &ex) {
      local_apply_failed = true;
      local_apply_failure_message =
          std::string("parallel greedy apply batch ") +
          std::to_string(batch_index) + ": " + ex.what();
    }

    const SyncedFailureState apply_failure =
        sync_failure_state(local_apply_failed, local_apply_failure_message,
                           mpi_rank, mpi_timings);
    if (apply_failure.failed) {
      result.failed = true;
      result.failure_message = apply_failure.message;
      return result;
    }
#else
    (void)mpi_rank;
    (void)mpi_size;

    if (!has_local_candidate || local_swarm_idx < 0 ||
        local_candidate_fitness >= current_greedy_fitness) {
      break;
    }

    try {
      apply_sparse_delta(current_greedy, swarms[local_swarm_idx].get_active_dims(),
                         local_candidate_delta_values);
      current_greedy_fitness = local_candidate_fitness;
    } catch (const std::exception &ex) {
      result.failed = true;
      result.failure_message = ex.what();
      return result;
    }
    has_local_candidate = false;
#endif
  }

  result.context = std::move(current_greedy);
  result.fitness = current_greedy_fitness;
  return result;
}

} // namespace

// Constructor for uniform topology
CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           NetworkType topology, int shuffle_freq,
                           int stagnation_patience, double w_start,
                           double w_end, double c1, double c2,
                           unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology,
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2,
               seed) {}

CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           const SubSwarmTopologyConfig &topology_config,
                           int shuffle_freq, int stagnation_patience,
                           double w_start, double w_end, double c1, double c2,
                           unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology_config,
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2,
               seed) {}

// Constructor for different types of topologies
CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           const std::vector<NetworkType> &topologies,
                           int shuffle_freq, int stagnation_patience,
                           double w_start, double w_end, double c1, double c2,
                           unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies,
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2,
               seed) {}

CPSOParallel::CPSOParallel(
    int k_subswarms, int num_particles_per_swarm,
    const std::vector<SubSwarmTopologyConfig> &topologies, int shuffle_freq,
    int stagnation_patience, double w_start, double w_end, double c1,
    double c2, unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2, seed) {}

CpsoRunArtifacts CPSOParallel::run_optimization_loop(
    const TestFunction &f, StoppingCriteriaManager &stop_manager,
    std::vector<SubSwarm> &swarms, std::vector<std::mt19937> &gens,
    ContextVector &context, std::mt19937 &global_gen) {

  int mpi_rank = 0, mpi_size = 1;
#if CPSO_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
  const bool disable_greedy_merge_fallback =
      skip_parallel_greedy_merge_fallback();

  // mpi_timings accumulates time spent in MPI primitives.
  // optimization_start anchors the wall-clock duration of the run.
  // fitness_history stores the best global fitness after each iteration.
  MpiTimingStats mpi_timings;
  const auto optimization_start = std::chrono::steady_clock::now();
  std::vector<double> fitness_history;
  int iter = 0;

  // Helper used to build an artifact for the numerical_failure
  auto build_numeric_failure = [&](const std::string &message) {
    return build_parallel_failure_artifacts(
        context, mpi_size, iter, fitness_history, message, mpi_timings,
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      optimization_start)
            .count());
  };

  try {
    // Derive the static ownership layout used during the whole run.Across MPI ranks,
    // sub-swarms themselves are also distributed contiguously, so every rank
    const int total_dim = f.dim;
    const int dims_per_swarm = total_dim / get_num_subswarms();
    const int remainder = total_dim % get_num_subswarms();
    const int max_swarm_dims = dims_per_swarm + (remainder > 0 ? 1 : 0);
    const std::vector<std::pair<int, int>> swarm_ranges =
        compute_all_swarm_ranges(get_num_subswarms(), mpi_size);
    const std::pair<int, int> local_range =
        (mpi_rank >= 0 && mpi_rank < static_cast<int>(swarm_ranges.size()))
            ? swarm_ranges[mpi_rank]
            : std::make_pair(0, 0);
    const int local_start_idx = local_range.first;
    const int local_end_idx = local_range.second;

    // Synchronize the starting context across all ranks.
    std::vector<double> init_vec = context.get_full_vector();
    double init_fitness = context.get_best_fitness();
    ensure_finite_vector(init_vec, "initial parallel context vector");
    ensure_finite_value(init_fitness, "initial parallel context fitness");

#if CPSO_HAVE_MPI
    time_mpi_call(mpi_timings.bcast_s, [&]() {
      MPI_Bcast(init_vec.data(), total_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    });
    time_mpi_call(mpi_timings.bcast_s, [&]() {
      MPI_Bcast(&init_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    });
#endif

    context.set_full_vector(init_vec, init_fitness);

    // The k-swarms are not equally distributed
    // global_max_local_swarms indicates the number of subswarm that every process has
    const int local_num_swarms = local_end_idx - local_start_idx;
    int global_max_local_swarms = 0;
#if CPSO_HAVE_MPI
    time_mpi_call(mpi_timings.allreduce_s, [&]() {
      MPI_Allreduce(&local_num_swarms, &global_max_local_swarms, 1, MPI_INT,
                    MPI_MAX, MPI_COMM_WORLD);
    });
#else
    global_max_local_swarms = local_num_swarms;
#endif

    // Initialize local sub-swarms and evaluates their first improvements.
    // The initialization itself is local, but any improved local gbest must be
    // folded back into a shared global context before the main loop starts.
    for (int b = 0; b < global_max_local_swarms; ++b) {

      // Snapshot the context at the beginning of this initialization batch.
      double current_base_fitness = context.get_best_fitness();
      std::vector<double> current_base_vec = context.get_full_vector();
      ensure_finite_vector(current_base_vec,
                           "parallel base context before initialization batch");
      ensure_finite_value(current_base_fitness,
                          "parallel base fitness before initialization batch");

      // Local contribution for this batch
      std::vector<double> local_init_delta(
          static_cast<size_t>(max_swarm_dims), 0.0);
      int local_has_init_candidate = 0;
      const int i = local_start_idx + b;

      // Only the rank that owns sub-swarm "i" initializes it here.
      bool local_init_failed = false;
      std::string local_init_failure_message;
      if (i < local_end_idx) {
        try {
          swarms[i].initialize(gens[i], context, f);
          if (swarms[i].get_gbest_val() < current_base_fitness) {
            local_init_delta =
                pack_swarm_delta(swarms[i], current_base_vec, max_swarm_dims);
            local_has_init_candidate = 1;
          }
        } catch (const std::exception &ex) {
          local_init_failed = true;
          local_init_failure_message =
              std::string("parallel initialization batch ") +
              std::to_string(b) + ": " + ex.what();
          local_has_init_candidate = 0;
          std::fill(local_init_delta.begin(), local_init_delta.end(), 0.0);
        }
      }

      const SyncedFailureState init_failure =
          sync_failure_state(local_init_failed, local_init_failure_message,
                             mpi_rank, mpi_timings);
      if (init_failure.failed) {
        return build_numeric_failure(init_failure.message);
      }

#if CPSO_HAVE_MPI
      // Gather from every rank a flag saying whether that rank found an improving initialized swarm and the delta improvements produced by the swarm
      std::vector<int> gathered_init_flags(static_cast<size_t>(mpi_size), 0);
      std::vector<double> gathered_init_deltas(
          static_cast<size_t>(mpi_size * max_swarm_dims), 0.0);

      time_mpi_call(mpi_timings.allgather_s, [&]() {
        MPI_Allgather(&local_has_init_candidate, 1, MPI_INT,
                      gathered_init_flags.data(), 1, MPI_INT,
                      MPI_COMM_WORLD);
      });
      time_mpi_call(mpi_timings.allgather_s, [&]() {
        MPI_Allgather(local_init_delta.data(), max_swarm_dims, MPI_DOUBLE,
                      gathered_init_deltas.data(), max_swarm_dims, MPI_DOUBLE,
                      MPI_COMM_WORLD);
      });

      if (any_rank_has_candidate(gathered_init_flags)) {
        // Firstly apply every improving local
        std::vector<double> synced_init_vector = build_batch_candidate_vector(
            current_base_vec, b, swarm_ranges, swarms, gathered_init_deltas,
            max_swarm_dims, &gathered_init_flags);

        double init_synced_best = std::numeric_limits<double>::infinity();
        bool init_eval_failed = false;
        std::string init_eval_failure_message;
        if (mpi_rank == 0) {
          init_eval_failed = !try_evaluate_candidate_fitness(
              f, synced_init_vector, "batched initialization context",
              init_synced_best, init_eval_failure_message);
          if (init_eval_failed) {
            init_eval_failure_message =
                std::string("parallel initialization merge batch ") +
                std::to_string(b) + ": " + init_eval_failure_message;
          }
        }
        const SyncedFailureState init_eval_failure =
            sync_failure_state(init_eval_failed, init_eval_failure_message,
                               mpi_rank, mpi_timings);
        if (init_eval_failure.failed) {
          return build_numeric_failure(init_eval_failure.message);
        }
        time_mpi_call(mpi_timings.bcast_s, [&]() {
          MPI_Bcast(&init_synced_best, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        });

        // If the merge helps, go on, otherwise fall back to the greedy merge
        if (init_synced_best < current_base_fitness) {
          context.set_full_vector(synced_init_vector, init_synced_best);
        } else {
          ParallelMergeResult merge_result = parallel_greedy_merge(
              f, current_base_vec, current_base_fitness, swarms, swarm_ranges,
              b, local_init_delta, local_has_init_candidate != 0,
              max_swarm_dims, mpi_rank, mpi_size, mpi_timings);
          if (merge_result.failed) {
            return build_numeric_failure(merge_result.failure_message);
          }
          context.set_full_vector(merge_result.context, merge_result.fitness);
        }
      }
#else
      // There is only a local improvement
      if (i < local_end_idx && local_has_init_candidate != 0) {
        std::vector<double> local_init_vector = current_base_vec;
        apply_sparse_delta(local_init_vector, swarms[i].get_active_dims(),
                           local_init_delta);
        double init_local_best = evaluate_candidate_fitness(
            f, local_init_vector, "local initialization context");
        if (init_local_best < current_base_fitness) {
          context.set_full_vector(local_init_vector, init_local_best);
        }
      }
#endif
    }

    // Main optimization loop
    int injection_count = 0;
    bool must_stop = false;
    double last_avg_distance = std::numeric_limits<double>::infinity();
    CpsoStopReason stop_reason = CpsoStopReason::UNKNOWN;

    while (!must_stop) {
      iter++;
      stop_manager.increment_iterations();

      // Periodically reshuffle the dimension assignment of the sub-swarms
      if (iter > 1 && iter % get_shuffle_freq() == 0) {
#if CPSO_HAVE_MPI
        time_mpi_call(mpi_timings.barrier_s, [&]() {
          MPI_Barrier(MPI_COMM_WORLD);
        });
#endif
        std::vector<int> permutation(total_dim);
        if (mpi_rank == 0) {
          std::iota(permutation.begin(), permutation.end(), 0);
          std::shuffle(permutation.begin(), permutation.end(), global_gen);
        }
#if CPSO_HAVE_MPI
        time_mpi_call(mpi_timings.bcast_s, [&]() {
          MPI_Bcast(permutation.data(), total_dim, MPI_INT, 0, MPI_COMM_WORLD);
        });
#endif

        int current_dim_start = 0;
        for (int swarm_idx = 0; swarm_idx < get_num_subswarms(); ++swarm_idx) {
          int swarm_dims = dims_per_swarm + (swarm_idx < remainder ? 1 : 0);
          std::vector<int> new_active_dims(static_cast<size_t>(swarm_dims));
          for (int d = 0; d < swarm_dims; ++d) {
            new_active_dims[d] = permutation[current_dim_start + d];
          }
          current_dim_start += swarm_dims;

          const bool is_owned =
              (swarm_idx >= local_start_idx && swarm_idx < local_end_idx);
          swarms[swarm_idx].update_active_dims(new_active_dims, context,
                                               gens[swarm_idx], is_owned);
        }
      }

      double progress_ratio =
          static_cast<double>(stop_manager.get_current_iters()) /
          stop_manager.get_max_iters();
      if (progress_ratio > 1.0) {
        progress_ratio = 1.0;
      }
      ensure_finite_value(progress_ratio, "parallel progress ratio");

      const double current_w =
          get_w_max() - (get_w_max() - get_w_min()) * progress_ratio;
      ensure_finite_value(current_w, "parallel inertia weight");

      // Update and Merge one batch slot at a time
      for (int b = 0; b < global_max_local_swarms; ++b) {
        // Save the shared context before any sub-swarm in this batch
        std::vector<double> base_context = context.get_full_vector();
        double base_fitness = context.get_best_fitness();
        ensure_finite_vector(base_context,
                             "parallel base context before update batch");
        ensure_finite_value(base_fitness,
                            "parallel base fitness before update batch");

        std::vector<double> local_full_delta(
            static_cast<size_t>(max_swarm_dims), 0.0);
        std::vector<double> local_salvaged_delta(
            static_cast<size_t>(max_swarm_dims), 0.0);
        int local_has_salvaged_candidate = 0;

        const int i = local_start_idx + b;
        bool local_update_failed = false;
        std::string local_update_failure_message;
        if (i < local_end_idx) {
          try {
            // Re-evaluate the local swarm against the latest global context vector apply one PSO update then checks its personal/global
            swarms[i].recalculate_fitness(context, f);
            swarms[i].update_velocities_and_positions(
                current_w, get_c1(), get_c2(), gens[i], progress_ratio);
            swarms[i].evaluate_and_update(context, f);

            local_full_delta =
                pack_swarm_delta(swarms[i], base_context, max_swarm_dims);
            if (swarms[i].get_gbest_val() < base_fitness) {
              local_salvaged_delta = local_full_delta;
              local_has_salvaged_candidate = 1;
            }
          } catch (const std::exception &ex) {
            local_update_failed = true;
            local_update_failure_message =
                std::string("parallel update batch ") + std::to_string(b) +
                ": " + ex.what();
            local_has_salvaged_candidate = 0;
            std::fill(local_full_delta.begin(), local_full_delta.end(), 0.0);
            std::fill(local_salvaged_delta.begin(), local_salvaged_delta.end(),
                      0.0);
          }
        }

        const SyncedFailureState update_failure =
            sync_failure_state(local_update_failed, local_update_failure_message,
                               mpi_rank, mpi_timings);
        if (update_failure.failed) {
          return build_numeric_failure(update_failure.message);
        }

#if CPSO_HAVE_MPI
        // Collect the raw local update from each rank, plus the stricter
        std::vector<double> gathered_full_deltas(
            static_cast<size_t>(mpi_size * max_swarm_dims), 0.0);
        std::vector<double> gathered_salvaged_deltas(
            static_cast<size_t>(mpi_size * max_swarm_dims), 0.0);
        std::vector<int> gathered_salvaged_flags(
            static_cast<size_t>(mpi_size), 0);
//@note A lot of all to all communication here, 
// Probably inevitable, but they slow down the parallization a lot.
// Maybe one can pack all data in a single allgather call to reduce the overhead?
        time_mpi_call(mpi_timings.allgather_s, [&]() {
          MPI_Allgather(local_full_delta.data(), max_swarm_dims, MPI_DOUBLE,
                        gathered_full_deltas.data(), max_swarm_dims,
                        MPI_DOUBLE, MPI_COMM_WORLD);
        });
        time_mpi_call(mpi_timings.allgather_s, [&]() {
          MPI_Allgather(local_salvaged_delta.data(), max_swarm_dims,
                        MPI_DOUBLE, gathered_salvaged_deltas.data(),
                        max_swarm_dims, MPI_DOUBLE, MPI_COMM_WORLD);
        });
        time_mpi_call(mpi_timings.allgather_s, [&]() {
          MPI_Allgather(&local_has_salvaged_candidate, 1, MPI_INT,
                        gathered_salvaged_flags.data(), 1, MPI_INT,
                        MPI_COMM_WORLD);
        });

        // Merge every updated local best in the batch and evaluate the resulting global context on rank 0.
        std::vector<double> new_full_vector = build_batch_candidate_vector(
            base_context, b, swarm_ranges, swarms, gathered_full_deltas,
            max_swarm_dims);

        double new_true_fitness = std::numeric_limits<double>::infinity();
        bool new_true_fitness_failed = false;
        std::string new_true_fitness_failure_message;
        if (mpi_rank == 0) {
          new_true_fitness_failed = !try_evaluate_candidate_fitness(
              f, new_full_vector, "batched parallel context",
              new_true_fitness, new_true_fitness_failure_message);
          if (new_true_fitness_failed) {
            new_true_fitness_failure_message =
                std::string("parallel context merge batch ") +
                std::to_string(b) + ": " + new_true_fitness_failure_message;
          }
        }
        const SyncedFailureState full_merge_failure =
            sync_failure_state(new_true_fitness_failed,
                               new_true_fitness_failure_message, mpi_rank,
                               mpi_timings);
        if (full_merge_failure.failed) {
          return build_numeric_failure(full_merge_failure.message);
        }
        time_mpi_call(mpi_timings.bcast_s, [&]() {
          MPI_Bcast(&new_true_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        });

        // Merge policy
        if (new_true_fitness < base_fitness) {
          context.set_full_vector(new_full_vector, new_true_fitness);
        } else if (any_rank_has_candidate(gathered_salvaged_flags)) {
          std::vector<double> final_salvaged_vector =
              build_batch_candidate_vector(base_context, b, swarm_ranges,
                                           swarms, gathered_salvaged_deltas,
                                           max_swarm_dims,
                                           &gathered_salvaged_flags);

          double final_fitness = std::numeric_limits<double>::infinity();
          bool final_fitness_failed = false;
          std::string final_fitness_failure_message;
          if (mpi_rank == 0) {
            final_fitness_failed = !try_evaluate_candidate_fitness(
                f, final_salvaged_vector, "salvaged parallel context",
                final_fitness, final_fitness_failure_message);
            if (final_fitness_failed) {
              final_fitness_failure_message =
                  std::string("parallel salvaged merge batch ") +
                  std::to_string(b) + ": " + final_fitness_failure_message;
            }
          }
          const SyncedFailureState salvaged_merge_failure =
              sync_failure_state(final_fitness_failed,
                                 final_fitness_failure_message, mpi_rank,
                                 mpi_timings);
          if (salvaged_merge_failure.failed) {
            return build_numeric_failure(salvaged_merge_failure.message);
          }
          time_mpi_call(mpi_timings.bcast_s, [&]() {
            MPI_Bcast(&final_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
          });

          if (final_fitness < base_fitness) {
            context.set_full_vector(final_salvaged_vector, final_fitness);
          } else if (!disable_greedy_merge_fallback) {
            ParallelMergeResult merge_result = parallel_greedy_merge(
                f, base_context, base_fitness, swarms, swarm_ranges, b,
                local_salvaged_delta, local_has_salvaged_candidate != 0,
                max_swarm_dims, mpi_rank, mpi_size, mpi_timings);
            if (merge_result.failed) {
              return build_numeric_failure(merge_result.failure_message);
            }
            context.set_full_vector(merge_result.context,
                                    merge_result.fitness);
          } else {
            // Experimental communication-light mode: keep the current shared
            // context when both full and salvaged merges fail, instead of
            // entering the extra greedy synchronization rounds.
          }
        }
#else
        // Apply only the local batch update and keep it if it improves the current context.
        if (i < local_end_idx) {
          std::vector<double> local_vector = base_context;
          apply_sparse_delta(local_vector, swarms[i].get_active_dims(),
                             local_full_delta);
          double local_best = evaluate_candidate_fitness(
              f, local_vector, "local parallel context");
          if (local_best < base_fitness) {
            context.set_full_vector(local_vector, local_best);
          }
        }
#endif
      }

      // After all batch merges, measure the new global state and update the stopping machinery.
      double current_best_fitness = context.get_best_fitness();
      const std::vector<double> &current_gbest_pos = context.get_full_vector();
      ensure_finite_value(current_best_fitness,
                          "current best parallel fitness");
      ensure_finite_vector(current_gbest_pos,
                           "current best parallel position");

      fitness_history.push_back(current_best_fitness);

      SwarmMetrics::compute_avg_distance(swarms, current_gbest_pos,
                                         last_avg_distance, local_start_idx,
                                         local_end_idx, true,
                                         &mpi_timings.allreduce_s);
      ensure_finite_value(last_avg_distance, "parallel average swarm distance");

      bool local_stop =
          stop_manager.should_stop(current_best_fitness, last_avg_distance);
      int current_stag_iters = stop_manager.get_current_stagnation_iters();
      const bool stop_for_max_iters = stop_manager.reached_max_iters();
      const bool stop_for_low_diversity =
          stop_manager.reached_diversity_limit(last_avg_distance);
      const bool stop_for_stagnation =
          local_stop && !stop_for_max_iters && !stop_for_low_diversity &&
          current_stag_iters >= stop_manager.get_max_stagnation_iters();
      const bool should_inject =
          !stop_for_max_iters && !stop_for_low_diversity &&
          current_stag_iters > 0 &&
          (current_stag_iters % get_stagnation_patience() == 0);

      // Injection velocities to avoid stopping in false positives
      if (should_inject) {
        injection_count++;
        bool hard_reset = (injection_count % 3 == 0);
        for (int swarm_idx = local_start_idx; swarm_idx < local_end_idx;
             ++swarm_idx) {
          swarms[swarm_idx].inject_velocities(gens[swarm_idx], hard_reset);
          if (hard_reset) {
            swarms[swarm_idx].reset_gbest_attraction();
          }
        }

        stop_manager.reset_stagnation();
        if (stop_for_stagnation) {
          local_stop = false;
        }
      }

#if CPSO_HAVE_MPI
      int local_stop_int = local_stop ? 1 : 0;
      int global_stop_int = 0;
      time_mpi_call(mpi_timings.allreduce_s, [&]() {
        MPI_Allreduce(&local_stop_int, &global_stop_int, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
      });
      must_stop = (global_stop_int > 0);
#else
      must_stop = local_stop;
#endif
      if (must_stop) {
        stop_reason = infer_cpso_stop_reason(stop_for_max_iters,
                                             stop_for_low_diversity,
                                             stop_for_stagnation);
      }
    }

    // Write into the artifact all the data for the benchmarks
    const double optimization_wall_time_s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      optimization_start)
            .count();

    CpsoRunArtifacts artifacts;
    artifacts.best_position = context.get_full_vector();
    artifacts.best_fitness = context.get_best_fitness();
    artifacts.best_fitness_history = fitness_history;
    artifacts.cores = mpi_size;
    artifacts.iterations = iter;
    fill_parallel_timing_artifacts(artifacts, mpi_timings,
                                   optimization_wall_time_s);
    artifacts.stop_reason = stop_reason;
    return artifacts;
  } catch (const std::exception &ex) {
    const SyncedFailureState synced_failure = sync_failure_state(
        true, std::string("parallel CPSO numeric failure: ") + ex.what(),
        mpi_rank, mpi_timings);
    return build_numeric_failure(synced_failure.message);
  }
}
