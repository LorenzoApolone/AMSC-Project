/**
 * @file CPSOBase.cpp
 * @brief Implements the shared setup logic used by the CPSO solvers.
 */
#include "CPSOBase.hpp"
#include "NumericValidation.hpp"
#include <chrono>
#include <stdexcept>

#if __has_include(<mpi.h>)
#include <mpi.h>
#define CPSO_HAVE_MPI 1
#else
#define CPSO_HAVE_MPI 0
#endif

namespace {

std::vector<SubSwarmTopologyConfig>
normalize_topology_configs(int k_subswarms,
                           const std::vector<SubSwarmTopologyConfig> &topologies,
                           int particles_per_swarm) {
  std::vector<SubSwarmTopologyConfig> normalized = topologies;
  if (normalized.size() > static_cast<size_t>(k_subswarms)) {
    normalized.resize(static_cast<size_t>(k_subswarms));
  }
  const SubSwarmTopologyConfig fallback =
      topologies.empty()
          ? SubSwarmTopologyConfig::from_type(NetworkType::SCALE_FREE)
          : topologies.front();

  while (normalized.size() < static_cast<size_t>(k_subswarms)) {
    normalized.push_back(fallback);
  }

  for (const auto &config : normalized) {
    validate_subswarm_topology_config(config, particles_per_swarm);
  }

  return normalized;
}

std::vector<SubSwarmTopologyConfig>
to_topology_configs(const std::vector<NetworkType> &topologies) {
  std::vector<SubSwarmTopologyConfig> configs;
  configs.reserve(topologies.size());
  for (NetworkType type : topologies) {
    configs.push_back(SubSwarmTopologyConfig::from_type(type));
  }
  return configs;
}

void validate_cpso_parameters(int k_subswarms, int num_particles_per_swarm,
                              int shuffle_freq, int stagnation_patience) {
  // Safety Checks
  if (k_subswarms <= 0) {
    throw std::invalid_argument("number of subswarms must be > 0");
  }
  if (num_particles_per_swarm <= 0) {
    throw std::invalid_argument("particles per swarm must be > 0");
  }
  if (shuffle_freq <= 0) {
    throw std::invalid_argument("shuffle frequency must be > 0");
  }
  if (stagnation_patience <= 0) {
    throw std::invalid_argument("stagnation patience must be > 0");
  }
}


bool try_build_finite_initial_fitness(const TestFunction &f,
                                      std::vector<double> &candidate,
                                      std::mt19937 &gen,
                                      double lower_bound,
                                      double upper_bound,
                                      double &fitness_out) {
  constexpr int max_attempts = 64;
  std::uniform_real_distribution<double> dist(lower_bound, upper_bound);

  // Retry a few random candidates before falling back to the domain midpoint.
  for (int i = 0; i < max_attempts; ++i) {
    const double fitness = sanitize_fitness(f.value(candidate));
    if (is_finite_value(fitness)) {
      fitness_out = fitness;
      return true;
    }

    for (double &coord : candidate) {
      coord = dist(gen);
    }
  }

  // Try the center of the domain as a deterministic safe candidate.
  const double midpoint = 0.5 * (lower_bound + upper_bound);
  for (double &coord : candidate) {
    coord = midpoint;
  }

  const double midpoint_fitness = sanitize_fitness(f.value(candidate));
  if (is_finite_value(midpoint_fitness)) {
    fitness_out = midpoint_fitness;
    return true;
  }

  return false;
}

} // namespace

CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   NetworkType topology, int shuffle_freq,
                   int stagnation_patience, double w_start,
                   double w_end, double coeff1, double coeff2,
                   unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm,
               SubSwarmTopologyConfig::from_type(topology), shuffle_freq,
               stagnation_patience, w_start, w_end, coeff1, coeff2, seed) {}

CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   const SubSwarmTopologyConfig &topology_config,
                   int shuffle_freq, int stagnation_patience, double w_start,
                   double w_end, double coeff1, double coeff2,
                   unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm,
               std::vector<SubSwarmTopologyConfig>(k_subswarms,
                                                   topology_config),
               shuffle_freq, stagnation_patience, w_start, w_end, coeff1,
               coeff2, seed) {}

CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   const std::vector<NetworkType> &topologies,
                   int shuffle_freq, int stagnation_patience,
                   double w_start, double w_end, double coeff1, double coeff2,
                   unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm,
               to_topology_configs(topologies), shuffle_freq,
               stagnation_patience, w_start, w_end, coeff1, coeff2, seed) {}

CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   const std::vector<SubSwarmTopologyConfig> &topologies,
                   int shuffle_freq, int stagnation_patience,
                   double w_start, double w_end, double coeff1, double coeff2,
                   unsigned int seed)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      shuffle_freq(shuffle_freq), stagnation_patience(stagnation_patience),
      w_max(w_start), w_min(w_end), c1(coeff1), c2(coeff2),
      master_seed(seed) {
  validate_cpso_parameters(k_subswarms, num_particles_per_swarm, shuffle_freq,
                           stagnation_patience);

  subswarm_topologies = normalize_topology_configs(
      k_subswarms, topologies, num_particles_per_swarm);
}

CpsoRunArtifacts CPSOBase::optimize_raw(
    const TestFunction &f, StoppingCriteriaManager &stop_manager) {
  if (num_subswarms > static_cast<int>(f.dim)) {
    throw std::invalid_argument(
        "number of subswarms must be <= function dimension");
  }

  unsigned int effective_seed = get_master_seed();
  int mpi_rank = 0;
  int mpi_size = 1;

#if CPSO_HAVE_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    MPI_Bcast(&effective_seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  }
#endif

  // Creating two different generators for velocity and topology (independents).
  std::vector<std::mt19937> gens(num_subswarms);
  std::vector<std::mt19937> topology_gens(num_subswarms);
  for (int i = 0; i < num_subswarms; ++i) {
    gens[i] = std::mt19937(effective_seed + 1u + static_cast<unsigned int>(i));
    topology_gens[i] =
        std::mt19937(effective_seed + 10001u + static_cast<unsigned int>(i));
  }
  std::mt19937 global_gen(effective_seed);

  int total_dim = f.dim;
  // Split the global dimension set into equal cooperative parts.
  auto bounds = f.get_domain();
  if (bounds.first > bounds.second) {
    std::swap(bounds.first, bounds.second);
  }
  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  // Prebuild the dimension split once. In MPI mode only the owned range is
  // advanced, but keeping the full list makes the batch indexing simple.
  std::vector<SubSwarm> swarms;
  swarms.reserve(num_subswarms);

  int current_dim_start = 0;
  for (int i = 0; i < num_subswarms; ++i) {
    int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);
    std::vector<int> active_dims;
    for (int d = 0; d < swarm_dims; ++d) {
      active_dims.push_back(current_dim_start + d);
    }
    current_dim_start += swarm_dims;

    std::vector<std::vector<int>> sub_adj_list = build_subswarm_topology(
        subswarm_topologies[i], particles_per_swarm, topology_gens[i]);

    swarms.emplace_back(particles_per_swarm, active_dims, bounds.first,
                        bounds.second, sub_adj_list);
  }

  ContextVector context(total_dim);

  // Start from a finite global context so every cooperative evaluation has a valid base.
  std::uniform_real_distribution<double> dist_init(bounds.first, bounds.second);
  std::vector<double> init_vec(total_dim);
  for (int i = 0; i < total_dim; ++i) {
    init_vec[i] = dist_init(global_gen);
  }

  double init_fitness = std::numeric_limits<double>::infinity();
  if (!try_build_finite_initial_fitness(f, init_vec, global_gen, bounds.first,
                                        bounds.second, init_fitness)) {
    CpsoRunArtifacts artifacts;
    artifacts.best_position = init_vec;
    artifacts.best_fitness = std::numeric_limits<double>::infinity();
    artifacts.seed_used = effective_seed;
    artifacts.cores = mpi_size;
    artifacts.iterations = 0;
    artifacts.stop_reason = CpsoStopReason::NUMERIC_FAILURE;
    artifacts.failure_message =
        "unable to build a finite initial CPSO context vector";
    return artifacts;
  }

  context.set_full_vector(init_vec, init_fitness);
  CpsoRunArtifacts artifacts =
      run_optimization_loop(f, stop_manager, swarms, gens, context, global_gen);
  artifacts.seed_used = effective_seed;
  return artifacts;
}

OutputObject CPSOBase::optimize(const TestFunction &f,
                                StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();
  CpsoRunArtifacts artifacts = optimize_raw(f, stop_manager);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;

  OutputObject out = build_cpso_output(f, artifacts, stop_manager);
  out.n_points = static_cast<unsigned int>(get_particles_per_swarm() *
                                           get_num_subswarms());
  out.execution_time = elapsed.count();
  return out;
}
