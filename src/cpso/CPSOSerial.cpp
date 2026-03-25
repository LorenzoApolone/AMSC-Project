/**
 * @file CPSOSerial.cpp
 * @brief Implementation of the CPSOSerial class methods
 */
#include "CPSOSerial.hpp"
#include "SwarmMetrics.hpp"
#include <algorithm>
#include <numeric>

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       NetworkType topology, int shuffle_freq,
                       int stagnation_patience, double w_start, double w_end,
                       double c1, double c2, unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2, seed) {}

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       const SubSwarmTopologyConfig &topology_config,
                       int shuffle_freq, int stagnation_patience,
                       double w_start, double w_end, double c1, double c2,
                       unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology_config,
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2,
               seed) {}

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       const std::vector<NetworkType> &topologies,
                       int shuffle_freq, int stagnation_patience,
                       double w_start, double w_end,
                       double c1, double c2, unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2, seed) {}

CPSOSerial::CPSOSerial(
    int k_subswarms, int num_particles_per_swarm,
    const std::vector<SubSwarmTopologyConfig> &topologies, int shuffle_freq,
    int stagnation_patience, double w_start, double w_end, double c1,
    double c2, unsigned int seed)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2, seed) {}

CpsoRunArtifacts CPSOSerial::run_optimization_loop(
    const TestFunction &f, StoppingCriteriaManager &stop_manager,
    std::vector<SubSwarm> &swarms, std::vector<std::mt19937> &gens,
    ContextVector &context, std::mt19937 &global_gen) {
  // Parameters kept during the exec
  std::vector<double> fitness_history;
  int iter = 0;

  // Lambda function used to take track of the failure of the algorithm
  auto build_numeric_failure = [&](const std::string &message) {
    CpsoRunArtifacts artifacts;
    artifacts.best_position = context.get_full_vector();
    artifacts.best_fitness = std::numeric_limits<double>::infinity();
    artifacts.best_fitness_history = fitness_history;
    artifacts.cores = 1;
    artifacts.iterations = iter;
    artifacts.stop_reason = CpsoStopReason::NUMERIC_FAILURE;
    artifacts.failure_message = message;
    return artifacts;
  };

  try {
    // Derive the composition of the funcions
    int total_dim = f.dim;
    bool must_stop = false;
    double last_avg_distance = std::numeric_limits<double>::infinity();
    CpsoStopReason stop_reason = CpsoStopReason::UNKNOWN;

    int dims_per_swarm = total_dim / get_num_subswarms();
    int remainder = total_dim % get_num_subswarms();
    int injection_count = 0;

    // Each sub-swarm initializes locally, and if its local gbest already improves the current shared fitness, its coordinates are copied into the global context immediately.
    std::vector<double> init_vector = context.get_full_vector();
    double init_fitness = context.get_best_fitness();

    for (int i = 0; i < get_num_subswarms(); ++i) {
      swarms[i].initialize(gens[i], context, f);
      if (swarms[i].get_gbest_val() < init_fitness) {
        const auto &active_dims = swarms[i].get_active_dims();
        const auto &best_pos = swarms[i].get_gbest_pos();
        for (size_t d = 0; d < active_dims.size(); ++d) {
          init_vector[active_dims[d]] = best_pos[d];
        }
        init_fitness = swarms[i].get_gbest_val();
        context.set_full_vector(init_vector, init_fitness);
      }
    }

    // Main optimization loop
    while (!must_stop) {
      iter++;
      stop_manager.increment_iterations();

      // Periodically reshuffle the dimension assignment of the sub-swarms.
      if (iter > 1 && iter % get_shuffle_freq() == 0) {
        std::vector<int> permutation(total_dim);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), global_gen);

        int current_dim_start = 0;
        for (int i = 0; i < get_num_subswarms(); ++i) {
          int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);
          std::vector<int> new_active_dims(swarm_dims);
          for (int d = 0; d < swarm_dims; ++d) {
            new_active_dims[d] = permutation[current_dim_start + d];
          }
          current_dim_start += swarm_dims;
          swarms[i].update_active_dims(new_active_dims, context, gens[i]);
        }
      }

      // Convert iteration progress into the inertia weight used by the PSO update rule for this step.
      double progress_ratio =
          static_cast<double>(stop_manager.get_current_iters()) /
          stop_manager.get_max_iters();
      if (progress_ratio > 1.0) {
        progress_ratio = 1.0;
      }

      double current_w =
          get_w_max() - (get_w_max() - get_w_min()) * progress_ratio;

      // Advance each sub-swarm and immediately fold any improving local gbest into the shared context.
      for (int i = 0; i < get_num_subswarms(); ++i) {
        // Evaluate the sub_swarm of the latest context and perform one PSO algorithm
        swarms[i].recalculate_fitness(context, f);
        swarms[i].update_velocities_and_positions(current_w, get_c1(), get_c2(),
                                                  gens[i], progress_ratio);
        swarms[i].evaluate_and_update(context, f);


        double current_fitness = context.get_best_fitness();
        if (swarms[i].get_gbest_val() < current_fitness) {
          double new_fitness = swarms[i].get_gbest_val();
          std::vector<double> current_vector = context.get_full_vector();
          const auto &active_dims = swarms[i].get_active_dims();
          const auto &best_pos = swarms[i].get_gbest_pos();
          for (size_t d = 0; d < active_dims.size(); ++d) {
            current_vector[active_dims[d]] = best_pos[d];
          }
          context.set_full_vector(current_vector, new_fitness);
        }
      }

      double current_best_fitness = context.get_best_fitness();
      const std::vector<double> &current_gbest_pos = context.get_full_vector();
      fitness_history.push_back(current_best_fitness);

      // The average distance is computed over all sub-swarms and given to the SC
      SwarmMetrics::compute_avg_distance(swarms, current_gbest_pos,
                                         last_avg_distance, 0,
                                         get_num_subswarms(), false);

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

      // Applying velocity injection toi prevent stagnation in a false positive position
      if (should_inject) {
        injection_count++;
        bool hard_reset = (injection_count % 3 == 0);
        for (int i = 0; i < get_num_subswarms(); ++i) {
          swarms[i].inject_velocities(gens[i], hard_reset);
          if (hard_reset) {
            swarms[i].reset_gbest_attraction();
          }
        }

        stop_manager.reset_stagnation();
        if (stop_for_stagnation) {
          local_stop = false;
        }
      }

      if (local_stop) {
        stop_reason = infer_cpso_stop_reason(stop_for_max_iters,
                                             stop_for_low_diversity,
                                             stop_for_stagnation);
        must_stop = true;
      }
    }

    
    CpsoRunArtifacts artifacts;
    artifacts.best_position = context.get_full_vector();
    artifacts.best_fitness = context.get_best_fitness();
    artifacts.best_fitness_history = fitness_history;
    artifacts.cores = 1;
    artifacts.iterations = iter;
    artifacts.stop_reason = stop_reason;
    return artifacts;
  } catch (const std::exception &ex) {
    return build_numeric_failure(std::string("serial CPSO numeric failure: ") + ex.what());
  }
}