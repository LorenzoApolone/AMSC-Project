/**
 * @file CPSOSerial.cpp
 * @brief Implementation of the CPSOSerial class methods
 */
#include "CPSOSerial.hpp"
#include "SwarmMetrics.hpp"
#include <algorithm>
#include <numeric>

// Constructor of the CPSOSerial class with uniform network topology for all sub-swarms
CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       NetworkType topology, int shuffle_freq,
                       int stagnation_patience, double w_start, double w_end,
                       double c1, double c2)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2) {}

// Constructor of the CPSOSerial class with specific network topologies for each sub-swarm
CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       const std::vector<NetworkType> &topologies,
                       int shuffle_freq, int stagnation_patience, double w_start, double w_end,
                       double c1, double c2)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2) {}


OutputObject CPSOSerial::run_optimization_loop(
    const TestFunction &f, StoppingCriteriaManager &stop_manager,
    std::vector<SubSwarm> &swarms, std::vector<std::mt19937> &gens,
    ContextVector &context, std::mt19937 &global_gen) {

  int total_dim = f.dim;
  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;
  double last_avg_distance = std::numeric_limits<double>::infinity();

  int dims_per_swarm = total_dim / get_num_subswarms();
  int remainder = total_dim % get_num_subswarms();

  int injection_count = 0;

  std::vector<double> init_vector = context.get_full_vector();
  double init_fitness = context.get_best_fitness();

  // Initialize the sub-swarms within the particles, propagating context progressively
  for (int i = 0; i < get_num_subswarms(); ++i) {
    swarms[i].initialize(gens[i], context, f);
    if (swarms[i].get_gbest_val() < init_fitness) {
        const auto &active_dims = swarms[i].get_active_dims();
        const auto &best_pos = swarms[i].get_gbest_pos();
        for (size_t d = 0; d < active_dims.size(); ++d){
            init_vector[active_dims[d]] = best_pos[d];
        }
        init_fitness = swarms[i].get_gbest_val();
        context.set_full_vector(init_vector, init_fitness); 
    }
  }

  while (!must_stop) {
    iter++;
    stop_manager.increment_iterations();

    // Periodically shuffle the dimensions among the sub-swarms
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

    double progress_ratio =
        (double)stop_manager.get_current_iters() / stop_manager.get_max_iters();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    // Update the inertia weight
    double current_w = get_w_max() - (get_w_max() - get_w_min()) * progress_ratio;
    
    // Sequential update: iterate over each sub-swarm and apply updates greedily
    for (int i = 0; i < get_num_subswarms(); ++i) {
      // Recalculate historical fitness using the *current* context
      swarms[i].recalculate_fitness(context, f);

      // Update velocities and positions
      swarms[i].update_velocities_and_positions(current_w, get_c1(), get_c2(), gens[i], progress_ratio);

      // Evaluate and update local/global bests for the subswarm
      swarms[i].evaluate_and_update(context, f);

      // Immediately merge if it improves the global context
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
    
    // Update the history
    double current_best_fitness = context.get_best_fitness();
    const std::vector<double> &current_gbest_pos = context.get_full_vector();
    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // Evaluate the distance
    SwarmMetrics::compute_avg_distance(iter, swarms, current_gbest_pos, last_avg_distance, 0, get_num_subswarms(), false);

    bool local_stop = stop_manager.should_stop(current_best_fitness, last_avg_distance);

    // Inject random velocities if stagnation patience is exceeded
    int current_stag_iters = stop_manager.get_current_stagnation_iters();
    if (current_stag_iters > 0 && (current_stag_iters % get_stagnation_patience() == 0)) {
      injection_count++;
      bool hard_reset = (injection_count % 3 == 0);
      for (int i = 0; i < get_num_subswarms(); ++i) {
        swarms[i].inject_velocities(gens[i], hard_reset);
        if (hard_reset) {
            swarms[i].reset_gbest_attraction();
        }
      }
    }

    if (local_stop) {
      must_stop = true;
    }
  }

  OutputObject out(f.get_name(), total_dim, get_particles_per_swarm() * get_num_subswarms(),
                   context.get_full_vector(), f.get_true_solution(),
                   context.get_best_fitness(), history, 1, 0.0,
                   iter, stop_manager);

  return out;
}
