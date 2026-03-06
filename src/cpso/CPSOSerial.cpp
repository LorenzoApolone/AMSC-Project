#include "CPSOSerial.hpp"
#include <algorithm>
#include <numeric>

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       NetworkType topology, int shuffle_freq,
                       int stagnation_patience, double w_start, double w_end,
                       double c1, double c2)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology, shuffle_freq,
               stagnation_patience, w_start, w_end, c1, c2) {}

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

  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  int stagnation_counter = 0;
  double previous_best_fitness = context.get_best_fitness();

  // Initialize the sub-swarms mathematically evaluating the context
  for (int i = 0; i < num_subswarms; ++i) {
    swarms[i].initialize(gens[i], context, f);
    context.update(swarms[i].get_gbest_pos(), swarms[i].get_active_dims(),
                   swarms[i].get_gbest_val());
  }

  while (!must_stop) {
    iter++;
    stop_manager.increment_iterations();

    // MISSING BUG FIX 1: Add Stagnation & Inject Velocities (from CPSOBase traits)
    // Periodically shuffle the dimensions among the sub-swarms
    if (iter > 1 && iter % shuffle_freq == 0) {
      std::vector<int> permutation(total_dim);
      std::iota(permutation.begin(), permutation.end(), 0);
      std::shuffle(permutation.begin(), permutation.end(), global_gen);

      int current_dim_start = 0;
      for (int i = 0; i < num_subswarms; ++i) {
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

    double current_w = w_max - (w_max - w_min) * progress_ratio;
    for (int i = 0; i < num_subswarms; ++i) {
      swarms[i].update_velocities_and_positions(current_w, c1, c2, gens[i]);
      swarms[i].evaluate_and_update(context, f);
    }

    for (auto &swarm : swarms) {
      if (swarm.get_gbest_val() < context.get_best_fitness()) {
        context.update(swarm.get_gbest_pos(), swarm.get_active_dims(),
                       swarm.get_gbest_val());
      }
    }

    double new_true_fitness = f.value(context.get_full_vector());
    context.set_full_vector(
        context.get_full_vector(),
        std::min(context.get_best_fitness(), new_true_fitness));

    double current_best_fitness = context.get_best_fitness();
    const std::vector<double> &current_gbest_pos = context.get_full_vector();
    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // MISSING BUG FIX 2: Added Stagnation loop
    if (previous_best_fitness - current_best_fitness < 1e-6) {
      stagnation_counter++;
    } else {
      stagnation_counter = 0;
      previous_best_fitness = current_best_fitness;
    }

    if (stagnation_counter >= stagnation_patience) {
      for (int i = 0; i < num_subswarms; ++i) {
        swarms[i].inject_velocities(gens[i]);
      }
      stagnation_counter = 0;
    }

    // Evaluate the distance avoiding duplications using the CPSOBase template routine
    compute_avg_distance(iter, swarms, current_gbest_pos, last_avg_distance, 0, num_subswarms, false);

    if (stop_manager.should_stop(current_best_fitness, last_avg_distance)) {
      must_stop = true;
    }
  }

  // CPSOBase optimize() driver handles Elapsed time tracking
  OutputObject out(f.get_name(), total_dim, particles_per_swarm * num_subswarms,
                   context.get_full_vector(), f.get_true_solution(),
                   context.get_best_fitness(), history, 1, 0.0,
                   iter, stop_manager);

  return out;
}
