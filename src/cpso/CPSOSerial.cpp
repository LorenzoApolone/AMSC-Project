#include "CPSOSerial.hpp"
#include "../topology/create_network.hpp"
#include <chrono>

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       double w_start, double w_end, double coeff1,
                       double coeff2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      w_max(w_start), w_min(w_end), c1(coeff1), c2(coeff2) {}

OutputObject CPSOSerial::optimize(const TestFunction &f,
                                  StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();
  std::mt19937 gen(42); // Random Generation of the particles

  int total_dim = f.dim;
  auto bounds = f.get_domain();

  // Dimensions Division and Sub-Swarms Initialization
  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  std::vector<SubSwarm> swarms;
  swarms.reserve(num_subswarms);

  int current_dim_start = 0;
  for (int i = 0; i < num_subswarms; ++i) {

    // Distribute remaining dimensions among sub-swarms
    int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);

    std::vector<int> active_dims;
    for (int d = 0; d < swarm_dims; ++d) {
      active_dims.push_back(current_dim_start + d);
    }
    current_dim_start += swarm_dims;

    // Scale-Free Network Generation for this Sub-Swarm
    std::vector<std::vector<int>> sub_adj_list;
    create_scale_free_network(particles_per_swarm, 2, sub_adj_list);

    swarms.emplace_back(particles_per_swarm, active_dims, bounds.first,
                        bounds.second, sub_adj_list);
  }

  // Context Vector Initialization
  ContextVector context(total_dim);
  context.set_full_vector(
      std::vector<double>(total_dim,
                          bounds.first + (bounds.second - bounds.first) / 2.0),
      std::numeric_limits<double>::infinity()); // Init in the middle

  // Initial vector evaluation (creating a fictitious random best)
  std::vector<double> init_vec(total_dim);
  std::uniform_real_distribution<double> dist_pos(bounds.first, bounds.second);
  for (int i = 0; i < total_dim; ++i)
    init_vec[i] = dist_pos(gen);
  context.set_full_vector(init_vec, f.value(init_vec));

  // Particles Initialization in Sub-Swarms
  for (auto &swarm : swarms) {
    swarm.initialize(gen, context, f);
    // IMMEDIATE context update with the best init of this swarm
    context.update(swarm.get_gbest_pos(), swarm.get_active_dims(),
                   swarm.get_gbest_val());
  }

  // History structure
  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;

  // Main Optimization Loop
  while (!must_stop) {
    iter++;

    // Linearly decreasing inertial weight dependent on FEvals
    double progress_ratio = (double)stop_manager.get_current_fevals() /
                            stop_manager.get_max_fevals();

    // Ensure ratio doesn't exceed 1.0 (failsafe)
    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    double current_w = w_max - (w_max - w_min) * progress_ratio;

    for (auto &swarm : swarms) {
      // Updates positions and velocities of the sub-swarm
      swarm.update_velocities_and_positions(current_w, c1, c2, gen);

      // Evaluates exploiting the Context Vector
      int fevals = swarm.evaluate_and_update(context, f);
      stop_manager.add_evaluations(fevals);

      // CPSO-S (Serial): We update the Context Vector after the
      // sub-swarm's turn
      context.update(swarm.get_gbest_pos(), swarm.get_active_dims(),
                     swarm.get_gbest_val());
    }

    // Stop Criteria
    double current_best_fitness = context.get_best_fitness();
    const std::vector<double> &current_gbest_pos = context.get_full_vector();

    // Save the current error
    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // Extracting positions for the manager (Swarm Diversity)
    std::vector<std::vector<double>> empty_diversity;

    if (stop_manager.should_stop(current_best_fitness, empty_diversity,
                                 current_gbest_pos)) {
      must_stop = true;
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;

  StopCriterion dummy_stop(iter, 0.0);

  OutputObject out(f.get_name(), total_dim, particles_per_swarm * num_subswarms,
                   context.get_full_vector(), f.get_true_solution(),
                   context.get_best_fitness(), history, 1, elapsed.count(),
                   iter, dummy_stop);

  return out;
}
