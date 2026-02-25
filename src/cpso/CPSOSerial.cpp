#include "CPSOSerial.hpp"
#include "../topology/create_network.hpp"
#include <chrono>

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       NetworkType topology, double w_start, double w_end,
                       double c1, double c2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      w_max(w_start), w_min(w_end), c1(c1), c2(c2) {
  subswarm_topologies.assign(k_subswarms, topology);
}

CPSOSerial::CPSOSerial(int k_subswarms, int num_particles_per_swarm,
                       const std::vector<NetworkType> &topologies,
                       double w_start, double w_end, double c1, double c2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      w_max(w_start), w_min(w_end), c1(c1), c2(c2) {
  subswarm_topologies = topologies;
  while (subswarm_topologies.size() < static_cast<size_t>(k_subswarms)) {
    subswarm_topologies.push_back(topologies.empty() ? NetworkType::SCALE_FREE
                                                     : topologies.front());
  }
}

OutputObject CPSOSerial::optimize(const TestFunction &f,
                                  StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();
  std::mt19937 gen(42);

  int total_dim = f.dim;
  auto bounds = f.get_domain();
  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

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

    std::vector<std::vector<int>> sub_adj_list;
    switch (subswarm_topologies[i]) {
    case NetworkType::SMALL_WORLD:
      create_network(particles_per_swarm, 0.3, sub_adj_list);
      break;
    case NetworkType::SCALE_FREE:
      create_scale_free_network(particles_per_swarm, 2, sub_adj_list);
      break;
    case NetworkType::RANDOM:
      create_random_network(particles_per_swarm, 0.5, sub_adj_list);
      break;
    case NetworkType::FULLY_CONNECTED:
      create_fully_connected_network(particles_per_swarm, sub_adj_list);
      break;
    }

    swarms.emplace_back(particles_per_swarm, active_dims, bounds.first,
                        bounds.second, sub_adj_list);
  }

  ContextVector context(total_dim);
  context.set_full_vector(
      std::vector<double>(total_dim,
                          bounds.first + (bounds.second - bounds.first) / 2.0),
      std::numeric_limits<double>::infinity());

  std::vector<double> init_vec(total_dim);
  std::uniform_real_distribution<double> dist_pos(bounds.first, bounds.second);
  for (int i = 0; i < total_dim; ++i)
    init_vec[i] = dist_pos(gen);
  context.set_full_vector(init_vec, f.value(init_vec));

  for (auto &swarm : swarms) {
    swarm.initialize(gen, context, f);
    context.update(swarm.get_gbest_pos(), swarm.get_active_dims(),
                   swarm.get_gbest_val());
  }
  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;

  while (!must_stop) {
    iter++;

    double progress_ratio = (double)stop_manager.get_current_fevals() /
                            stop_manager.get_max_fevals();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    double current_w = w_max - (w_max - w_min) * progress_ratio;
    for (auto &swarm : swarms) {
      swarm.update_velocities_and_positions(current_w, c1, c2, gen);
      int fevals = swarm.evaluate_and_update(context, f);
      stop_manager.add_evaluations(fevals);
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
