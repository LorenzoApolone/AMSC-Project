#include "CPSOParallel.hpp"
#include "../topology/create_network.hpp"
#include <algorithm>
#include <chrono>
#include <mpi.h>
#include <numeric>

CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           NetworkType topology, double w_start, double w_end,
                           double c1, double c2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      w_max(w_start), w_min(w_end), c1(c1), c2(c2) {
  subswarm_topologies.assign(k_subswarms, topology);
}

CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
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

OutputObject CPSOParallel::optimize(const TestFunction &f,
                                    StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();

  int mpi_rank = 0, mpi_size = 1;
#ifdef MPI_VERSION
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  std::vector<std::mt19937> gens(num_subswarms);
  for (int i = 0; i < num_subswarms; ++i) {
    gens[i] = std::mt19937(1337 + i);
  }

  std::mt19937 global_gen(1337);
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

  std::vector<int> swarms_per_proc(mpi_size, num_subswarms / mpi_size);
  for (int i = 0; i < num_subswarms % mpi_size; ++i) {
    swarms_per_proc[i]++;
  }

  int local_start_idx = 0;
  for (int i = 0; i < mpi_rank; ++i) {
    local_start_idx += swarms_per_proc[i];
  }
  int local_end_idx = local_start_idx + swarms_per_proc[mpi_rank];

  ContextVector context(total_dim);
  std::vector<double> init_vec(total_dim);
  double init_fitness = 0.0;

  if (mpi_rank == 0) {
    std::uniform_real_distribution<double> dist_pos(bounds.first,
                                                    bounds.second);
    for (int i = 0; i < total_dim; ++i)
      init_vec[i] = dist_pos(global_gen);
    init_fitness = f.value(init_vec);
  }

#ifdef MPI_VERSION
  MPI_Bcast(init_vec.data(), total_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&init_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  context.set_full_vector(init_vec, init_fitness);

  for (int i = 0; i < num_subswarms; ++i) {
    swarms[i].initialize(gens[i], context, f);
    context.update(swarms[i].get_gbest_pos(), swarms[i].get_active_dims(),
                   swarms[i].get_gbest_val());
  }

  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;

  while (!must_stop) {
    iter++;
    std::vector<std::vector<double>> subswarm_bests(num_subswarms);
    std::vector<double> subswarm_best_vals(num_subswarms);
    double progress_ratio = (double)stop_manager.get_current_fevals() /
                            stop_manager.get_max_fevals();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    double current_w = w_max - (w_max - w_min) * progress_ratio;

    int local_fevals = 0;

    for (int i = local_start_idx; i < local_end_idx; ++i) {

      swarms[i].update_velocities_and_positions(current_w, c1, c2, gens[i]);

      int fevals = swarms[i].evaluate_and_update(context, f);
      local_fevals += fevals;

      subswarm_bests[i] = swarms[i].get_gbest_pos();
      subswarm_best_vals[i] = swarms[i].get_gbest_val();
    }

#ifdef MPI_VERSION
    int global_fevals = 0;
    MPI_Allreduce(&local_fevals, &global_fevals, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    stop_manager.add_evaluations(global_fevals);

    std::vector<double> local_vals;
    std::vector<double> local_pos;
    for (int i = local_start_idx; i < local_end_idx; ++i) {
      local_vals.push_back(subswarm_best_vals[i]);
      for (double p : subswarm_bests[i]) {
        local_pos.push_back(p);
      }
    }

    std::vector<int> recvcounts_vals(mpi_size);
    std::vector<int> displs_vals(mpi_size, 0);
    std::vector<int> recvcounts_pos(mpi_size);
    std::vector<int> displs_pos(mpi_size, 0);

    int total_active_dims = 0;
    for (int i = 0; i < mpi_size; ++i) {
      recvcounts_vals[i] = swarms_per_proc[i];

      int elements = 0;
      int p_start = 0;
      for (int j = 0; j < i; ++j)
        p_start += swarms_per_proc[j];
      int p_end = p_start + swarms_per_proc[i];

      for (int k = p_start; k < p_end; ++k) {
        elements += swarms[k].get_active_dims().size();
      }
      recvcounts_pos[i] = elements;
      total_active_dims += elements;

      if (i > 0) {
        displs_vals[i] = displs_vals[i - 1] + recvcounts_vals[i - 1];
        displs_pos[i] = displs_pos[i - 1] + recvcounts_pos[i - 1];
      }
    }

    std::vector<double> global_vals(num_subswarms);
    std::vector<double> global_pos(total_active_dims);

    MPI_Allgatherv(local_vals.data(), local_vals.size(), MPI_DOUBLE,
                   global_vals.data(), recvcounts_vals.data(),
                   displs_vals.data(), MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Allgatherv(local_pos.data(), local_pos.size(), MPI_DOUBLE,
                   global_pos.data(), recvcounts_pos.data(), displs_pos.data(),
                   MPI_DOUBLE, MPI_COMM_WORLD);

    int pos_idx = 0;
    for (int i = 0; i < num_subswarms; ++i) {
      subswarm_best_vals[i] = global_vals[i];

      int dims = swarms[i].get_active_dims().size();
      subswarm_bests[i].resize(dims);
      for (int d = 0; d < dims; ++d) {
        subswarm_bests[i][d] = global_pos[pos_idx++];
      }
    }

#else
    stop_manager.add_evaluations(local_fevals);
#endif

    for (int i = 0; i < num_subswarms; ++i) {
      if (subswarm_best_vals[i] < context.get_best_fitness()) {
        context.update(subswarm_bests[i], swarms[i].get_active_dims(),
                       subswarm_best_vals[i]);
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
