#include "CPSOParallel.hpp"
#include "../topology/create_network.hpp"
#include <algorithm>
#include <chrono>
#include <mpi.h>
#include <numeric>

CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           NetworkType topology, int shuffle_freq,
                           int stagnation_patience, double w_start,
                           double w_end, double c1, double c2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      shuffle_freq(shuffle_freq), stagnation_patience(stagnation_patience),
      w_max(w_start), w_min(w_end), c1(c1), c2(c2) {
  subswarm_topologies.assign(k_subswarms, topology);
}

CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           const std::vector<NetworkType> &topologies,
                           int shuffle_freq, int stagnation_patience,
                           double w_start, double w_end, double c1, double c2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      shuffle_freq(shuffle_freq), stagnation_patience(stagnation_patience),
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

  int stagnation_counter = 0;
  double previous_best_fitness = context.get_best_fitness();

  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;

  while (!must_stop) {
    iter++;
    stop_manager.increment_iterations();

    if (iter > 1 && iter % shuffle_freq == 0) {
      std::vector<int> permutation(total_dim);
      if (mpi_rank == 0) {
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), global_gen);
      }
#ifdef MPI_VERSION
      MPI_Bcast(permutation.data(), total_dim, MPI_INT, 0, MPI_COMM_WORLD);
#endif

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

    std::vector<std::vector<double>> subswarm_bests(num_subswarms);
    std::vector<double> subswarm_best_vals(num_subswarms);
    double progress_ratio =
        (double)stop_manager.get_current_iters() / stop_manager.get_max_iters();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    double current_w = w_max - (w_max - w_min) * progress_ratio;

    for (int i = local_start_idx; i < local_end_idx; ++i) {

      swarms[i].update_velocities_and_positions(current_w, c1, c2, gens[i]);

      swarms[i].evaluate_and_update(context, f);

      subswarm_bests[i] = swarms[i].get_gbest_pos();
      subswarm_best_vals[i] = swarms[i].get_gbest_val();
    }

#ifdef MPI_VERSION
    std::vector<double> send_context = context.get_full_vector();
    std::vector<double> recv_context(total_dim);

    int right_nbr = (mpi_rank + 1) % mpi_size;
    int left_nbr = (mpi_rank - 1 + mpi_size) % mpi_size;

    MPI_Sendrecv(send_context.data(), total_dim, MPI_DOUBLE, right_nbr, 0,
                 recv_context.data(), total_dim, MPI_DOUBLE, left_nbr, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = local_start_idx; i < local_end_idx; ++i) {
      const auto &active_dims = swarms[i].get_active_dims();
      for (size_t d = 0; d < active_dims.size(); ++d) {
        recv_context[active_dims[d]] = subswarm_bests[i][d];
      }
    }

    double new_true_fitness = f.value(recv_context);

    if (new_true_fitness < context.get_best_fitness()) {
      context.set_full_vector(recv_context, new_true_fitness);
    }
#else
    for (int i = 0; i < num_subswarms; ++i) {
      if (subswarm_best_vals[i] < context.get_best_fitness()) {
        context.update(subswarm_bests[i], swarms[i].get_active_dims(),
                       subswarm_best_vals[i]);
      }
    }

    double new_true_fitness = f.value(context.get_full_vector());

    if (new_true_fitness < context.get_best_fitness()) {
      context.set_full_vector(context.get_full_vector(), new_true_fitness);
    }
#endif

    double current_best_fitness = context.get_best_fitness();

    if (previous_best_fitness - current_best_fitness < 1e-6) {
      stagnation_counter++;
    } else {
      stagnation_counter = 0;
      previous_best_fitness = current_best_fitness;
    }

    if (stagnation_counter >= stagnation_patience) {
      for (int i = local_start_idx; i < local_end_idx; ++i) {
        swarms[i].inject_velocities(gens[i]);
      }
      stagnation_counter = 0;
    }

    const std::vector<double> &current_gbest_pos = context.get_full_vector();

    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // Calculate explicitly the average distance
    std::vector<double> local_dist_sq(particles_per_swarm, 0.0);
    for (int p = 0; p < particles_per_swarm; ++p) {
      for (int i = local_start_idx; i < local_end_idx; ++i) {
        const auto &particle = swarms[i].get_particles()[p];
        const auto &active_dims = swarms[i].get_active_dims();
        for (size_t d = 0; d < active_dims.size(); ++d) {
          double diff =
              particle.position[d] - current_gbest_pos[active_dims[d]];
          local_dist_sq[p] += diff * diff;
        }
      }
    }

    std::vector<double> global_dist_sq(particles_per_swarm, 0.0);
#ifdef MPI_VERSION
    MPI_Allreduce(local_dist_sq.data(), global_dist_sq.data(),
                  particles_per_swarm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    global_dist_sq = local_dist_sq;
#endif

    double total_distance = 0.0;
    for (int p = 0; p < particles_per_swarm; ++p) {
      total_distance += std::sqrt(global_dist_sq[p]);
    }
    double avg_distance = total_distance / particles_per_swarm;

    bool local_stop =
        stop_manager.should_stop(current_best_fitness, avg_distance);

#ifdef MPI_VERSION
    int local_stop_int = local_stop ? 1 : 0;
    int global_stop_int = 0;
    MPI_Allreduce(&local_stop_int, &global_stop_int, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    must_stop = (global_stop_int > 0);
#else
    must_stop = local_stop;
#endif
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
