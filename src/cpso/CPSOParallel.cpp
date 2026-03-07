/**
 * @file CPSOParallel.cpp
 * @brief Implementation of the CPSOParallel class methods inheriting from CPSOBase.
 */
#include "CPSOParallel.hpp"
#include <algorithm>
#include <numeric>

#if __has_include(<mpi.h>)
#include <mpi.h>
#endif

// Constructor for uniform topology
CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           NetworkType topology, int shuffle_freq,
                           int stagnation_patience, double w_start,
                           double w_end, double c1, double c2)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topology, 
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2) {}


// Constructor for different types of topologies
CPSOParallel::CPSOParallel(int k_subswarms, int num_particles_per_swarm,
                           const std::vector<NetworkType> &topologies,
                           int shuffle_freq, int stagnation_patience,
                           double w_start, double w_end, double c1, double c2)
    : CPSOBase(k_subswarms, num_particles_per_swarm, topologies, 
               shuffle_freq, stagnation_patience, w_start, w_end, c1, c2) {}


OutputObject CPSOParallel::run_optimization_loop(const TestFunction &f, StoppingCriteriaManager &stop_manager,
                                                 std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                                 ContextVector& context, std::mt19937& global_gen) {
  int mpi_rank = 0, mpi_size = 1;
#ifdef MPI_VERSION
  // Initialize variables for MPI rank and size
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  int total_dim = f.dim;
  
  // Distribute the sub-swarms across the available MPI processes
  std::vector<int> swarms_per_proc(mpi_size, num_subswarms / mpi_size);
  for (int i = 0; i < num_subswarms % mpi_size; ++i) {
    swarms_per_proc[i]++;
  }

  // Determine the start and end indices of the sub-swarms for the current MPI process
  int local_start_idx = 0;
  for (int i = 0; i < mpi_rank; ++i) {
    local_start_idx += swarms_per_proc[i];
  }
  int local_end_idx = local_start_idx + swarms_per_proc[mpi_rank];

  std::vector<double> init_vec = context.get_full_vector();
  double init_fitness = context.get_best_fitness();
  
  // Rank 0 initializes the vector and fitness
  if (mpi_rank == 0) {
      std::uniform_real_distribution<double> dist_pos(f.get_domain().first, f.get_domain().second);
      for (int i = 0; i < total_dim; ++i)
          init_vec[i] = dist_pos(global_gen);
      init_fitness = f.value(init_vec);
  }

  // Broadcast the initial vector and its fitness to all processes
#ifdef MPI_VERSION
  MPI_Bcast(init_vec.data(), total_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&init_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  context.set_full_vector(init_vec, init_fitness);

  // Initialize the particles in each sub-swarm
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
  double last_avg_distance = std::numeric_limits<double>::infinity();

  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  // Main optimization loop
  while (!must_stop) {
    iter++;
    stop_manager.increment_iterations();

    // Periodically shuffle the dimensions among the sub-swarms
    if (iter > 1 && iter % shuffle_freq == 0) {
#ifdef MPI_VERSION
      // Synchronize all processes before shuffling
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      // Rank 0 creates a permutation of the dimensions
      std::vector<int> permutation(total_dim);
      if (mpi_rank == 0) {
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), global_gen);
      }
#ifdef MPI_VERSION
      // Send the permutation to all the MPI ranks
      MPI_Bcast(permutation.data(), total_dim, MPI_INT, 0, MPI_COMM_WORLD);
#endif

      // Update the active dimensions for each sub-swarm
      int current_dim_start = 0;
      for (int i = 0; i < num_subswarms; ++i) {
        int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);
        std::vector<int> new_active_dims(swarm_dims);

        // Fill the new active dimensions
        for (int d = 0; d < swarm_dims; ++d) {
          new_active_dims[d] = permutation[current_dim_start + d];
        }
        current_dim_start += swarm_dims;

        // Update the active dimensions for the current sub-swarm
        swarms[i].update_active_dims(new_active_dims, context, gens[i]);
      }
    }

    std::vector<std::vector<double>> subswarm_bests(num_subswarms);
    std::vector<double> subswarm_best_vals(num_subswarms);
    double progress_ratio =
        (double)stop_manager.get_current_iters() / stop_manager.get_max_iters();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    // Calculate the current inertia weight
    double current_w = w_max - (w_max - w_min) * progress_ratio;

    // Update velocities, positions, and evaluate fitness for each local sub-swarm
    for (int i = local_start_idx; i < local_end_idx; ++i) {

      swarms[i].update_velocities_and_positions(current_w, c1, c2, gens[i]);

      swarms[i].evaluate_and_update(context, f);

      subswarm_bests[i] = swarms[i].get_gbest_pos();
      subswarm_best_vals[i] = swarms[i].get_gbest_val();
    }

    // Synchronize the context vector across processes
#ifdef MPI_VERSION
    std::vector<double> send_context = context.get_full_vector();
    std::vector<double> recv_context(total_dim);

    // Using a ring topology to send and receive the context vector
    int right_nbr = (mpi_rank + 1) % mpi_size;
    int left_nbr = (mpi_rank - 1 + mpi_size) % mpi_size;

    MPI_Sendrecv(send_context.data(), total_dim, MPI_DOUBLE, right_nbr, 0,
                 recv_context.data(), total_dim, MPI_DOUBLE, left_nbr, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Update the received context vector with the best local positions
    for (int i = local_start_idx; i < local_end_idx; ++i) {
      const auto &active_dims = swarms[i].get_active_dims();
      for (size_t d = 0; d < active_dims.size(); ++d) {
        recv_context[active_dims[d]] = subswarm_bests[i][d];
      }
    }

    // Evaluate the new full context vector
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

    // Check for stagnation
    if (previous_best_fitness - current_best_fitness < 1e-6) {
      stagnation_counter++;
    } else {
      stagnation_counter = 0;
      previous_best_fitness = current_best_fitness;
    }

    // Inject random velocities if stagnation patience is exceeded
    if (stagnation_counter >= stagnation_patience) {
      for (int i = local_start_idx; i < local_end_idx; ++i) {
        swarms[i].inject_velocities(gens[i]);
      }
      stagnation_counter = 0;
    }

    const std::vector<double> &current_gbest_pos = context.get_full_vector();

    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // Compute average distance between particles and global best position
    compute_avg_distance(iter, swarms, current_gbest_pos, last_avg_distance, local_start_idx, local_end_idx, true);

    bool local_stop = stop_manager.should_stop(current_best_fitness, last_avg_distance);

    // Synchronize the stop condition across all processes
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

  OutputObject out(f.get_name(), total_dim, particles_per_swarm * num_subswarms,
                   context.get_full_vector(), f.get_true_solution(),
                   context.get_best_fitness(), history, 1, 0.0,
                   iter, stop_manager);

  return out;
}
