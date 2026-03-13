/**
 * @file CPSOParallel.cpp
 * @brief Implementation of the CPSOParallel class methods inheriting from CPSOBase
 */
#include "CPSOParallel.hpp"
#include "SwarmMetrics.hpp"
#include <algorithm>
#include <numeric>

#if __has_include(<mpi.h>)
#include <mpi.h>
#define CPSO_HAVE_MPI 1
#else
#define CPSO_HAVE_MPI 0
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
#if CPSO_HAVE_MPI
  // Initialize variables for MPI rank and size
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  int total_dim = f.dim;

  int effective_mpi_size = std::min(mpi_size, get_num_subswarms());

  // Distribute the sub-swarms across the available MPI processes
  std::vector<int> swarms_per_proc(mpi_size, 0);
  int base_count = get_num_subswarms() / effective_mpi_size;
  int remainder_swarms = get_num_subswarms() % effective_mpi_size;
  
  for (int i = 0; i < effective_mpi_size; ++i) {
    swarms_per_proc[i] = base_count + (i < remainder_swarms ? 1 : 0);
  }

  // Determine the start and end indices of the sub-swarms for the current MPI process
  int local_start_idx = 0;
  for (int i = 0; i < mpi_rank; ++i) {
    local_start_idx += swarms_per_proc[i];
  }
  int local_end_idx = local_start_idx + swarms_per_proc[mpi_rank];

  std::vector<double> init_vec = context.get_full_vector();
  double init_fitness = context.get_best_fitness();

  // Broadcast the initial vector and its fitness to all processes
#if CPSO_HAVE_MPI
  MPI_Bcast(init_vec.data(), total_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&init_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  context.set_full_vector(init_vec, init_fitness);

  int local_num_swarms = local_end_idx - local_start_idx;
  int global_max_local_swarms = 0;
#if CPSO_HAVE_MPI
  MPI_Allreduce(&local_num_swarms, &global_max_local_swarms, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
  global_max_local_swarms = local_num_swarms;
#endif

  // Initialization in micro-batches
  for (int b = 0; b < global_max_local_swarms; ++b) {
      std::vector<double> init_delta(total_dim, 0.0);
      double current_base_fitness = context.get_best_fitness();
      std::vector<double> current_base_vec = context.get_full_vector();
      int i = local_start_idx + b;
      
      if (i < local_end_idx) {
          swarms[i].initialize(gens[i], context, f);
          if (swarms[i].get_gbest_val() < current_base_fitness) {
              const auto &active_dims = swarms[i].get_active_dims();
              for (size_t d = 0; d < active_dims.size(); ++d) {
                  init_delta[active_dims[d]] = swarms[i].get_gbest_pos()[d] - current_base_vec[active_dims[d]];
              }
          }
      }

#if CPSO_HAVE_MPI
      // Gather all local improvements from different MPI ranks
      std::vector<double> init_global_delta(total_dim, 0.0);
      MPI_Allreduce(init_delta.data(), init_global_delta.data(), total_dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // Apply all improvements at once
      std::vector<double> synced_init_vector = current_base_vec;
      for (int dim = 0; dim < total_dim; ++dim) {
          synced_init_vector[dim] += init_global_delta[dim];
      }

      double init_synced_best = 0.0;
      if (mpi_rank == 0) {
          init_synced_best = f.value(synced_init_vector);
      }
      MPI_Bcast(&init_synced_best, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (init_synced_best < current_base_fitness) {
          //The combination of all updates improves the solution
          context.set_full_vector(synced_init_vector, init_synced_best);
      } else {
          // The combined updates degrade the solution
          std::vector<double> current_greedy = current_base_vec;
          double current_greedy_fitness = current_base_fitness;

          std::vector<double> local_imps(get_num_subswarms(), 0.0);
          if (i < local_end_idx && swarms[i].get_gbest_val() < current_base_fitness) {
              local_imps[i] = current_base_fitness - swarms[i].get_gbest_val();
          }

          // Share improvements across all ranks
          std::vector<double> global_imps(get_num_subswarms(), 0.0);
          MPI_Allreduce(local_imps.data(), global_imps.data(), get_num_subswarms(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

          // Sort subswarms by improvement
          std::vector<int> sorted_swarms(get_num_subswarms());
          std::iota(sorted_swarms.begin(), sorted_swarms.end(), 0);
          std::sort(sorted_swarms.begin(), sorted_swarms.end(), [&](int a, int b) {
              return global_imps[a] > global_imps[b];
          });

          // Apply updates one by one and keep only those that actually improve the solution
          for (int swarm_idx : sorted_swarms) {
              if (global_imps[swarm_idx] <= 0.0) continue; 

              std::vector<double> test_greedy = current_greedy;
              const auto& active_dims = swarms[swarm_idx].get_active_dims();
              for (int d : active_dims) {
                 test_greedy[d] += init_global_delta[d];
              }

              double test_fitness = 0.0;
              if (mpi_rank == 0) {
                  test_fitness = f.value(test_greedy);
              }
              MPI_Bcast(&test_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              
              if (test_fitness < current_greedy_fitness) {
                  current_greedy = test_greedy;
                  current_greedy_fitness = test_fitness;
              }
          }
          context.set_full_vector(current_greedy, current_greedy_fitness);
      }
#else
      if (i < local_end_idx) {
          std::vector<double> local_init_vector = current_base_vec;
          for (int dim = 0; dim < total_dim; ++dim) {
              local_init_vector[dim] += init_delta[dim];
          }
          double init_local_best = f.value(local_init_vector);
          
          if (init_local_best < current_base_fitness) {
              context.set_full_vector(local_init_vector, init_local_best);
          }
      }
#endif
  }

  int injection_count = 0;

  std::vector<double> history;
  int iter = 0;
  bool must_stop = false;
  double last_avg_distance = std::numeric_limits<double>::infinity();

  int dims_per_swarm = total_dim / get_num_subswarms();
  int remainder = total_dim % get_num_subswarms();

  // Main optimization loop
  while (!must_stop) {
    iter++;
    stop_manager.increment_iterations();

    // Periodically shuffle the dimensions among the sub-swarms
    if (iter > 1 && iter % get_shuffle_freq() == 0) {
#if CPSO_HAVE_MPI
      // Synchronize all processes before shuffling
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      // Rank 0 creates a permutation of the dimensions
      std::vector<int> permutation(total_dim);
      if (mpi_rank == 0) {
        std::iota(permutation.begin(), permutation.end(), 0);
        std::shuffle(permutation.begin(), permutation.end(), global_gen);
      }
#if CPSO_HAVE_MPI
      // Send the permutation to all the MPI ranks
      MPI_Bcast(permutation.data(), total_dim, MPI_INT, 0, MPI_COMM_WORLD);
#endif

      // Update the active dimensions for each sub-swarm
      int current_dim_start = 0;
      for (int i = 0; i < get_num_subswarms(); ++i) {
        int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);
        std::vector<int> new_active_dims(swarm_dims);

        // Fill the new active dimensions
        for (int d = 0; d < swarm_dims; ++d) {
          new_active_dims[d] = permutation[current_dim_start + d];
        }
        current_dim_start += swarm_dims;

        // Update the active dimensions for the current sub-swarm
        bool is_owned = (i >= local_start_idx && i < local_end_idx);
        swarms[i].update_active_dims(new_active_dims, context, gens[i], is_owned);
      }
    }

    std::vector<std::vector<double>> subswarm_bests(get_num_subswarms());
    std::vector<double> subswarm_best_vals(get_num_subswarms(), std::numeric_limits<double>::infinity());
    double progress_ratio =
        (double)stop_manager.get_current_iters() / stop_manager.get_max_iters();

    if (progress_ratio > 1.0)
      progress_ratio = 1.0;

    // Calculate the current inertia weight
    double current_w = get_w_max() - (get_w_max() - get_w_min()) * progress_ratio;

    // Main optimization loop in micro-batches
    for (int b = 0; b < global_max_local_swarms; ++b) {
        std::vector<double> base_context = context.get_full_vector();
        double base_fitness = context.get_best_fitness();
        std::vector<double> delta_context(total_dim, 0.0);
        
        int i = local_start_idx + b;
        if (i < local_end_idx) {
            swarms[i].recalculate_fitness(context, f);
            swarms[i].update_velocities_and_positions(current_w, get_c1(), get_c2(), gens[i], progress_ratio);
            swarms[i].evaluate_and_update(context, f);

            const auto &active_dims = swarms[i].get_active_dims();
            for (size_t d = 0; d < active_dims.size(); ++d) {
                delta_context[active_dims[d]] = swarms[i].get_gbest_pos()[d] - base_context[active_dims[d]];
            }
        }

#if CPSO_HAVE_MPI
        // Gather all local delta updates from all sub-swarms across MPI ranks
        std::vector<double> global_delta(total_dim, 0.0);
        MPI_Allreduce(delta_context.data(), global_delta.data(), total_dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Apply all deltas together
        std::vector<double> new_full_vector = base_context;
        for (int dim = 0; dim < total_dim; ++dim) {
            new_full_vector[dim] += global_delta[dim];
        }
        
        double new_true_fitness = 0.0;
        if (mpi_rank == 0) {
            new_true_fitness = f.value(new_full_vector);
        }
        MPI_Bcast(&new_true_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if (new_true_fitness < base_fitness) {
            // All merged updates improved the function
            context.set_full_vector(new_full_vector, new_true_fitness);
        } else {
            // The full combination failed
            std::vector<double> salvaged_delta(total_dim, 0.0);
            if (i < local_end_idx && swarms[i].get_gbest_val() < base_fitness) {
                const auto &active_dims = swarms[i].get_active_dims();
                const auto &best_pos = swarms[i].get_gbest_pos();
                for (size_t d = 0; d < active_dims.size(); ++d) {
                    salvaged_delta[active_dims[d]] = best_pos[d] - base_context[active_dims[d]];
                }
            }
          
            std::vector<double> global_salvaged_delta(total_dim, 0.0);
            MPI_Allreduce(salvaged_delta.data(), global_salvaged_delta.data(), total_dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          
            std::vector<double> final_salvaged_vector = base_context;
            for (int dim = 0; dim < total_dim; ++dim) {
                final_salvaged_vector[dim] += global_salvaged_delta[dim];
            }
          
            double final_fitness = 0.0;
            if (mpi_rank == 0) {
                final_fitness = f.value(final_salvaged_vector);
            }
            MPI_Bcast(&final_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
            if (final_fitness < base_fitness) {
                // Merging only the proven "best" positions works
                context.set_full_vector(final_salvaged_vector, final_fitness); 
            } else {
                // Even the proven bests conflict when merged
                std::vector<double> current_greedy = base_context;
                double current_greedy_fitness = base_fitness;

                std::vector<double> local_imps(get_num_subswarms(), 0.0);
                if (i < local_end_idx && swarms[i].get_gbest_val() < base_fitness) {
                    local_imps[i] = base_fitness - swarms[i].get_gbest_val();
                }
                
                std::vector<double> global_imps(get_num_subswarms(), 0.0);
                MPI_Allreduce(local_imps.data(), global_imps.data(), get_num_subswarms(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

                std::vector<int> sorted_swarms(get_num_subswarms());
                std::iota(sorted_swarms.begin(), sorted_swarms.end(), 0);
                std::sort(sorted_swarms.begin(), sorted_swarms.end(), [&](int a, int b) {
                    return global_imps[a] > global_imps[b];
                });

                for (int swarm_idx : sorted_swarms) {
                    if (global_imps[swarm_idx] <= 0.0) continue; 

                    std::vector<double> test_greedy = current_greedy;
                    const auto& active_dims = swarms[swarm_idx].get_active_dims();
                    for (int d : active_dims) {
                       test_greedy[d] += global_delta[d];
                    }

                    double test_fitness = 0.0;
                    if (mpi_rank == 0) {
                        test_fitness = f.value(test_greedy);
                    }
                    MPI_Bcast(&test_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    
                    if (test_fitness < current_greedy_fitness) {
                        current_greedy = test_greedy;
                        current_greedy_fitness = test_fitness;
                    }
                }
                context.set_full_vector(current_greedy, current_greedy_fitness);
            }
        }
#else
        if (i < local_end_idx) {
            std::vector<double> local_vector = base_context;
            for (int dim = 0; dim < total_dim; ++dim) {
                local_vector[dim] += delta_context[dim];
            }
            double local_best = f.value(local_vector);
            
            if (local_best < base_fitness) {
                context.set_full_vector(local_vector, local_best);
            }
        }
#endif
    }

    double current_best_fitness = context.get_best_fitness();
    const std::vector<double> &current_gbest_pos = context.get_full_vector();

    double current_normalized_error = f.error(current_gbest_pos);
    history.push_back(current_normalized_error);

    // Compute average distance between particles and global best position
    SwarmMetrics::compute_avg_distance(iter, swarms, current_gbest_pos, last_avg_distance, local_start_idx, local_end_idx, true);

    bool local_stop = stop_manager.should_stop(current_best_fitness, last_avg_distance);

    // Synchronize the stop condition across all processes
#if CPSO_HAVE_MPI
    int local_stop_int = local_stop ? 1 : 0;
    int global_stop_int = 0;
    MPI_Allreduce(&local_stop_int, &global_stop_int, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    must_stop = (global_stop_int > 0);
#else
    must_stop = local_stop;
#endif

    // Inject random velocities if stagnation patience is exceeded
    int current_stag_iters = stop_manager.get_current_stagnation_iters();
    if (!must_stop && current_stag_iters > 0 && (current_stag_iters % get_stagnation_patience() == 0)) {
      injection_count++;
      bool hard_reset = (injection_count % 3 == 0);
      for (int i = local_start_idx; i < local_end_idx; ++i) {
        swarms[i].inject_velocities(gens[i], hard_reset);
        if (hard_reset) {
            swarms[i].reset_gbest_attraction();
        }
      }
    }
  }

  OutputObject out(f.get_name(), total_dim, get_particles_per_swarm() * get_num_subswarms(),
                   context.get_full_vector(), f.get_true_solution(),
                   context.get_best_fitness(), history, 1, 0.0,
                   iter, stop_manager);

  return out;
}
