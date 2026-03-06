#include "CPSOBase.hpp"
#include <cmath>
#include <chrono>

#if __has_include(<mpi.h>)
#include <mpi.h>
#endif

// Constructor for uniform topology
CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   NetworkType topology, int shuffle_freq,
                   int stagnation_patience, double w_start,
                   double w_end, double coeff1, double coeff2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      shuffle_freq(shuffle_freq), stagnation_patience(stagnation_patience),
      w_max(w_start), w_min(w_end), c1(coeff1), c2(coeff2) {
  subswarm_topologies.assign(k_subswarms, topology);
}

// Constructor for different types of topologies
CPSOBase::CPSOBase(int k_subswarms, int num_particles_per_swarm,
                   const std::vector<NetworkType> &topologies,
                   int shuffle_freq, int stagnation_patience,
                   double w_start, double w_end, double coeff1, double coeff2)
    : num_subswarms(k_subswarms), particles_per_swarm(num_particles_per_swarm),
      shuffle_freq(shuffle_freq), stagnation_patience(stagnation_patience),
      w_max(w_start), w_min(w_end), c1(coeff1), c2(coeff2) {
  subswarm_topologies = topologies;
  // Padding for the remaining topologies with SCALE_FREE
  while (subswarm_topologies.size() < static_cast<size_t>(k_subswarms)) {
    subswarm_topologies.push_back(topologies.empty() ? NetworkType::SCALE_FREE
                                                     : topologies.front());
  }
}

// Method that computes the average distance between particles and the global best position
void CPSOBase::compute_avg_distance(int iter, const std::vector<SubSwarm>& swarms, 
                                    const std::vector<double>& current_gbest_pos, 
                                    double& last_avg_distance, 
                                    int local_start_idx, int local_end_idx, bool use_mpi) {

  // Compute average distance on the first iter, and then every 10 iterations
  if (iter == 1 || iter % 10 == 0) {
    std::vector<double> local_dist_sq(particles_per_swarm, 0.0);
    for (int p = 0; p < particles_per_swarm; ++p) {
      for (int i = local_start_idx; i < local_end_idx; ++i) {
        const auto &particle = swarms[i].get_particles()[p];
        const auto &active_dims = swarms[i].get_active_dims();
        for (size_t d = 0; d < active_dims.size(); ++d) {
          double diff = particle.position[d] - current_gbest_pos[active_dims[d]];
          local_dist_sq[p] += diff * diff;
        }
      }
    }

    std::vector<double> global_dist_sq(particles_per_swarm, 0.0);
    
    // Use MPI in the parallel case, otherwise it uses the local distances
    if (use_mpi) {
#ifdef MPI_VERSION
      MPI_Allreduce(local_dist_sq.data(), global_dist_sq.data(),
                    particles_per_swarm, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    } else {
      global_dist_sq = local_dist_sq;
    }

    double total_distance = 0.0;
    for (int p = 0; p < particles_per_swarm; ++p) {
      total_distance += std::sqrt(global_dist_sq[p]);
    }
    last_avg_distance = total_distance / particles_per_swarm;
  }
}

// Optimization loop
OutputObject CPSOBase::optimize(const TestFunction &f, StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();

  // Generating random number generators for each subswarm
  std::vector<std::mt19937> gens(num_subswarms);
  for (int i = 0; i < num_subswarms; ++i) {
    gens[i] = std::mt19937(1337 + i);
  }
  std::mt19937 global_gen(1337);

  // Domain Decomposition depending on how many subswarm we have
  int total_dim = f.dim;
  auto bounds = f.get_domain();
  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  std::vector<SubSwarm> swarms;
  swarms.reserve(num_subswarms);

  // Creating the subswarms with the active dimensions and the topology
  int current_dim_start = 0;
  for (int i = 0; i < num_subswarms; ++i) {
    int swarm_dims = dims_per_swarm + (i < remainder ? 1 : 0);
    std::vector<int> active_dims;
    for (int d = 0; d < swarm_dims; ++d) {
      active_dims.push_back(current_dim_start + d);
    }
    current_dim_start += swarm_dims;

    // Creating the adjacency list for the current subswarm topology
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

    swarms.emplace_back(particles_per_swarm, active_dims, bounds.first, bounds.second, sub_adj_list);
  }

  // Creating the Global Context Vector
  ContextVector context(total_dim);
  std::vector<double> init_vec(total_dim);
  
  // Set the initial global best position to the mid-point of the domain
  context.set_full_vector(std::vector<double>(total_dim, bounds.first + (bounds.second - bounds.first) / 2.0),
                          std::numeric_limits<double>::infinity());

  OutputObject out = run_optimization_loop(f, stop_manager, swarms, gens, context, global_gen);

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  out.execution_time = elapsed.count();
  
  return out;
}
