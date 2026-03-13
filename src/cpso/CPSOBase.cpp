#include "CPSOBase.hpp"
#include <chrono>

#if __has_include(<mpi.h>)
#include <mpi.h>
#define CPSO_HAVE_MPI 1
#else
#define CPSO_HAVE_MPI 0
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


// Optimization loop
OutputObject CPSOBase::optimize(const TestFunction &f, StoppingCriteriaManager &stop_manager) {
  auto start_time = std::chrono::high_resolution_clock::now();

  // Validation checks
  if (num_subswarms > static_cast<int>(f.dim)) {
      throw std::invalid_argument("number of subswarms must be <= function dimension");
  }

  // Generating random number generators for each subswarm
  std::random_device rd;
  unsigned int master_seed = rd();

#if CPSO_HAVE_MPI
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    MPI_Bcast(&master_seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }
#endif

  std::vector<std::mt19937> gens(num_subswarms);
  for (int i = 0; i < num_subswarms; ++i) {
    gens[i] = std::mt19937(master_seed + 1 + i);
  }
  std::mt19937 global_gen(master_seed);

  // Domain Decomposition depending on how many subswarm we have
  int total_dim = f.dim;
  auto bounds = f.get_domain();
  int dims_per_swarm = total_dim / num_subswarms;
  int remainder = total_dim % num_subswarms;

  // Creating the subswarms array
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
    }

    swarms.emplace_back(particles_per_swarm, active_dims, bounds.first, bounds.second, sub_adj_list);
  }

  // Creating the Global Context Vector
  ContextVector context(total_dim);
  
  // Random Initialization instead of Mid-Point
  std::uniform_real_distribution<double> dist_init(bounds.first, bounds.second);
  std::vector<double> init_vec(total_dim);
  for (int i = 0; i < total_dim; ++i) {
      init_vec[i] = dist_init(global_gen);
  }
  
  context.set_full_vector(init_vec, f.value(init_vec));

  OutputObject out = run_optimization_loop(f, stop_manager, swarms, gens, context, global_gen);

  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  out.execution_time = elapsed.count();
  
  return out;
}
