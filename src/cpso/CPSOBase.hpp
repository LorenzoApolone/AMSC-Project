#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "../topology/create_network.hpp"
#include "ContextVector.hpp"
#include "SubSwarm.hpp"
#include <random>
#include <vector>

/**
 * @class CPSOBase
 * @brief Abstract base class for the CPSO algorithm
 */
class CPSOBase {
protected:

  // Number of sub-swarms
  int num_subswarms;

  // Number of particles per sub-swarm
  int particles_per_swarm;

  // Frequency at which dimensions are shuffled among sub-swarms
  int shuffle_freq;

  // Number of iterations without improvement before injecting random velocities
  int stagnation_patience;

  // Maximum and minimum inertia weight
  double w_max, w_min;

  // Cognitive and social learning factors
  double c1, c2;

  // Network topologies
  std::vector<NetworkType> subswarm_topologies;

  // Method to compute the average distance between particles in the swarm
  void compute_avg_distance(int iter, const std::vector<SubSwarm>& swarms, 
                            const std::vector<double>& current_gbest_pos, 
                            double& last_avg_distance, 
                            int local_start_idx, int local_end_idx, 
                            bool use_mpi);

public:

  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           NetworkType topology, int shuffle_freq,
           int stagnation_patience, double w_start,
           double w_end, double coeff1, double coeff2);

  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           const std::vector<NetworkType> &topologies,
           int shuffle_freq, int stagnation_patience,
           double w_start, double w_end, double coeff1, double coeff2);

  virtual ~CPSOBase() = default;

  /**
   * @brief Executes the optimization process
   */
  OutputObject optimize(const TestFunction &f, StoppingCriteriaManager &stop_manager);

protected:
  /**
   * @brief Run method overridable by Parallel or Serial implementations
   */
  virtual OutputObject run_optimization_loop(const TestFunction &f, StoppingCriteriaManager &stop_manager,
                                             std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                             ContextVector& context, std::mt19937& global_gen) = 0;
};
