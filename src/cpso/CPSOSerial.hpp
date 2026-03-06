/**
 * @file CPSOSerial.hpp
 * @brief Manages the execution of the serial CPSO algorithm
 */
#pragma once

#include "CPSOBase.hpp"

/**
 * @class CPSOSerial
 * @brief Manages the execution of the serial CPSO algorithm
 */
class CPSOSerial : public CPSOBase {
public:
  /**
   * @brief Constructs a CPSOSerial algorithm with a uniform network topology for all sub-swarms
   * 
   * @param k_subswarms Number of sub-swarms the problem is divided into
   * @param num_particles_per_swarm Number of particles residing in each sub-swarm
   * @param topology The network topology uniformly applied to all sub-swarms
   * @param shuffle_freq The iteration frequency at which dimensions are randomly re-assigned across sub-swarms
   * @param stagnation_patience Number of iterations without global improvement before velocities are injected
   * @param w_start Initial inertia weight maximum value for the particles
   * @param w_end Final inertia weight minimum value for the particles
   * @param coeff1 Cognitive learning factor
   * @param coeff2 Social learning factor
   */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm, NetworkType topology,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

  /**
   * @brief Constructs a CPSOSerial algorithm with specific network topologies for each sub-swarm
   * 
   * @param k_subswarms Number of sub-swarms the problem is divided into
   * @param num_particles_per_swarm Number of particles residing in each sub-swarm
   * @param topologies Vector containing the network topologies for each sub-swarm
   * @param shuffle_freq The iteration frequency at which dimensions are randomly re-assigned across sub-swarms
   * @param stagnation_patience Number of iterations without global improvement before velocities are injected
   * @param w_start Initial inertia weight maximum value for the particles
   * @param w_end Final inertia weight minimum value for the particles
   * @param coeff1 Cognitive learning factor
   * @param coeff2 Social learning factor
   */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const std::vector<NetworkType> &topologies,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

protected:
  OutputObject run_optimization_loop(const TestFunction &f,
                                     StoppingCriteriaManager &stop_manager,
                                     std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                     ContextVector& context, std::mt19937& global_gen) override;
};
