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
private:

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

protected:

  /** @brief Gets the number of sub-swarms */
  int get_num_subswarms() const { return num_subswarms; }
  
  /** @brief Gets the number of particles per sub-swarm */
  int get_particles_per_swarm() const { return particles_per_swarm; }
  
  /** @brief Gets the frequency at which dimensions are shuffled among sub-swarms */
  int get_shuffle_freq() const { return shuffle_freq; }
  
  /** @brief Gets the number of stagnation iterations before random injection */
  int get_stagnation_patience() const { return stagnation_patience; }
  
  /** @brief Gets the maximum inertia weight */
  double get_w_max() const { return w_max; }
  
  /** @brief Gets the minimum inertia weight */
  double get_w_min() const { return w_min; }
  
  /** @brief Gets the cognitive learning factor */
  double get_c1() const { return c1; }
  
  /** @brief Gets the social learning factor */
  double get_c2() const { return c2; }
  
  /** @brief Gets the network topologies for all the sub-swarms */
  const std::vector<NetworkType>& get_subswarm_topologies() const { return subswarm_topologies; }

public:

  /**
   * @brief Constructs the CPSOBase with a uniform network topology for all sub-swarms
   * 
   * @param k_subswarms The number of sub-swarms
   * @param num_particles_per_swarm The number of particles in each sub-swarm
   * @param topology The network topology to apply to all sub-swarms
   * @param shuffle_freq Iteration frequency for re-shuffling active dimensions
   * @param stagnation_patience Iterations without improvement before injecting random velocities
   * @param w_start The starting (maximum) inertia weight
   * @param w_end The final (minimum) inertia weight
   * @param coeff1 The cognitive learning factor (C1)
   * @param coeff2 The social learning factor (C2)
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           NetworkType topology, int shuffle_freq,
           int stagnation_patience, double w_start,
           double w_end, double coeff1, double coeff2);

  /**
   * @brief Constructs the CPSOBase with specific network topologies for each sub-swarm
   * 
   * @param k_subswarms The number of sub-swarms
   * @param num_particles_per_swarm The number of particles in each sub-swarm
   * @param topologies A vector defining the specific network topology for each sub-swarm
   * @param shuffle_freq Iteration frequency for re-shuffling active dimensions
   * @param stagnation_patience Iterations without improvement before injecting random velocities
   * @param w_start The starting (maximum) inertia weight
   * @param w_end The final (minimum) inertia weight
   * @param coeff1 The cognitive learning factor (C1)
   * @param coeff2 The social learning factor (C2)
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           const std::vector<NetworkType> &topologies,
           int shuffle_freq, int stagnation_patience,
           double w_start, double w_end, double coeff1, double coeff2);

  /**
   * @brief Default virtual destructor
   */
  virtual ~CPSOBase() = default;

  /**
   * @brief Executes the optimization process
   * 
   * @param f The test function to minimize
   * @param stop_manager The criteria manager that signals when to stop the algorithm
   * @return An OutputObject mapping the full execution metrics and results
   */
  OutputObject optimize(const TestFunction &f, StoppingCriteriaManager &stop_manager);

protected:
  /**
   * @brief Protected abstract method defining the specific optimization loop logic
   * 
   * @param f The test function to minimize
   * @param stop_manager Manager to track stopping limits
   * @param swarms The populated array of configured sub-swarms to be processed
   * @param gens Random Number Generators associated independently with each sub-swarm
   * @param context The shared ContextVector caching the aggregated global solution bests
   * @param global_gen Master RNG used for generic global operations (like dimensional shuffling)
   * @return The OutputObject holding execution outcomes
   */
  virtual OutputObject run_optimization_loop(const TestFunction &f, StoppingCriteriaManager &stop_manager,
                                             std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                             ContextVector& context, std::mt19937& global_gen) = 0;
};
