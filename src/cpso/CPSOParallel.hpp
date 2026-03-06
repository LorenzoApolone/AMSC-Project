#pragma once

#include "CPSOBase.hpp"

/**
 * @class CPSOParallel
 * @brief Manages the execution of the parallel CPSO algorithm
 */
class CPSOParallel : public CPSOBase {

public:

  /**
   * @brief Constructs a CPSOParallel with an uniform topology
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               NetworkType topology, int shuffle_freq = 50,
               int stagnation_patience = 50, double w_start = 0.9,
               double w_end = 0.4, double coeff1 = 1.49618,
               double coeff2 = 1.49618);

  /**
   * @brief Constructs a CPSOParallel optimizer with heterogeneous topologies
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               const std::vector<NetworkType> &topologies,
               int shuffle_freq = 50, int stagnation_patience = 50,
               double w_start = 0.9, double w_end = 0.4,
               double coeff1 = 1.49618, double coeff2 = 1.49618);

protected:
  OutputObject run_optimization_loop(const TestFunction &f, StoppingCriteriaManager &stop_manager,
                                     std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                     ContextVector& context, std::mt19937& global_gen) override;
};
