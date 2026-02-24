#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "ContextVector.hpp"
#include "SubSwarm.hpp"
#include <random>
#include <vector>

/**
 * @class CPSOParallel
 * @brief Logical Parallel Implementation of Cooperative Particle Swarm
 * Optimization (CPSO-P).
 *
 * In CPSO-P, the K sub-swarms are all evaluated against the SAME
 * Context Vector (snapshotted at the beginning of the iteration). The update
 * of the Context Vector happens in bulk only at the END of the loop over all
 * sub-swarms (logical parallelism).
 */
class CPSOParallel {
private:
  int num_subswarms;
  int particles_per_swarm;
  double w_max, w_min, c1, c2;

public:
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               double w_start = 0.9, double w_end = 0.4,
               double coeff1 = 1.49618, double coeff2 = 1.49618);

  /**
   * @brief Executes the CPSO-P optimization (Logical Parallelism)
   */
  OutputObject optimize(const TestFunction &f,
                        StoppingCriteriaManager &stop_manager);
};
