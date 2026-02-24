#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "ContextVector.hpp"
#include "SubSwarm.hpp"
#include <random>
#include <vector>

/**
 * @class CPSOSerial
 * @brief Serial Implementation of Cooperative Particle Swarm Optimization
 * (CPSO-S).
 *
 * CPSO-S divides the D-dimensional problem into K
 * one-dimensional (or reduced-dimensionality) sub-swarms. Each sub-swarm
 * is updated and evaluated sequentially. After each sub-swarm's turn,
 * the partial best found immediately contributes to
 * updating the Context Vector.
 */
class CPSOSerial {
private:
  int num_subswarms;
  int particles_per_swarm;
  double w_max, w_min, c1, c2;

public:
  CPSOSerial(int k_subswarms, int num_particles_per_swarm, double w_start = 0.9,
             double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

  /**
   * @brief Executes CPSO-S optimization
   *
   * @param f The objective function to minimize
   * @param stop_manager Stopping criteria manager (FEvals, Stagnation,
   * etc.)
   * @return OutputObject containing the results, global gBest, and history
   */
  OutputObject optimize(const TestFunction &f,
                        StoppingCriteriaManager &stop_manager);
};
