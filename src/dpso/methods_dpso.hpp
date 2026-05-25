#ifndef METHODS_HPP
#define METHODS_HPP

#include "../interfaces/interfaces.hpp"

/**
 * @brief Parameters for DMS-PSO-HS algorithms
 */
struct DPSOParameters {
  double w = 0.729;
  double c1 = 1.49445;
  double c2 = 1.49445;
  int regrouping_period = 5;
  int sub_swarm_size = 5;
  double hmcr = 0.98;
  double par_min = 0.01;
  double par_max = 0.99;
};

/**
 * @brief Distributed DMS-PSO-HS via MPI.
 *
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles across all ranks.
 * @param max_iter Maximum number of iterations.
 * @param params Hyperparameters for the algorithm.
 * @return OutputObject containing optimization results.
 */
OutputObject dpso(const TestFunction &f, unsigned int dim,
                  unsigned int n_points_total, int max_iter,
                  const DPSOParameters &params = DPSOParameters(),
                  double convergence_tol = 1e-4, unsigned int seed = 0);

/**
 * @brief Serial DMS-PSO-HS algorithm.
 *
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles in the swarm.
 * @param max_iter Maximum number of iterations.
 * @param params Hyperparameters for the algorithm.
 * @return OutputObject containing optimization results.
 */
OutputObject dpso_serial(const TestFunction &f, unsigned int dim,
                         unsigned int n_points_total, int max_iter,
                         const DPSOParameters &params = DPSOParameters(),
                         double convergence_tol = 1e-4, unsigned int seed = 0);

#endif
