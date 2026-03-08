#ifndef METHODS_HPP
#define METHODS_HPP

#include "interfaces.hpp"

/**
 * @brief Distributed DMS-PSO-HS via MPI.
 * 
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles across all ranks.
 * @param max_iter Maximum number of iterations.
 * @return OutputObject containing optimization results.
 */
OutputObject dpso(const TestFunction& f, 
                     unsigned int dim, 
                     unsigned int n_points_total, 
                     int max_iter);

/**
 * @brief Serial DMS-PSO-HS algorithm.
 * 
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles in the swarm.
 * @param max_iter Maximum number of iterations.
 * @return OutputObject containing optimization results.
 */
OutputObject dpso_serial(const TestFunction& f, 
                         unsigned int dim, 
                         unsigned int n_points_total, 
                         int max_iter);

#endif