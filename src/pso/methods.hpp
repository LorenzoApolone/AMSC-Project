#ifndef METHODS_HPP
#define METHODS_HPP
/**
 * @file methods.hpp
 * @brief Public entry points for the standard serial and MPI PSO solvers.
 */

#include "../interfaces/interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include <vector>

/**
 * @brief Run the standard serial Particle Swarm Optimization solver.
 *
 * @param f Objective function to optimize.
 * @param d Dimension of the search space.
 * @param stop Stopping criteria manager shared with the benchmark driver.
 * @param n_points Number of particles in the swarm.
 * @param seed Random seed used to initialize particles and velocities.
 * @return OutputObject containing best solution, history, timing and metadata.
 */
OutputObject pso_serial(const TestFunction &f, int d,
                        StoppingCriteriaManager &stop, int n_points,
                        unsigned int seed = 12345);

/**
 * @brief Run the standard MPI Particle Swarm Optimization solver.
 *
 * @param f Objective function to optimize.
 * @param d Dimension of the search space.
 * @param stop Stopping criteria manager shared with the benchmark driver.
 * @param n_points Total number of particles distributed across MPI ranks.
 * @param seed Random seed used to initialize particles and velocities.
 * @return OutputObject containing best solution, history, timing and metadata.
 */
OutputObject pso_mpi(const TestFunction &f, int d,
                     StoppingCriteriaManager &stop, int n_points,
                     unsigned int seed = 12345);

#endif
