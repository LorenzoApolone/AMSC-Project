#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include <vector>

/**
 * @file pso_topology.hpp
 * @brief Parallel Particle Swarm Optimization with network topology.
 *
 * This module execute the PSO algorithm with a network topology where
 * particles communicate according to a predefined network topology 
 * passed by reference. The algorithm run with MPI and distributes particles
 * across processes.
 *
 * Each particle updates its velocity and position based on:
 * - its personal best (pbest)
 * - the best particle in its neighborhood (lbest)
 *
 * The neighborhood structure is defined by an adjacency list representing
 * the communication topology.
 */


/**
 * @brief Execute the topology-based parallel PSO algorithm.
 *
 * This function runs a Particle Swarm Optimization (PSO) algorithm in parallel
 * using MPI. Particles are distributed across processes and exchange their
 * personal-best information in order to compute neighborhood-best solutions
 * according to a given communication topology.
 *
 * The algorithm iteratively updates particle velocities and positions until
 * a stopping criterion is satisfied.
 *
 * @param f Test function.
 * @param d Dimensionality of the search space.
 * @param stop Manager handling the stopping criteria of the algorithm.
 * @param n_points Total number of particles in the swarm.
 * @param adjacency_list Adjacency list defining the communication topology
 *        between particles. Each index corresponds to a particle and stores
 *        the list of its neighbors.
 * @param t_allgatherv_tot Accumulator for the total communication time spent
 *        in MPI_Allgatherv operations.
 *
 * @return OutputObject containing the results of the optimization,
 * including the best solution found, its fitness value, iteration history,
 * and execution statistics.
 */

OutputObject pso_topology(const TestFunction &f,
                       int d,
                        StoppingCriteriaManager &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list,  double &t_allgatherv_tot);
