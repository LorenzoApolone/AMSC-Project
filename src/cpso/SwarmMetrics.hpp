/**
 * @file SwarmMetrics.hpp
 * @brief Header file for swarm-related metric computations
 */

#pragma once

#include "SubSwarm.hpp"
#include <vector>

/**
 * @namespace SwarmMetrics
 * @brief Namespace containing metric evaluation functions for particle swarms
 */
namespace SwarmMetrics {

  /**
   * @brief Computes the average Euclidean distance between particles and the current global best position
   * 
   * @param swarms The list of sub-swarms
   * @param current_gbest_pos The current global best position to compute distance against
   * @param last_avg_distance Output reference to store the computed average distance
   * @param local_start_idx The starting index of sub-swarms assigned to this execution (for MPI)
   * @param local_end_idx The ending index of sub-swarms assigned to this execution (for MPI)
   * @param use_mpi Boolean flag indicating if MPI reduction should be used to aggregate data globally
   * @param mpi_allreduce_time_s Optional accumulator for MPI_Allreduce timing in the parallel solver
   */
  void compute_avg_distance(const std::vector<SubSwarm>& swarms, 
                            const std::vector<double>& current_gbest_pos, 
                            double& last_avg_distance, 
                            int local_start_idx, int local_end_idx, 
                            bool use_mpi,
                            double *mpi_allreduce_time_s = nullptr);

}
