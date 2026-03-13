#include "SwarmMetrics.hpp"
#include <cmath>

#if __has_include(<mpi.h>)
#include <mpi.h>
#define CPSO_HAVE_MPI 1
#else
#define CPSO_HAVE_MPI 0
#endif

namespace SwarmMetrics {

void compute_avg_distance(int iter, const std::vector<SubSwarm>& swarms, 
                          const std::vector<double>& current_gbest_pos, 
                          double& last_avg_distance, 
                          int local_start_idx, int local_end_idx, bool use_mpi) {

  // Compute average distance on the first iter, and then every 10 iterations
  if (iter == 1 || iter % 10 == 0) {
    double local_sum_dist = 0.0;
    int local_particle_count = 0;

    for (int i = local_start_idx; i < local_end_idx; ++i) {
      // Extract the active dimensions and positions of the sub-swarm
      const auto &active_dims = swarms[i].get_active_dims();
      const auto &positions = swarms[i].get_positions();
      int current_particles_count = swarms[i].get_num_particles();
      int dim = active_dims.size();
      for (int p = 0; p < current_particles_count; ++p) {
        double dist_sq = 0.0;
        for (int d = 0; d < dim; ++d) {
          // Calculate the distance in the active dimensions
          double diff = positions[p * dim + d] - current_gbest_pos[active_dims[d]];
          dist_sq += diff * diff;
        }
        // Add the distance to the local sum
        local_sum_dist += std::sqrt(dist_sq);
        // Add the number of particles to the local count
        local_particle_count++;
      }
    }

    double global_sum_dist = 0.0;
    int global_particle_count = 0;
    

    if (use_mpi && !CPSO_HAVE_MPI) {
      use_mpi = false;
    }

    // use MPI in the parallel case, otherwise it uses the local distances
    if (use_mpi) {
#if CPSO_HAVE_MPI
      MPI_Allreduce(&local_sum_dist, &global_sum_dist, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_particle_count, &global_particle_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    } else {
      global_sum_dist = local_sum_dist;
      global_particle_count = local_particle_count;
    }

    last_avg_distance = (global_particle_count > 0) ? (global_sum_dist / global_particle_count) : 0.0;
  }
}

} // namespace SwarmMetrics
