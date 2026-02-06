/**
 * @file pso_mpi.cpp
 * @brief Implementation of the Parallel Particle Swarm Optimization (PSO) using MPI.
 */

#include "interfaces.hpp"
#include <algorithm>
#include <limits>
#include <mpi.h>
#include <random>
#include <vector>

/**
 * @struct PSOHyperparameters
 * @brief Constants and configuration parameters for the PSO algorithm
 */
struct PSOHyperparameters {
  /// @brief Cognitive coefficient (local best influence)
  static constexpr double C1 = 1.49618;
  
  /// @brief Social coefficient (global best influence)
  static constexpr double C2 = 1.49618;

  // Dynamic inertia configuration
  /// @brief Maximum inertia weight (start of exploration)
  static constexpr double W_MAX = 0.9;
  /// @brief Minimum inertia weight (end of exploitation)
  static constexpr double W_MIN = 0.4;

  /// @brief Factor used to scale initial particle velocities
  static constexpr double V_INIT_FACTOR = 0.1;
};

/**
 * @brief Solves the optimization problem using parallel PSO
 * * @param f The test function (so the problem) to solve
 * @param d The dimensionality of the problem
 * @param stop The stopping criterion configuration
 * @param n_points Total number of particles in the swarm
 * @return OutputObject containing results and metrics
 */
OutputObject pso_mpi(const TestFunction &f, int d, const StopCriterion &stop,
                     int n_points, bool &converged) {

  // --- MPI Setup ---
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double start_time = MPI_Wtime();

  // Particle Distribution
  /// @brief Calculate how many particles this specific core will manage
  int local_n = n_points / size;
  int remainder = n_points % size;
  if (rank < remainder)
    local_n++;

  // Local Data Structures
  std::vector<std::vector<double>> pos(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> vel(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> pbest_pos = pos;
  std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

  // Global Best (Shared Knowledge)
  std::vector<double> gbest_pos(d);
  double gbest_val = std::numeric_limits<double>::max();

  // History (Rank 0 only)
  std::vector<double> history;

  // Initialization
  std::pair<double, double> bounds = f.get_domain();
  double LB = bounds.first;
  double UB = bounds.second;

  /// @brief Seed for random number generator (using different seeds per rank is crucial)
  std::mt19937 gen(rank + 42);
  std::uniform_real_distribution<> dis(LB, UB);
  std::uniform_real_distribution<> dis_01(0.0, 1.0);

  for (int i = 0; i < local_n; ++i) {
    for (int j = 0; j < d; ++j) {
      pos[i][j] = dis(gen);
      /// @brief Initializes velocity as a fraction of search space to prevent explosions
      vel[i][j] = (dis(gen) - dis(gen)) * PSOHyperparameters::V_INIT_FACTOR;
      pbest_pos[i][j] = pos[i][j];
    }

    // Initial Evaluation
    double fitness = f.value(pos[i]);
    pbest_val[i] = fitness;

    if (fitness < gbest_val) {
      gbest_val = fitness;
      gbest_pos = pos[i];
    }
  }

  // Initial Synchronization
  /// @brief Ensures everyone starts with the same global best
  struct {
    double val;
    int rank;
  } loc_data, glob_data;
  loc_data.val = gbest_val;
  loc_data.rank = rank;

  MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                MPI_COMM_WORLD);
  gbest_val = glob_data.val;
  MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

  // Main Optimization Loop
  int iter = 0;
  bool must_stop = false;

  /// @brief We use max_iter to calculate the inertia fraction
  int max_iter_limit = stop.get_max_iter();

  while (!must_stop) {

    /// @brief We calculatee dynamic inertia weight (which is linearlly decreasing)
    double current_w =
        PSOHyperparameters::W_MAX -
        ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) * iter /
         max_iter_limit);

    // Update Particles (this is the local Work)
    for (int i = 0; i < local_n; ++i) {
      for (int j = 0; j < d; ++j) {
        double r1 = dis_01(gen);
        double r2 = dis_01(gen);

        /// @brief Update velocity using current_w, C1, and C2
        vel[i][j] =
            current_w * vel[i][j] +
            PSOHyperparameters::C1 * r1 * (pbest_pos[i][j] - pos[i][j]) +
            PSOHyperparameters::C2 * r2 * (gbest_pos[j] - pos[i][j]);

        // Update position
        pos[i][j] += vel[i][j];

        /// @brief Boundary handling (it clamps to the edges)
        if (pos[i][j] < LB)
          pos[i][j] = LB;
        if (pos[i][j] > UB)
          pos[i][j] = UB;
      }

      // Evaluate fitness
      double current_fit = f.value(pos[i]);

      // Update personal best
      if (current_fit < pbest_val[i]) {
        pbest_val[i] = current_fit;
        pbest_pos[i] = pos[i];
      }

      // Update local view of global best
      if (pbest_val[i] < gbest_val) {
        gbest_val = pbest_val[i];
        gbest_pos = pbest_pos[i];
      }
    }

    // Synchronization (Communication step)

    /// @brief Find the absolute best value among all ranks
    loc_data.val = gbest_val;
    loc_data.rank = rank;
    MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  MPI_COMM_WORLD);

    /// @brief Share the coordinates of that best value
    if (glob_data.val < gbest_val) {
      gbest_val = glob_data.val;
      MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank,
                MPI_COMM_WORLD);
    } else {
      MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank,
                MPI_COMM_WORLD);
    }

    // Check Stopping Criterion (Rank 0 decides)
    int stop_signal = 0;
    if (rank == 0) {
      double current_normalized_error = f.error(gbest_pos);
      history.push_back(current_normalized_error);
      
      /// @brief Check tolerance using patience
      if (stop.should_stop(iter, current_normalized_error)) {
        stop_signal = 1;
        if (iter < stop.get_max_iter())
          converged = true;
      }
    }
    MPI_Bcast(&stop_signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (stop_signal)
      must_stop = true;

    iter++;
  }

  double end_time = MPI_Wtime();

  // Return Result
  OutputObject out(f.get_name(),          // function_name
                   d,                     // dimension
                   n_points,              // number of particles
                   gbest_pos,             // best position found
                   f.get_true_solution(), // true solution
                   gbest_val,             // best fitness value
                   history,               // convergence history
                   size,                  // number of cores (serial)
                   end_time - start_time, // execution time
                   iter,                  // number of iterations
                   stop                   // stop criterion
  );

  return out;
}