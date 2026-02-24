#pragma once

#include <cmath>
#include <limits>
#include <vector>

/**
 * @class StoppingCriteriaManager
 * @brief Manages the termination criteria for the PSO/CPSO algorithms.
 */
class StoppingCriteriaManager {
private:
  int max_fevals;
  int max_stagnation_iters;

  double stagnation_tolerance;

  double diversity_tolerance;

  int current_fevals;
  int current_stagnation_iters;
  double last_best_fitness;

public:
  StoppingCriteriaManager(int max_fevals_limit, int stagnation_iters_limit = 50,
                          double stagnation_tol = 1e-6,
                          double diversity_tol = 1e-4)
      : max_fevals(max_fevals_limit),
        max_stagnation_iters(stagnation_iters_limit),
        stagnation_tolerance(stagnation_tol),
        diversity_tolerance(diversity_tol), current_fevals(0),
        current_stagnation_iters(0),
        last_best_fitness(std::numeric_limits<double>::infinity()) {}

  /**
   * @brief Increments the function evaluations (FEvals) counter
   */
  void add_evaluations(int num_evals) { current_fevals += num_evals; }

  int get_current_fevals() const { return current_fevals; }

  int get_max_fevals() const { return max_fevals; }

  /**
   * @brief Evaluates whether the algorithm should stop.
   *
   * @param current_best_fitness The current global fitness (Context Vector
   * or gBest).
   * @param swarm_positions The current positions of the particles in the swarm
   * (to calculate the diversity). In Parallel/Serial CPSO, this can be
   * the active portion of the dimensions.
   * @param target_pos The position of the gBest or Context Vector associated
   * with the fitness.
   * @return true if at least one of the stopping criteria is reached, false
   * otherwise.
   */
  bool should_stop(double current_best_fitness,
                   const std::vector<std::vector<double>> &swarm_positions,
                   const std::vector<double> &target_pos) {

    // Maximum Function Evaluations
    if (current_fevals >= max_fevals) {
      return true;
    }

    // Stagnation Control
    if (std::abs(last_best_fitness - current_best_fitness) <
        stagnation_tolerance) {
      current_stagnation_iters++;
    } else {
      // Significant improvement, reset the counter if there is a significant
      // improvement
      current_stagnation_iters = 0;
    }

    last_best_fitness = current_best_fitness;

    if (current_stagnation_iters >= max_stagnation_iters) {
      return true;
    }

    // Swarm Diversity / Convergence
    // Calculation of the average distance of the particles from the target_pos
    if (!swarm_positions.empty() && !target_pos.empty()) {
      double total_distance = 0.0;
      int num_particles = swarm_positions.size();
      int dims = std::min(swarm_positions[0].size(),
                          target_pos.size()); // Prevent out-of-bounds

      for (const auto &pos : swarm_positions) {
        double dist_sq = 0.0;
        for (int j = 0; j < dims; ++j) {
          double diff = pos[j] - target_pos[j];
          dist_sq += diff * diff;
        }
        total_distance += std::sqrt(dist_sq);
      }

      double avg_distance = total_distance / num_particles;
      if (avg_distance < diversity_tolerance) {
        return true;
      }
    }

    return false;
  }
};
