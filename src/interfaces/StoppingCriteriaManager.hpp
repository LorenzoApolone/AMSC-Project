#pragma once

#include <cmath>
#include <limits>

/**
 * @class StoppingCriteriaManager
 * @brief Manages the termination criteria for the PSO/CPSO algorithms.
 */
class StoppingCriteriaManager {
private:
  int max_iters;
  int max_stagnation_iters;

  double stagnation_tolerance;

  double diversity_tolerance;

  int current_iters;
  int current_stagnation_iters;
  double last_best_fitness;

public:
  StoppingCriteriaManager(int max_iters_limit, int stagnation_iters_limit = 50,
                          double stagnation_tol = 1e-6,
                          double diversity_tol = 1e-4)
      : max_iters(max_iters_limit),
        max_stagnation_iters(stagnation_iters_limit),
        stagnation_tolerance(stagnation_tol),
        diversity_tolerance(diversity_tol), current_iters(0),
        current_stagnation_iters(0),
        last_best_fitness(std::numeric_limits<double>::infinity()) {}

  /**
   * @brief Base evaluations for stopping criteria (fevals and stagnation).
   */
  bool should_stop_base(double current_best_fitness) {
    // Maximum Iterations
    if (current_iters >= max_iters) {
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

    return false;
  }

  /**
   * @brief Increments the iterations counter by 1
   */
  void increment_iterations() { current_iters++; }

  int get_current_iters() const { return current_iters; }

  int get_max_iters() const { return max_iters; }

  /**
   * @brief Evaluates whether the algorithm should stop using an explicitly
   * calculated average distance.
   *
   * @param current_best_fitness The current global fitness.
   * @param explicit_avg_distance The explicitly calculated average distance of
   * the swarm particles.
   * @return true if at least one of the stopping criteria is reached, false
   * otherwise.
   */
  bool should_stop(double current_best_fitness, double explicit_avg_distance) {

    if (should_stop_base(current_best_fitness)) {
      return true;
    }

    if (explicit_avg_distance >= 0.0 &&
        explicit_avg_distance < diversity_tolerance) {
      return true;
    }

    return false;
  }
};
