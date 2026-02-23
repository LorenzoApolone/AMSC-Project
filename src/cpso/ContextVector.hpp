#pragma once

#include "../interfaces.hpp"
#include <vector>

/**
 * @class ContextVector
 * @brief Manages a global n-dimensional vector for Cooperative PSO variants.
 *
 * The context vector allows sub-swarms (which explore subsets of dimensions)
 * to evaluate their fitness by merging their partial positions with the best
 * known state of the entire multi-dimensional problem.
 */
class ContextVector {
private:
  std::vector<double> full_vector;
  double best_fitness;

public:
  // Constructs a vector of size total_dims initialized to 0
  ContextVector(int total_dims);

  // Getters
  const std::vector<double> &get_full_vector() const;
  double get_best_fitness() const;

  /**
   * @brief Evaluates a partial position by merging it with the Context Vector.
   * @param f The test function to evaluate on.
   * @param partial_pos The position of the particle in K dimensions
   * (Sub-swarm).
   * @param active_dims The global indices of the K dimensions optimized by the
   * sub-swarm.
   * @return The calculated fitness on the entire n-dimensional pseudo-vector.
   */
  double evaluate_particle(const TestFunction &f,
                           const std::vector<double> &partial_pos,
                           const std::vector<int> &active_dims) const;

  /**
   * @brief Updates the Context Vector with a new partial best if it brings
   * improvements.
   */
  void update(const std::vector<double> &partial_best_pos,
              const std::vector<int> &active_dims, double new_fitness);

  /**
   * @brief Completely sets the vector and fitness (mostly used during
   * initialization).
   */
  void set_full_vector(const std::vector<double> &new_vector, double fitness);
};
