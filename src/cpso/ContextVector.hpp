/**
 * @file ContextVector.hpp
 * @brief Defines the ContextVector class and methods
 */

#pragma once

#include "../interfaces/interfaces.hpp"
#include <vector>

/**
 * @class ContextVector
 * @brief Manages the full context vector that contains combinations of partial solutions from sub-swarms
 */
class ContextVector {
private:

  // The complete vector representing the global position
  std::vector<double> full_vector;

  // The fitness value associated with the full_vector
  double best_fitness;

  // Scratch buffer for evaluating particles without loop allocations
  mutable std::vector<double> eval_scratch;

public:

  /**
   * @brief Constructs a new ContextVector
   * @param total_dims The total dimensionality of the optimization problem
   */
  ContextVector(int total_dims);

  /**
   * @brief Gets the full context vector
   * @return A constant reference to the entire vector
   */
  const std::vector<double> &get_full_vector() const;

  /**
   * @brief Gets the best fitness value known by the context vector
   * @return The best fitness value
   */
  double get_best_fitness() const;

  /**
   * @brief Evaluates a partial position by replacing the active dimensions the full vector
   * @param f The test function to be evaluated
   * @param partial_pos The partial position to evaluate
   * @param active_dims The indices of the dimensions corresponding to the partial position
   * @return The fitness value of the combined vector
   */
  double evaluate_particle(const TestFunction &f,
                           const std::vector<double> &partial_pos,
                           const std::vector<int> &active_dims) const;



  /**
   * @brief Sets the full context vector and its corresponding fitness limit
   * @param new_vector The new full position vector
   * @param fitness The fitness value of the new position vector
   */
  void set_full_vector(const std::vector<double> &new_vector, double fitness);
};
