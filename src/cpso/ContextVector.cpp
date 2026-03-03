/**
 * @file ContextVector.cpp
 * @brief Implementation of the ContextVector class methods.
 */

#include "ContextVector.hpp"
#include <limits>
#include <stdexcept>

// Constructor for the ContextVector class.
ContextVector::ContextVector(int total_dims)
    : full_vector(total_dims, 0.0),
      best_fitness(std::numeric_limits<double>::infinity()) {}


double ContextVector::evaluate_particle(const TestFunction &f,
                                 const std::vector<double> &partial_pos,
                                 const std::vector<int> &active_dims) const {
  
  // Safety Check 
  if (partial_pos.size() != active_dims.size()) {
    throw std::invalid_argument(
        "Size mismatch between sizes");
  }

  // Create a temporary copy of the full vector
  std::vector<double> temp_vector = full_vector;

  // Overwrite the vector with target dimensions
  for (size_t i = 0; i < active_dims.size(); ++i) {
    temp_vector[active_dims[i]] = partial_pos[i];
  }

  // Find the fitness the temporary vector
  return f.value(temp_vector);
}

void ContextVector::update(const std::vector<double> &partial_best_pos,
                           const std::vector<int> &active_dims,
                           double new_fitness) {
  
  // We only replace if the new fitness is better than the best one
  if (new_fitness <= best_fitness) {

    for (size_t i = 0; i < active_dims.size(); ++i) {
      full_vector[active_dims[i]] = partial_best_pos[i];
    }

    best_fitness = new_fitness;

  }
}

// Setter method
void ContextVector::set_full_vector(const std::vector<double> &new_vector,
                                    double fitness) {
  full_vector = new_vector;
  best_fitness = fitness;
}


// Getters
const std::vector<double> &ContextVector::get_full_vector() const {
  return full_vector;
}
double ContextVector::get_best_fitness() const { return best_fitness; }