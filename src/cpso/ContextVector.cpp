/**
 * @file ContextVector.cpp
 * @brief Implementation of the ContextVector class methods
 */

#include "ContextVector.hpp"
#include <limits>
#include <stdexcept>

// Constructor for the ContextVector class
ContextVector::ContextVector(int total_dims)
    : full_vector(total_dims, 0.0),
      best_fitness(std::numeric_limits<double>::infinity()),
      eval_scratch(total_dims, 0.0) {}


double ContextVector::evaluate_particle(const TestFunction &f,
                                 const std::vector<double> &partial_pos,
                                 const std::vector<int> &active_dims) const {
   
  if (partial_pos.size() != active_dims.size()) {
    throw std::invalid_argument(
        "Size mismatch between sizes");
  }

  for (size_t i = 0; i < active_dims.size(); ++i) {
    eval_scratch[active_dims[i]] = partial_pos[i];
  }

  // Evaluate the modified vector
  double fitness = f.value(eval_scratch);

  for (size_t i = 0; i < active_dims.size(); ++i) {
    eval_scratch[active_dims[i]] = full_vector[active_dims[i]];
  }


  return fitness;
}



// Setter method
void ContextVector::set_full_vector(const std::vector<double> &new_vector,
                                    double fitness) {
  full_vector = new_vector;
  eval_scratch = new_vector;
  best_fitness = fitness;
}


// Getters
const std::vector<double> &ContextVector::get_full_vector() const {
  return full_vector;
}
double ContextVector::get_best_fitness() const { return best_fitness; }