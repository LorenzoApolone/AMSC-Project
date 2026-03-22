/**
 * @file ContextVector.cpp
 * @brief Implementation of the ContextVector class methods
 */

#include "ContextVector.hpp"
#include "NumericValidation.hpp"
#include <limits>
#include <stdexcept>

namespace {

void validate_active_dims(const std::vector<int> &active_dims,
                          size_t total_dims) {
  std::vector<bool> seen(total_dims, false);
  for (int dim : active_dims) {
    if (dim < 0 || static_cast<size_t>(dim) >= total_dims) {
      throw std::out_of_range("active dimension index out of bounds");
    }
    if (seen[dim]) {
      throw std::invalid_argument("active dimensions must be unique");
    }
    seen[dim] = true;
  }
}

} // namespace

ContextVector::ContextVector(int total_dims)
    : full_vector(total_dims > 0 ? static_cast<size_t>(total_dims) : 0u, 0.0),
      best_fitness(std::numeric_limits<double>::infinity()),
      eval_scratch(total_dims > 0 ? static_cast<size_t>(total_dims) : 0u, 0.0) {
  if (total_dims <= 0) {
    throw std::invalid_argument("context vector dimension must be > 0");
  }
}

double ContextVector::evaluate_particle(const TestFunction &f,
                                        const std::vector<double> &partial_pos,
                                        const std::vector<int> &active_dims) const {
  if (partial_pos.size() != active_dims.size()) {
    throw std::invalid_argument(
        "size mismatch between partial position and active dimensions");
  }

  ensure_finite_vector(partial_pos, "partial particle position");
  validate_active_dims(active_dims, eval_scratch.size());

  // Only the active coordinates are replaced, the rest keeps the current context value.
  for (size_t i = 0; i < active_dims.size(); ++i) {
    eval_scratch[active_dims[i]] = partial_pos[i];
  }

  double fitness = sanitize_fitness(f.value(eval_scratch));

  // Restore the buffer so the next evaluation starts from the unchanged full context.
  for (size_t i = 0; i < active_dims.size(); ++i) {
    eval_scratch[active_dims[i]] = full_vector[active_dims[i]];
  }

  return fitness;
}

void ContextVector::set_full_vector(const std::vector<double> &new_vector,
                                    double fitness) {
  if (new_vector.size() != full_vector.size()) {
    throw std::invalid_argument(
        "context vector size mismatch in set_full_vector");
  }

  ensure_finite_vector(new_vector, "context vector");
  ensure_finite_value(fitness, "context fitness");

  full_vector = new_vector;
  // Keep the buffer aligned with the latest accepted full context.
  eval_scratch = new_vector;
  best_fitness = fitness;
}

const std::vector<double> &ContextVector::get_full_vector() const {
  return full_vector;
}
double ContextVector::get_best_fitness() const { return best_fitness; }
