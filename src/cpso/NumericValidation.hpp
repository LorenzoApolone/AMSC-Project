/**
 * @file NumericValidation.hpp
 * @brief Small helpers used to guard CPSO state against non-finite values.
 */
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Checks whether a scalar value is finite.
 */
inline bool is_finite_value(double value) { return std::isfinite(value); }

/**
 * @brief Throws when a scalar value is not finite.
 */
inline void ensure_finite_value(double value, const char *context) {
  if (!is_finite_value(value)) {
    throw std::runtime_error(std::string(context) + " must be finite");
  }
}

/**
 * @brief Throws when any element of a vector is not finite.
 */
inline void ensure_finite_vector(const std::vector<double> &values,
                                 const char *context) {
  for (double value : values) {
    if (!is_finite_value(value)) {
      throw std::runtime_error(std::string(context) +
                               " contains a non-finite value");
    }
  }
}

/**
 * @brief Replaces non-finite fitness values with positive infinity.
 */
inline double sanitize_fitness(double value) {
  return is_finite_value(value) ? value
                                : std::numeric_limits<double>::infinity();
}
