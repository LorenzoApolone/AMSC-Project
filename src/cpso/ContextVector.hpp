#pragma once

#include "../interfaces.hpp"
#include <vector>

class ContextVector {
private:
  std::vector<double> full_vector;
  double best_fitness;

public:
  ContextVector(int total_dims);

  const std::vector<double> &get_full_vector() const;
  double get_best_fitness() const;
  double evaluate_particle(const TestFunction &f,
                           const std::vector<double> &partial_pos,
                           const std::vector<int> &active_dims) const;

  void update(const std::vector<double> &partial_best_pos,
              const std::vector<int> &active_dims, double new_fitness);

  void set_full_vector(const std::vector<double> &new_vector, double fitness);
};
