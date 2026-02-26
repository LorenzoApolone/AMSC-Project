#pragma once

#include "../interfaces.hpp"
#include "ContextVector.hpp"
#include <limits>
#include <random>
#include <vector>

struct CPSOParticle {
  std::vector<double> position;
  std::vector<double> velocity;
  std::vector<double> best_position;
  double current_value;
  double best_value;

  CPSOParticle(int dim) {
    position.resize(dim);
    velocity.resize(dim);
    best_position.resize(dim);
    current_value = std::numeric_limits<double>::infinity();
    best_value = std::numeric_limits<double>::infinity();
  }
};

class SubSwarm {
private:
  std::vector<CPSOParticle> particles;
  std::vector<int> active_dims;

  std::vector<double> gbest_pos;
  double gbest_val;

  std::vector<std::vector<int>> adjacency_list;

  double bounds_lower;
  double bounds_upper;

public:
  SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
           double lower_bound, double upper_bound,
           const std::vector<std::vector<int>> &adj_list = {});

  void initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f);

  void update_active_dims(const std::vector<int> &new_dims);

  void update_velocities_and_positions(double w, double c1, double c2,
                                       std::mt19937 &gen);

  int evaluate_and_update(ContextVector &ctx, const TestFunction &f);

  const std::vector<double> &get_gbest_pos() const;
  double get_gbest_val() const;
  const std::vector<int> &get_active_dims() const;
  const std::vector<CPSOParticle> &get_particles() const;
};
