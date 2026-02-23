#pragma once

#include "../interfaces.hpp"
#include "ContextVector.hpp"
#include <limits>
#include <random>
#include <vector>

/**
 * @struct CPSOParticle
 * Lightweight structure differentiated from the original Particle.hpp.
 * Maintains position/velocity only for the assigned subset of dimensions.
 */
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

/**
 * @class SubSwarm
 * Abstraction of the sub-swarm in CPSO. Manages and updates an array of
 * CPSOParticles on a specific subset of dimensions (active_dims).
 */
class SubSwarm {
private:
  std::vector<CPSOParticle> particles;
  std::vector<int> active_dims; // Optimized global indices

  std::vector<double> gbest_pos; // Class best memory (local to the subset)
  double gbest_val;

  std::vector<std::vector<int>> adjacency_list; // Scale-Free or other Topology

  double bounds_lower;
  double bounds_upper;

public:
  SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
           double lower_bound, double upper_bound,
           const std::vector<std::vector<int>> &adj_list = {});

  /**
   * @brief Uniformly initializes the positions and velocities of the particles
   * and performs the first evaluation using the Context Vector.
   */
  void initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Updates velocities and positions based on the standard PSO equation.
   * Uses the provided network topology to calculate the lBest (Local Best) of
   * each particle.
   */
  void update_velocities_and_positions(double w, double c1, double c2,
                                       std::mt19937 &gen);

  /**
   * @brief Formally evaluates the points via the Context Vector to
   * calculate the new fitness, thus updating the sub-swarm's local
   * pBest and gBest logic.
   * @return Returns the number of fevals executed to track the computational
   * budget.
   */
  int evaluate_and_update(ContextVector &ctx, const TestFunction &f);

  // Getters
  const std::vector<double> &get_gbest_pos() const;
  double get_gbest_val() const;
  const std::vector<int> &get_active_dims() const;
  const std::vector<CPSOParticle> &get_particles() const;
};
