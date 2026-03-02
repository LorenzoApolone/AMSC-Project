/**
 * @file SubSwarm.hpp
 * @brief Defines the SubSwarm class and CPSOParticle struct used in Cooperative
 * Particle Swarm Optimization
 */

#pragma once

#include "../interfaces.hpp"
#include "ContextVector.hpp"
#include <limits>
#include <random>
#include <vector>

/**
 * @struct CPSOParticle
 * @brief Represents a single particle for the CPSO algorithm
 */
struct CPSOParticle {

  // Position of the particle in the sub-swarm
  std::vector<double> position; 

  // Velocity of the particle
  std::vector<double> velocity; 

  // Best position found by this particle
  std::vector<double> best_position;
  
  // Fitness value of the current position
  double current_value; 

  // Fitness value of the best position found
  double best_value; 

  /**
   * @brief Constructs a CPSOParticle with a given dimensionality
   * @param dim The number of dimensions for this particle
   */
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
 * @brief Manages a swarm of particles operating on a specific subset of the problem
 */
class SubSwarm {
private:

  // The set of particles in this sub-swarm
  std::vector<CPSOParticle> particles; 

  // The indices of dimensions that this sub-swarm optimizes
  std::vector<int> active_dims;  

  // The best position found by the entire sub-swarm in the active dimensions
  std::vector<double> gbest_pos; 

  // The fitness value of the sub-swarm's global best position
  double gbest_val; 

  // Network topology defining particle neighbors for local-best calculation
  std::vector<std::vector<int>> adjacency_list; 

  // Lower bound of the search space
  double bounds_lower; 
  
  // Upper bound of the search space
  double bounds_upper; 

public:

  /**
   * @brief Constructs a new SubSwarm
   * @param num_particles Number of particles in the sub-swarm
   * @param active_dimensions The global dimension indices this sub-swarm will
   * optimize
   * @param lower_bound The lower boundary of the space
   * @param upper_bound The upper boundary of the space
   * @param adj_list Network topology for social interactions
   */
  SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
           double lower_bound, double upper_bound,
           const std::vector<std::vector<int>> &adj_list = {});

  /**
   * @brief Initializes the particles' positions and velocities randomly within
   * bounds, and evaluates their initial fitness.
   * @param gen Random number generator.
   * @param ctx ContextVector used to provide fixed values for inactive
   * dimensions during evaluation.
   * @param f The test function to be optimized.
   */
  void initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Updates the dimensions optimized by this sub-swarm, updating
   * internal states to match new dimensions.
   * @param new_dims The new set of dimension indices to optimize.
   * @param ctx The ContextVector to extract current values for the new
   * dimensions.
   * @param gen Random number generator for injecting small initial velocities.
   */
  void update_active_dims(const std::vector<int> &new_dims,
                          const ContextVector &ctx, std::mt19937 &gen);

  /**
   * @brief Injects random velocities into particles to prevent stagnation.
   * @param gen Random number generator.
   */
  void inject_velocities(std::mt19937 &gen);

  /**
   * @brief Updates the velocity and position of each particle using PSO
   * equations.
   * @param w Inertia weight.
   * @param c1 Cognitive learning factor.
   * @param c2 Social learning factor.
   * @param gen Random number generator.
   */
  void update_velocities_and_positions(double w, double c1, double c2,
                                       std::mt19937 &gen);

  /**
   * @brief Evaluates the fitness of all particles and updates their personal
   * bests and the sub-swarm's global best.
   * @param ctx ContextVector used to combine active and inactive dimensions for
   * full evaluation.
   * @param f The test function to be optimized.
   */
  void evaluate_and_update(ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Gets the global best position found by this sub-swarm.
   * @return A constant reference to the global best position vector
   * (sub-dimensional).
   */
  const std::vector<double> &get_gbest_pos() const;

  /**
   * @brief Gets the fitness value of the sub-swarm's global best position.
   * @return The global best fitness value.
   */
  double get_gbest_val() const;

  /**
   * @brief Gets the list of active dimensions optimized by this sub-swarm.
   * @return A constant reference to the active dimension indices vector.
   */
  const std::vector<int> &get_active_dims() const;

  /**
   * @brief Gets the collection of particles in the sub-swarm.
   * @return A constant reference to the vector of CPSOParticle objects.
   */
  const std::vector<CPSOParticle> &get_particles() const;
};
