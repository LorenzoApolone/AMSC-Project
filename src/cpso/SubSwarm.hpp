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
 * @class SubSwarm
 * @brief Manages a swarm of particles operating on a specific subset of the problem
 */
class SubSwarm {
private:

  // Number of particles in this sub-swarm
  int num_particles;

  // Positions of all particles in the sub-swarm
  std::vector<double> positions; 

  // Velocities of all particles in the sub-swarm
  std::vector<double> velocities; 

  // Best positions found by all particles
  std::vector<double> best_positions;
  
  // Fitness values of current positions
  std::vector<double> current_values; 

  // Fitness values of best positions found
  std::vector<double> best_values;
  // The indices of dimensions that this sub-swarm optimizes
  std::vector<int> active_dims;  

  // The best position found by the entire sub-swarm in the active dimensions
  std::vector<double> gbest_pos; 

  // The fitness value of the sub-swarm's global best position
  double gbest_val; 

  // Index of the global best particle
  int gbest_particle_idx;

  // Flag to temporarily disable gbest social attraction during stagnation resets
  bool ignore_gbest_this_iter;

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

  // Default copy constructor
  SubSwarm(const SubSwarm&) = default;

  // Default copy assignment operator
  SubSwarm& operator=(const SubSwarm&) = default;

  // Move Constructor
  SubSwarm(SubSwarm&&) noexcept = default;

  // Move assignment operator
  SubSwarm& operator=(SubSwarm&&) noexcept = default;

  // Default destructor
  ~SubSwarm() = default;

  /**
   * @brief Initializes the particles' positions and velocities randomly within
   * bounds and evaluates their initial fitness
   * @param gen Random number generator
   * @param ctx ContextVector used to provide fixed values for inactive dimensions
   * @param f The test function to be optimized
   */
  void initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Recalculates the historical best fitness of the particles and the sub-swarm's global best
   * based on the current ContextVector
   * @param ctx The ContextVector containing the current state of inactive dimensions
   * @param f The test function to be optimized
   */
  void recalculate_fitness(ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Updates the dimensions optimized by this sub-swarm, updating
   * internal states to match new dimensions
   * @param new_dims The new set of dimension indices to optimize
   * @param ctx The ContextVector to extract current values for the new dimensions
   * @param gen Random number generator for injecting small initial velocities
   * @param is_owned Whether the current MPI rank owns this sub-swarm 
   */
  void update_active_dims(const std::vector<int> &new_dims,
                          const ContextVector &ctx, std::mt19937 &gen, bool is_owned = true);

  /**
   * @brief Injects random velocities into particles to prevent stagnation
   * @param gen Random number generator
   * @param hard_reset Whether to reset velocities to zero
   */
  void inject_velocities(std::mt19937 &gen, bool hard_reset = false);

  /**
   * @brief Sets a flag to ignore global best attraction for one iter,
   * used during a hard jump/injection to avoid collapse
   */
  void reset_gbest_attraction();

  /**
   * @brief Updates the velocity and position of each particle using PSO
   * equations
   * @param w Inertia weight
   * @param c1 Cognitive learning factor
   * @param c2 Social learning factor
   * @param gen Random number generator
   */
  void update_velocities_and_positions(double w, double c1, double c2,
                                       std::mt19937 &gen, double progress_ratio = 0.5);

  /**
   * @brief Evaluates the fitness of all particles and updates their personal
   * bests and the sub-swarm's global best
   * @param ctx ContextVector used to combine active and inactive dimensions for
   * full evaluation
   * @param f The test function to be optimized
   */
  void evaluate_and_update(ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Gets the global best position found by this sub-swarm 
   * @return A constant reference to the global best position vector
   * (sub-dimensional)
   */
  const std::vector<double> &get_gbest_pos() const;

  /**
   * @brief Gets the fitness value of the sub-swarm's global best position
   * @return The global best fitness value
   */
  double get_gbest_val() const;

  /**
   * @brief Gets the list of active dimensions optimized by this sub-swarm
   * @return A constant reference to the active dimension indices vector
   */
  const std::vector<int> &get_active_dims() const;

  /**
   * @brief Gets the number of particles in the sub-swarm
   * @return The number of particles
   */
  int get_num_particles() const;

  /**
   * @brief Gets the positions of all particles in the sub-swarm
   * @return A constant reference to the flat vector of particle positions
   */
  const std::vector<double> &get_positions() const;
};
