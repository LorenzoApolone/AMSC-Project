/**
 * @file SubSwarm.hpp
 * @brief Defines the CPSO sub-swarm and its particle state.
 */
#pragma once

#include "../interfaces/interfaces.hpp"
#include "ContextVector.hpp"
#include <limits>
#include <random>
#include <vector>

/**
 * @class SubSwarm
 * @brief Represents one subswarm inside CPSO.
 */
class SubSwarm {
private:
  int num_particles;
  std::vector<double> positions;
  std::vector<double> velocities;
  std::vector<double> best_positions;
  std::vector<double> current_values;
  std::vector<double> best_values;
  std::vector<int> active_dims;
  std::vector<double> gbest_pos;
  double gbest_val;
  int gbest_particle_idx;
  bool ignore_gbest_this_iter;
  std::vector<std::vector<int>> adjacency_list;
  double bounds_lower;
  double bounds_upper;

public:
  /**
   * @brief Builds a fully allocated subswarm.
   * @param num_particles Number of particles stored in the subswarm.
   * @param active_dimensions Global dimensions assigned to the subswarm.
   * @param lower_bound Lower domain bound for every active dimension.
   * @param upper_bound Upper domain bound for every active dimension.
   * @param adj_list Optional neighborhood graph over the local particles.
   */
  SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
           double lower_bound, double upper_bound,
           const std::vector<std::vector<int>> &adj_list = {});

  SubSwarm(const SubSwarm &) = default;
  SubSwarm &operator=(const SubSwarm &) = default;
  SubSwarm(SubSwarm &&) noexcept = default;
  SubSwarm &operator=(SubSwarm &&) noexcept = default;
  ~SubSwarm() = default;

  /**
   * @brief Initializes particle positions, velocities and local best memories.
   * @param gen Random generator.
   * @param ctx Cooperative context vector used to evaluate particles.
   * @param f Benchmark function optimized by the subswarm.
   */
  void initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Evaluates current and personal-best particles after a context change.
   * @param ctx Cooperative context vector used to evaluate particles.
   * @param f Benchmark function optimized by the subswarm.
   */
  void recalculate_fitness(ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Updates the dimensions handled by the subswarm after a shuffle.
   * @param new_dims New global dimensions assigned to the subswarm.
   * @param ctx Current cooperative context used to seed reshuffled coordinates.
   * @param gen Random generator used when new coordinates must be sampled.
   * @param is_owned true when the local rank owns the subswarm storage after the shuffle.
   */
  void update_active_dims(const std::vector<int> &new_dims,
                          const ContextVector &ctx, std::mt19937 &gen,
                          bool is_owned = true);

  /**
   * @brief Perturbs particle velocities to escape prolonged stagnation.
   * @param gen Random generator used to sample the perturbation.
   * @param hard_reset true to also reset personal-best memories.
   */
  void inject_velocities(std::mt19937 &gen, bool hard_reset = false);

  /**
   * @brief Disables attraction toward the current subswarm best for one iteration.
   */
  void reset_gbest_attraction();

  /**
   * @brief Applies one PSO velocity/position update to the local particles.
   * @param w Inertia coefficient.
   * @param c1 Cognitive coefficient.
   * @param c2 Social coefficient.
   * @param gen Random generator used to sample the stochastic factors.
   * @param progress_ratio Normalized run progress used to adapt the velocity clamp.
   */
  void update_velocities_and_positions(double w, double c1, double c2,
                                       std::mt19937 &gen,
                                       double progress_ratio = 0.5);

  /**
   * @brief Evaluates the updated particles and refreshes local/global best memories.
   * @param ctx Cooperative context used to evaluate particles.
   * @param f Benchmark function optimized by the subswarm.
   */
  void evaluate_and_update(ContextVector &ctx, const TestFunction &f);

  /**
   * @brief Returns the best position currently known by the subswarm.
   * @return The subswarm best position in local coordinates.
   */
  const std::vector<double> &get_gbest_pos() const;

  /**
   * @brief Returns the fitness of the current subswarm best.
   * @return Best fitness currently stored by the subswarm.
   */
  double get_gbest_val() const;

  /**
   * @brief Returns the global dimensions currently assigned to the subswarm.
   * @return The active global dimensions handled by the subswarm.
   */
  const std::vector<int> &get_active_dims() const;

  /**
   * @brief Returns the number of locally stored particles.
   * @return Number of particles stored in the subswarm.
   */
  int get_num_particles() const;

  /**
   * @brief Returns the flattened particle-position buffer.
   * @return Flat storage containing the current particle positions.
   */
  const std::vector<double> &get_positions() const;
};
