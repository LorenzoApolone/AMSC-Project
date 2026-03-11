/**
 * @file dpso_serial.cpp
 * @brief Serial DMS-PSO-HS (Zhao et al., 2011).
 *
 * Single-thread version of the algorithm. The swarm is divided into
 * fixed-size sub-swarms; every regrouping_period iterations the particles
 * are globally reshuffled to promote diversity.
 * After each PSO update, a Harmony Search phase is applied to each
 * sub-swarm to intensify the local search.
 */

#include "interfaces.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include "methods_dpso.hpp"
#include "particle.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

/**
 * @brief Generates a random double within [min, max).
 *
 * @param min Minimum bound.
 * @param max Maximum bound.
 * @param gen PRNG instance.
 * @return Random double.
 */
static double random_double_serial(double min, double max, std::mt19937 &gen) {
  std::uniform_real_distribution<> dis(min, max);
  return dis(gen);
}

/**
 * @brief Generates a random integer within [min, max].
 *
 * @param min Minimum integer (inclusive).
 * @param max Maximum integer (inclusive).
 * @param gen PRNG instance.
 * @return Random integer.
 */
static int random_int_serial(int min, int max, std::mt19937 &gen) {
  std::uniform_int_distribution<> dis(min, max);
  return dis(gen);
}

/**
 * @brief Computes the squared Euclidean distance between two vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The squared Euclidean distance.
 */
static double euclidean_dist_squared_serial(const std::vector<double> &v1,
                                            const std::vector<double> &v2) {
  double sum = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    double diff = v1[i] - v2[i];
    sum += diff * diff;
  }
  return sum;
}

/**
 * @brief Applies Harmony Search to a sub-swarm.
 *
 * @param swarm Reference to the swarm.
 * @param start_idx Start index of the sub-swarm.
 * @param end_idx End index of the sub-swarm (exclusive).
 * @param f The test function to optimize.
 * @param gen PRNG instance.
 * @param lower_bound Vector of lower bounds.
 * @param upper_bound Vector of upper bounds.
 * @param current_iter Current iteration number.
 * @param max_iter Maximum number of iterations.
 * @param params Algorithm parameters containing HS configurations.
 * @param new_harmony Reusable buffer vector for the generated local harmony.
 */
static void apply_harmony_search_serial(
    std::vector<Particle> &swarm, int start_idx, int end_idx,
    const TestFunction &f, std::mt19937 &gen,
    const std::vector<double> &lower_bound,
    const std::vector<double> &upper_bound, int current_iter, int max_iter,
    const DPSOParameters &params, std::vector<double> &new_harmony) {
  int dim = lower_bound.size();
  int sub_swarm_size = end_idx - start_idx;
  if (sub_swarm_size <= 0)
    return;

  double PAR = params.par_min +
               ((params.par_max - params.par_min) / max_iter) * current_iter;

  double iter_ratio = (double)current_iter / max_iter;

  for (int d = 0; d < dim; ++d) {
    double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
    double bw_min = 0.0001;
    // Optimization: avoid slow exp/log combination, use std::pow directly
    double bw = bw_max * std::pow(bw_min / bw_max, iter_ratio);

    if (random_double_serial(0.0, 1.0, gen) < params.hmcr) {
      int idx = start_idx + random_int_serial(0, sub_swarm_size - 1, gen);
      new_harmony[d] = swarm[idx].best_position[d];
      if (random_double_serial(0.0, 1.0, gen) < PAR)
        new_harmony[d] += random_double_serial(-1.0, 1.0, gen) * bw;
    } else {
      new_harmony[d] =
          random_double_serial(lower_bound[d], upper_bound[d], gen);
    }
    new_harmony[d] =
        std::max(lower_bound[d], std::min(upper_bound[d], new_harmony[d]));
  }

  double new_val = f.value(new_harmony);
  int nearest_idx = -1;
  double min_dist = std::numeric_limits<double>::max();
  for (int i = start_idx; i < end_idx; ++i) {
    double d_sq =
        euclidean_dist_squared_serial(new_harmony, swarm[i].best_position);
    if (d_sq < min_dist) {
      min_dist = d_sq;
      nearest_idx = i;
    }
  }
  if (nearest_idx != -1 && new_val < swarm[nearest_idx].best_value) {
    swarm[nearest_idx].best_position = new_harmony;
    swarm[nearest_idx].best_value = new_val;
  }
}

/**
 * @brief Executes one local best (lbest) PSO + HS iteration on a sub-swarm.
 *
 * @param swarm Reference to the swarm.
 * @param start Start index of the sub-swarm.
 * @param end End index of the sub-swarm (exclusive).
 * @param dim Number of dimensions of the search space.
 * @param lb Vector of lower bounds.
 * @param ub Vector of upper bounds.
 * @param v_max Vector of maximum velocities.
 * @param params Algorithm parameters.
 * @param f The test function to optimize.
 * @param gen PRNG instance.
 * @param iter Current iteration number.
 * @param max_iter Maximum number of iterations.
 * @param hs_buffer Reusable buffer for harmony search.
 */
static void process_sub_swarm_serial(
    std::vector<Particle> &swarm, int start, int end, unsigned int dim,
    const std::vector<double> &lb, const std::vector<double> &ub,
    const std::vector<double> &v_max, const DPSOParameters &params,
    const TestFunction &f, std::mt19937 &gen, int iter, int max_iter,
    std::vector<double> &hs_buffer) {
  int lbest_idx = -1;
  double lbest_val = std::numeric_limits<double>::max();
  for (int i = start; i < end; ++i) {
    if (swarm[i].best_value < lbest_val) {
      lbest_val = swarm[i].best_value;
      lbest_idx = i;
    }
  }
  if (lbest_idx == -1)
    return;
  const std::vector<double> &lbest_pos = swarm[lbest_idx].best_position;

  for (int i = start; i < end; ++i) {
    Particle &p = swarm[i];
    for (unsigned int d = 0; d < dim; ++d) {
      double r1 = random_double_serial(0.0, 1.0, gen);
      double r2 = random_double_serial(0.0, 1.0, gen);
      p.velocity[d] = params.w * p.velocity[d] +
                      params.c1 * r1 * (p.best_position[d] - p.position[d]) +
                      params.c2 * r2 * (lbest_pos[d] - p.position[d]);
      p.velocity[d] = std::max(-v_max[d], std::min(v_max[d], p.velocity[d]));

      p.position[d] += p.velocity[d];

      // Boundary treatment: Clamp to domain boundaries and zero velocity
      if (p.position[d] < lb[d]) {
        p.position[d] = lb[d];
        p.velocity[d] = 0.0;
      } else if (p.position[d] > ub[d]) {
        p.position[d] = ub[d];
        p.velocity[d] = 0.0;
      }
    }

    p.current_value = f.value(p.position);
    if (p.current_value < p.best_value) {
      p.best_value = p.current_value;
      p.best_position = p.position;
    }
  }
  apply_harmony_search_serial(swarm, start, end, f, gen, lb, ub, iter, max_iter,
                              params, hs_buffer);
}

/**
 * @brief Shuffles the particles in the swarm to regroup sub-swarms.
 *
 * @param swarm Reference to the swarm.
 * @param g PRNG instance.
 */
static void regroup_particles_serial(std::vector<Particle> &swarm,
                                     std::mt19937 &g) {
  std::shuffle(swarm.begin(), swarm.end(), g);
}

/**
 * @brief Executes the Serial DMS-PSO-HS algorithm.
 *
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles in the swarm.
 * @param max_iter Maximum number of iterations.
 * @param params Algorithm parameters.
 * @param convergence_tol The tolerance used by the stopping criteria manager.
 * @return OutputObject containing the best position, value, and execution
 * metrics.
 */
OutputObject dpso_serial(const TestFunction &f, unsigned int dim,
                         unsigned int n_points_total, int max_iter,
                         const DPSOParameters &params, double convergence_tol) {
  StoppingCriteriaManager stop_manager(max_iter, 2000, convergence_tol, 1e-3);

  if (n_points_total < (unsigned int)params.sub_swarm_size) {
    std::cerr << "Error: total particles (" << n_points_total
              << ") less than sub-swarm size (" << params.sub_swarm_size
              << ").\n";
    return OutputObject(f.get_name(), dim, n_points_total, {},
                        f.get_true_solution(), 0.0, {}, 1, 0.0, 0,
                        stop_manager);
  }

  const auto &domain = f.get_domain();
  std::vector<double> lb(dim, domain.first);
  std::vector<double> ub(dim, domain.second);
  std::vector<double> v_max(dim);
  for (unsigned int d = 0; d < dim; ++d) {
    if (lb[d] > ub[d]) {
      std::swap(lb[d], ub[d]);
    }
  }
  for (unsigned int d = 0; d < dim; ++d)
    v_max[d] = 0.2 * (ub[d] - lb[d]);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::vector<Particle> swarm;
  swarm.reserve(n_points_total);
  for (unsigned int i = 0; i < n_points_total; ++i) {
    Particle p(dim);
    for (unsigned int d = 0; d < dim; ++d) {
      p.position[d] = random_double_serial(lb[d], ub[d], gen);
      p.velocity[d] = random_double_serial(-v_max[d], v_max[d], gen);
      p.best_position[d] = p.position[d];
    }
    p.current_value = f.value(p.position);
    p.best_value = p.current_value;
    swarm.push_back(p);
  }

  OutputObject results(f.get_name(), dim, n_points_total, {},
                       f.get_true_solution(), 0.0, {}, 1, 0.0, 0, stop_manager);
  results.x_best.resize(dim);
  double global_best_val = std::numeric_limits<double>::max();
  auto start_time = std::chrono::high_resolution_clock::now();
  int iter = 0;
  std::vector<double> hs_buffer(dim);

  while (true) {
    if (iter > 0 && iter % params.regrouping_period == 0)
      regroup_particles_serial(swarm, gen);

    int num_sub_swarms = swarm.size() / params.sub_swarm_size;
    int remainder = swarm.size() % params.sub_swarm_size;

    for (int s = 0; s < num_sub_swarms; ++s) {
      int start = s * params.sub_swarm_size;
      int end = std::min(start + params.sub_swarm_size, (int)swarm.size());
      process_sub_swarm_serial(swarm, start, end, dim, lb, ub, v_max, params, f,
                               gen, iter, max_iter, hs_buffer);
    }
    if (remainder > 0) {
      int start = num_sub_swarms * params.sub_swarm_size;
      process_sub_swarm_serial(swarm, start, (int)swarm.size(), dim, lb, ub,
                               v_max, params, f, gen, iter, max_iter,
                               hs_buffer);
    }

    double current_global_min = std::numeric_limits<double>::max();
    int best_idx = -1;
    for (int i = 0; i < (int)swarm.size(); ++i) {
      if (swarm[i].best_value < current_global_min) {
        current_global_min = swarm[i].best_value;
        best_idx = i;
      }
    }

    std::vector<double> global_best_pos(dim);
    if (best_idx != -1)
      global_best_pos = swarm[best_idx].best_position;

    double sum_dist = 0.0;
    for (const auto &p : swarm)
      sum_dist +=
          std::sqrt(euclidean_dist_squared_serial(p.position, global_best_pos));
    double avg_dist = swarm.empty() ? 0.0 : sum_dist / swarm.size();

    results.conv_history.push_back(current_global_min);
    global_best_val = current_global_min;

    stop_manager.increment_iterations();
    iter++;

    if (stop_manager.should_stop(global_best_val, avg_dist))
      break;
  }

  double best_val = std::numeric_limits<double>::max();
  int best_idx_final = -1;
  for (int i = 0; i < (int)swarm.size(); ++i) {
    if (swarm[i].best_value < best_val) {
      best_val = swarm[i].best_value;
      best_idx_final = i;
    }
  }
  if (best_idx_final != -1) {
    results.x_best = swarm[best_idx_final].best_position;
    results.f_val = best_val;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  results.execution_time =
      std::chrono::duration<double>(end_time - start_time).count();
  results.iterations = iter;
  return results;
}
