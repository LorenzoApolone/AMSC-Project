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
#include <stdexcept>

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
 * @param pos Swarm positions (SoA).
 * @param pbest_pos Swarm personal best positions (SoA).
 * @param pbest_val Swarm personal best values (SoA).
 * @param start_idx Start index of the sub-swarm.
 * @param end_idx End index of the sub-swarm (exclusive).
 * @param dim Number of dimensions of the search space.
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
    std::vector<double> &pos, std::vector<double> &current_val,
    std::vector<double> &pbest_pos, std::vector<double> &pbest_val,
    int start_idx, int end_idx, int dim, const TestFunction &f,
    std::mt19937 &gen, const std::vector<double> &lower_bound,
    const std::vector<double> &upper_bound, int current_iter, int max_iter,
    const DPSOParameters &params, std::vector<double> &new_harmony) {
  int sub_swarm_size = end_idx - start_idx;
  if (sub_swarm_size <= 0)
    return;

  double PAR = params.par_min +
               ((params.par_max - params.par_min) / max_iter) * current_iter;

  double iter_ratio = (double)current_iter / max_iter;

  for (int d = 0; d < dim; ++d) {
    double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
    double bw_min = 0.0001;
    double bw = bw_max * std::pow(bw_min / bw_max, iter_ratio);

    if (random_double_serial(0.0, 1.0, gen) < params.hmcr) {
      int p_idx = start_idx + random_int_serial(0, sub_swarm_size - 1, gen);
      new_harmony[d] = pbest_pos[p_idx * dim + d];
      if (random_double_serial(0.0, 1.0, gen) < PAR)
        new_harmony[d] += random_double_serial(-1.0, 1.0, gen) * bw;
    } else {
      new_harmony[d] = random_double_serial(lower_bound[d], upper_bound[d], gen);
    }
    new_harmony[d] =
        std::max(lower_bound[d], std::min(upper_bound[d], new_harmony[d]));
  }

  double new_val = f.value(new_harmony);
  int nearest_idx = -1;
  double min_dist = std::numeric_limits<double>::max();
  for (int i = start_idx; i < end_idx; ++i) {
    double d_sq = 0.0;
    for (int d = 0; d < dim; ++d) {
      double diff = new_harmony[d] - pbest_pos[i * dim + d];
      d_sq += diff * diff;
    }
    if (d_sq < min_dist) {
      min_dist = d_sq;
      nearest_idx = i;
    }
  }
  if (nearest_idx != -1 && new_val < pbest_val[nearest_idx]) {
    for (int d = 0; d < dim; ++d) {
      pbest_pos[nearest_idx * dim + d] = new_harmony[d];
      pos[nearest_idx * dim + d] = new_harmony[d];
    }
    pbest_val[nearest_idx] = new_val;
    current_val[nearest_idx] = new_val;
  }
}

/**
 * @brief Executes one local best (lbest) PSO + HS iteration on a sub-swarm.
 *
 * @param pos Swarm positions (SoA).
 * @param vel Swarm velocities (SoA).
 * @param pbest_pos Swarm personal best positions (SoA).
 * @param pbest_val Swarm personal best values (SoA).
 * @param current_val Swarm current values (SoA).
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
    std::vector<double> &pos, std::vector<double> &vel,
    std::vector<double> &pbest_pos, std::vector<double> &pbest_val,
    std::vector<double> &current_val, int start, int end, unsigned int dim,
    const std::vector<double> &lb, const std::vector<double> &ub,
    const std::vector<double> &v_max, const DPSOParameters &params,
    const TestFunction &f, std::mt19937 &gen, int iter, int max_iter,
    std::vector<double> &hs_buffer) {
  int lbest_idx = -1;
  double lbest_val = std::numeric_limits<double>::max();
  for (int i = start; i < end; ++i) {
    if (pbest_val[i] < lbest_val) {
      lbest_val = pbest_val[i];
      lbest_idx = i;
    }
  }
  if (lbest_idx == -1)
    return;

  for (int i = start; i < end; ++i) {
    for (unsigned int d = 0; d < dim; ++d) {
      double r1 = random_double_serial(0.0, 1.0, gen);
      double r2 = random_double_serial(0.0, 1.0, gen);
      int idx = i * dim + d;
      int lb_idx = lbest_idx * dim + d;

      vel[idx] = params.w * vel[idx] +
                 params.c1 * r1 * (pbest_pos[idx] - pos[idx]) +
                 params.c2 * r2 * (pbest_pos[lb_idx] - pos[idx]);
      vel[idx] = std::max(-v_max[d], std::min(v_max[d], vel[idx]));

      pos[idx] += vel[idx];

      if (pos[idx] < lb[d]) {
        pos[idx] = lb[d];
        vel[idx] = 0.0;
      } else if (pos[idx] > ub[d]) {
        pos[idx] = ub[d];
        vel[idx] = 0.0;
      }
    }

    std::vector<double> p_pos(dim);
    for (unsigned int d = 0; d < dim; ++d)
      p_pos[d] = pos[i * dim + d];

    current_val[i] = f.value(p_pos);
    if (current_val[i] < pbest_val[i]) {
      pbest_val[i] = current_val[i];
      for (unsigned int d = 0; d < dim; ++d)
        pbest_pos[i * dim + d] = pos[i * dim + d];
    }
  }
  apply_harmony_search_serial(pos, current_val, pbest_pos, pbest_val, start, end, dim, f, gen, lb,
                              ub, iter, max_iter, params, hs_buffer);
}

/**
 * @brief Shuffles the particles in the swarm to regroup sub-swarms in SoA.
 *
 * @param pos Swarm positions.
 * @param vel Swarm velocities.
 * @param pbest_pos Swarm personal best positions.
 * @param pbest_val Swarm personal best values.
 * @param current_val Swarm current values.
 * @param n_points Total number of particles.
 * @param dim Dimensions.
 * @param g PRNG instance.
 */
static void regroup_particles_serial(std::vector<double> &pos,
                                     std::vector<double> &vel,
                                     std::vector<double> &pbest_pos,
                                     std::vector<double> &pbest_val,
                                     std::vector<double> &current_val,
                                     int n_points, int dim, std::mt19937 &g) {
  std::vector<int> indices(n_points);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), g);

  std::vector<double> new_pos(n_points * dim);
  std::vector<double> new_vel(n_points * dim);
  std::vector<double> new_pbest_pos(n_points * dim);
  std::vector<double> new_pbest_val(n_points);
  std::vector<double> new_current_val(n_points);

  for (int i = 0; i < n_points; ++i) {
    int old_idx = indices[i];
    std::copy_n(pos.begin() + old_idx * dim, dim, new_pos.begin() + i * dim);
    std::copy_n(vel.begin() + old_idx * dim, dim, new_vel.begin() + i * dim);
    std::copy_n(pbest_pos.begin() + old_idx * dim, dim,
                new_pbest_pos.begin() + i * dim);
    new_pbest_val[i] = pbest_val[old_idx];
    new_current_val[i] = current_val[old_idx];
  }

  pos = std::move(new_pos);
  vel = std::move(new_vel);
  pbest_pos = std::move(new_pbest_pos);
  pbest_val = std::move(new_pbest_val);
  current_val = std::move(new_current_val);
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
                         const DPSOParameters &params, double convergence_tol, unsigned int seed) {
  // Stopping criteria aligned with CPSO:
  //   iterations_stagnation = max(100, max_iter / 4)
  //   stagnation_tol        = 1e-8
  //   diversity_tol         = 1e-6 (default)
  int scaled_stagnation = std::max(100, max_iter / 4);
  StoppingCriteriaManager stop_manager(max_iter, scaled_stagnation, 1e-8, 1e-6);

  if (n_points_total < (unsigned int)params.sub_swarm_size) {
    throw std::invalid_argument("Error: total particles (" + std::to_string(n_points_total) +
                                ") less than sub-swarm size (" + std::to_string(params.sub_swarm_size) +
                                ").");
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
  std::mt19937 gen(seed == 0 ? rd() : seed);

  std::vector<double> pos(n_points_total * dim);
  std::vector<double> vel(n_points_total * dim);
  std::vector<double> pbest_pos(n_points_total * dim);
  std::vector<double> pbest_val(n_points_total,
                                std::numeric_limits<double>::max());
  std::vector<double> current_val(n_points_total);

  for (unsigned int i = 0; i < n_points_total; ++i) {
    std::vector<double> p_pos(dim);
    for (unsigned int d = 0; d < dim; ++d) {
      pos[i * dim + d] = random_double_serial(lb[d], ub[d], gen);
      vel[i * dim + d] = random_double_serial(-v_max[d], v_max[d], gen);
      pbest_pos[i * dim + d] = pos[i * dim + d];
      p_pos[d] = pos[i * dim + d];
    }
    current_val[i] = f.value(p_pos);
    pbest_val[i] = current_val[i];
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
      regroup_particles_serial(pos, vel, pbest_pos, pbest_val, current_val,
                               n_points_total, dim, gen);

    int num_sub_swarms = n_points_total / params.sub_swarm_size;
    int remainder = n_points_total % params.sub_swarm_size;

    for (int s = 0; s < num_sub_swarms; ++s) {
      int start = s * params.sub_swarm_size;
      int end = start + params.sub_swarm_size;
      process_sub_swarm_serial(pos, vel, pbest_pos, pbest_val, current_val,
                               start, end, dim, lb, ub, v_max, params, f, gen,
                               iter, max_iter, hs_buffer);
    }
    if (remainder > 0) {
      int start = num_sub_swarms * params.sub_swarm_size;
      process_sub_swarm_serial(pos, vel, pbest_pos, pbest_val, current_val,
                               start, n_points_total, dim, lb, ub, v_max, params,
                               f, gen, iter, max_iter, hs_buffer);
    }

    double current_global_min = std::numeric_limits<double>::max();
    int best_idx = -1;
    for (int i = 0; i < (int)n_points_total; ++i) {
      if (pbest_val[i] < current_global_min) {
        current_global_min = pbest_val[i];
        best_idx = i;
      }
    }

    std::vector<double> global_best_pos(dim);
    if (best_idx != -1) {
      for (unsigned int d = 0; d < dim; ++d)
        global_best_pos[d] = pbest_pos[best_idx * dim + d];
    }

    double sum_dist = 0.0;
    for (unsigned int i = 0; i < n_points_total; ++i) {
      double d_sq = 0.0;
      for (unsigned int d = 0; d < dim; ++d) {
        double diff = pos[i * dim + d] - global_best_pos[d];
        d_sq += diff * diff;
      }
      sum_dist += std::sqrt(d_sq);
    }
    double avg_dist = n_points_total == 0 ? 0.0 : sum_dist / n_points_total;

    results.conv_history.push_back(current_global_min);
    global_best_val = current_global_min;

    stop_manager.increment_iterations();
    iter++;

    if (stop_manager.should_stop(global_best_val, avg_dist))
      break;
  }

  double best_val_final = std::numeric_limits<double>::max();
  int best_idx_final = -1;
  for (int i = 0; i < (int)n_points_total; ++i) {
    if (pbest_val[i] < best_val_final) {
      best_val_final = pbest_val[i];
      best_idx_final = i;
    }
  }
  if (best_idx_final != -1) {
    for (unsigned int d = 0; d < dim; ++d)
      results.x_best[d] = pbest_pos[best_idx_final * dim + d];
    results.f_val = best_val_final;
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  results.execution_time =
      std::chrono::duration<double>(end_time - start_time).count();
  results.iterations = iter;
  return results;
}
