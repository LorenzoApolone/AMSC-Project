/**
 * @file dpso.cpp
 * @brief Distributed DMS-PSO-HS via MPI (Zhao et al., 2011).
 *
 * Each rank manages a local sub-population. Every regrouping_period
 * iterations the global swarm is reshuffled across ranks via
 * MPI_Allgather + shuffle, ensuring information exchange and diversity.
 * Within each sub-swarm a Harmony Search phase is applied to the local
 * pbests, as described in the original paper.
 */

#include "interfaces.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include "methods_dpso.hpp"
#include "particle.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <random>
#include <vector>

/**
 * @brief Generates a random double within [min, max) using a rank-dependent
 * seed.
 *
 * Per-rank generator: different seed to ensure independent trajectories
 * across MPI processes during initialization and PSO updates.
 *
 * @param min Minimum bound.
 * @param max Maximum bound.
 * @param rank MPI rank of the process.
 * @return A random double.
 */
double random_double(double min, double max, int rank) {
  static std::mt19937 gen(rank * 10000 + 12345);
  std::uniform_real_distribution<> dis(min, max);
  return dis(gen);
}

/**
 * @brief Generates a random integer within [min, max] using a rank-dependent
 * seed.
 *
 * Avoids precision issues resulting from casting double variables.
 *
 * @param min Minimum integer (inclusive).
 * @param max Maximum integer (inclusive).
 * @param rank MPI rank of the process.
 * @return A random integer.
 */
int random_int(int min, int max, int rank) {
  static std::mt19937 gen(rank * 10000 + 54321);
  std::uniform_int_distribution<> dis(min, max);
  return dis(gen);
}

/**
 * @brief Computes the Euclidean distance between two vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The Euclidean distance.
 */
double euclidean_dist(const std::vector<double> &v1,
                      const std::vector<double> &v2) {
  double sum = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    double diff = v1[i] - v2[i];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

/**
 * HS phase on a sub-swarm [start_idx, end_idx).
 *
 * Generates a new harmony from the sub-swarm pbests with adaptive HMCR
 * and PAR (Eq. 20-21 of the paper). If the harmony beats the nearest
 * pbest (Euclidean distance), it replaces it.
 * PAR and bw increase/decrease linearly/exponentially with the iteration
 * to balance exploration in the early stages and intensification at the end.
 */
void apply_harmony_search(std::vector<Particle> &swarm, int start_idx,
                          int end_idx, const TestFunction &f, int rank,
                          const std::vector<double> &lower_bound,
                          const std::vector<double> &upper_bound,
                          int current_iter, int max_iter,
                          const DPSOParameters &params) {
  int dim = lower_bound.size();
  int sub_swarm_size = end_idx - start_idx;
  if (sub_swarm_size <= 0)
    return;

  double PAR = params.par_min +
               ((params.par_max - params.par_min) / max_iter) * current_iter;

  std::vector<double> new_harmony(dim);
  double iter_ratio = (double)current_iter / max_iter;

  for (int d = 0; d < dim; ++d) {
    double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
    double bw_min = 0.0001;
    double bw = bw_max * std::exp(std::log(bw_min / bw_max) * iter_ratio);

    if (random_double(0.0, 1.0, rank) < params.hmcr) {
      int idx = start_idx + random_int(0, sub_swarm_size - 1, rank);
      new_harmony[d] = swarm[idx].best_position[d];
      if (random_double(0.0, 1.0, rank) < PAR)
        new_harmony[d] += random_double(-1.0, 1.0, rank) * bw;
    } else {
      new_harmony[d] = random_double(lower_bound[d], upper_bound[d], rank);
    }
    new_harmony[d] =
        std::max(lower_bound[d], std::min(upper_bound[d], new_harmony[d]));
  }

  double new_val = f.value(new_harmony);
  int nearest_idx = -1;
  double min_dist = std::numeric_limits<double>::max();
  for (int i = start_idx; i < end_idx; ++i) {
    double d = euclidean_dist(new_harmony, swarm[i].best_position);
    if (d < min_dist) {
      min_dist = d;
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
 * Used for both complete sub-swarms and the possible remainder group.
 *
 * @param swarm Reference to the rank's local swarm.
 * @param start Start index of the sub-swarm.
 * @param end End index of the sub-swarm (exclusive).
 * @param dim Number of dimensions of the search space.
 * @param lb Vector of lower bounds.
 * @param ub Vector of upper bounds.
 * @param v_max Vector of maximum velocities.
 * @param w Inertia weight.
 * @param c1 Cognitive coefficient.
 * @param c2 Social coefficient.
 * @param f The test function to optimize.
 * @param rank MPI rank of the process.
 * @param iter Current iteration number.
 * @param max_iter Maximum number of iterations.
 */
static void process_sub_swarm(std::vector<Particle> &swarm, int start, int end,
                              unsigned int dim, const std::vector<double> &lb,
                              const std::vector<double> &ub,
                              const std::vector<double> &v_max,
                              const DPSOParameters &params,
                              const TestFunction &f, int rank, int iter,
                              int max_iter) {
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
  std::vector<double> lbest_pos = swarm[lbest_idx].best_position;

  for (int i = start; i < end; ++i) {
    Particle &p = swarm[i];
    for (unsigned int d = 0; d < dim; ++d) {
      double r1 = random_double(0.0, 1.0, rank);
      double r2 = random_double(0.0, 1.0, rank);
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
  apply_harmony_search(swarm, start, end, f, rank, lb, ub, iter, max_iter,
                       params);
}

/**
 * Global regrouping: gathers all particles with MPI_Allgather,
 * then each rank selects its own portion from a random permutation
 * with a rank-dependent seed, so that ranks obtain disjoint and
 * mutually different partitions.
 */
void regroup_particles(std::vector<Particle> &local_swarm, int dim, int rank,
                       int size) {
  int local_n = local_swarm.size();
  int p_data_size = 3 * dim + 2;

  // Optimization: avoid memory reallocation inside optimization loop
  static std::vector<double> send_buffer;
  send_buffer.clear();
  send_buffer.reserve(local_n * p_data_size);
  for (const auto &p : local_swarm) {
    send_buffer.insert(send_buffer.end(), p.position.begin(), p.position.end());
    send_buffer.insert(send_buffer.end(), p.velocity.begin(), p.velocity.end());
    send_buffer.insert(send_buffer.end(), p.best_position.begin(),
                       p.best_position.end());
    send_buffer.push_back(p.best_value);
    send_buffer.push_back(p.current_value);
  }

  static std::vector<double> recv_buffer;
  recv_buffer.resize(local_n * size * p_data_size);
  MPI_Allgather(send_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                recv_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                MPI_COMM_WORLD);

  static std::mt19937 g(19349663);
  static std::vector<int> indices;
  indices.resize(local_n * size);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), g);

  int global_idx_start = rank * local_n;
  for (int i = 0; i < local_n; ++i) {
    int picked = indices[global_idx_start + i];
    int base = picked * p_data_size;
    int off = 0;
    std::copy(recv_buffer.begin() + base + off,
              recv_buffer.begin() + base + off + dim,
              local_swarm[i].position.begin());
    off += dim;
    std::copy(recv_buffer.begin() + base + off,
              recv_buffer.begin() + base + off + dim,
              local_swarm[i].velocity.begin());
    off += dim;
    std::copy(recv_buffer.begin() + base + off,
              recv_buffer.begin() + base + off + dim,
              local_swarm[i].best_position.begin());
    off += dim;
    local_swarm[i].best_value = recv_buffer[base + off++];
    local_swarm[i].current_value = recv_buffer[base + off++];
  }
}

OutputObject dpso(const TestFunction &f, unsigned int dim,
                  unsigned int n_points_total, int max_iter,
                  const DPSOParameters &params) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  StoppingCriteriaManager stop_manager(max_iter, 2000, 1e-8, 1e-3);

  unsigned int n_points_per_rank = n_points_total / size;
  if (rank == 0 && n_points_total % size != 0)
    std::cerr << "Warning: total particles (" << n_points_total
              << ") not divisible by number of ranks (" << size << ").\n";

  if (n_points_per_rank < (unsigned int)params.sub_swarm_size) {
    if (rank == 0) {
      std::cerr << "Error: particles per rank (" << n_points_per_rank
                << ") less than sub-swarm size (" << params.sub_swarm_size
                << ").\n";
    }
    return OutputObject(f.get_name(), dim, n_points_per_rank * size, {},
                        f.get_true_solution(), 0.0, {}, size, 0.0, 0,
                        stop_manager);
  }
  if (n_points_per_rank % params.sub_swarm_size != 0 && rank == 0)
    std::cerr << "Warning: particles per rank (" << n_points_per_rank
              << ") not divisible by sub-swarm size (" << params.sub_swarm_size
              << ").\n";

  const auto &domain = f.get_domain();
  std::vector<double> lb(dim, domain.first);
  std::vector<double> ub(dim, domain.second);
  std::vector<double> v_max(dim);
  if (lb[0] > ub[0])
    for (unsigned int d = 0; d < dim; ++d)
      lb[d] = -ub[d];
  for (unsigned int d = 0; d < dim; ++d)
    v_max[d] = 0.2 * (ub[d] - lb[d]);

  std::vector<Particle> swarm;
  swarm.reserve(n_points_per_rank);
  for (unsigned int i = 0; i < n_points_per_rank; ++i) {
    Particle p(dim);
    for (unsigned int d = 0; d < dim; ++d) {
      p.position[d] = random_double(lb[d], ub[d], rank);
      p.velocity[d] = random_double(-v_max[d], v_max[d], rank);
      p.best_position[d] = p.position[d];
    }
    p.current_value = f.value(p.position);
    p.best_value = p.current_value;
    swarm.push_back(p);
  }

  OutputObject results(f.get_name(), dim, n_points_per_rank * size, {},
                       f.get_true_solution(), 0.0, {}, size, 0.0, 0,
                       stop_manager);
  results.x_best.resize(dim);
  double global_best_val = std::numeric_limits<double>::max();
  double start_time = MPI_Wtime();
  int iter = 0;

  while (true) {
    if (iter > 0 && iter % params.regrouping_period == 0)
      regroup_particles(swarm, dim, rank, size);

    int num_sub_swarms = swarm.size() / params.sub_swarm_size;
    int remainder = swarm.size() % params.sub_swarm_size;

    for (int s = 0; s < num_sub_swarms; ++s) {
      int start = s * params.sub_swarm_size;
      int end = std::min(start + params.sub_swarm_size, (int)swarm.size());
      process_sub_swarm(swarm, start, end, dim, lb, ub, v_max, params, f, rank,
                        iter, max_iter);
    }
    if (remainder > 0) {
      int start = num_sub_swarms * params.sub_swarm_size;
      process_sub_swarm(swarm, start, (int)swarm.size(), dim, lb, ub, v_max,
                        params, f, rank, iter, max_iter);
    }

    // Global best via MPI_MINLOC: obtains value and rank of the minimum
    struct {
      double val;
      int rank;
    } local_min, global_min;
    local_min.val = std::numeric_limits<double>::max();
    local_min.rank = rank;
    int local_best_idx = -1;
    for (int i = 0; i < (int)swarm.size(); ++i) {
      if (swarm[i].best_value < local_min.val) {
        local_min.val = swarm[i].best_value;
        local_best_idx = i;
      }
    }
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  MPI_COMM_WORLD);

    std::vector<double> global_best_pos(dim);
    if (rank == global_min.rank && local_best_idx != -1)
      global_best_pos = swarm[local_best_idx].best_position;
    MPI_Bcast(global_best_pos.data(), dim, MPI_DOUBLE, global_min.rank,
              MPI_COMM_WORLD);

    // Diversity: average distance from the global best, averaged across ranks
    double local_sum = 0.0;
    for (const auto &p : swarm)
      local_sum += euclidean_dist(p.position, global_best_pos);
    double local_avg = swarm.empty() ? 0.0 : local_sum / swarm.size();
    double global_avg_dist = 0.0;
    MPI_Allreduce(&local_avg, &global_avg_dist, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    global_avg_dist /= size;

    if (rank == 0)
      results.conv_history.push_back(global_min.val);
    global_best_val = global_min.val;

    stop_manager.increment_iterations();
    iter++;

    // The stop flag is decided only by rank 0 and then broadcast,
    // so that all ranks exit at the same iteration (avoids MPI deadlocks).
    int stop_flag = 0;
    if (rank == 0)
      stop_flag =
          stop_manager.should_stop(global_best_val, global_avg_dist) ? 1 : 0;
    MPI_Bcast(&stop_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (stop_flag)
      break;
  }

  struct {
    double val;
    int rank;
  } loc, glob;
  loc.val = std::numeric_limits<double>::max();
  loc.rank = rank;
  int best_local_idx = -1;
  for (int i = 0; i < (int)swarm.size(); ++i) {
    if (swarm[i].best_value < loc.val) {
      loc.val = swarm[i].best_value;
      best_local_idx = i;
    }
  }
  MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  if (rank == glob.rank && best_local_idx != -1)
    results.x_best = swarm[best_local_idx].best_position;
  MPI_Bcast(results.x_best.data(), dim, MPI_DOUBLE, glob.rank, MPI_COMM_WORLD);
  results.f_val = glob.val;
  results.execution_time = MPI_Wtime() - start_time;
  results.iterations = iter;
  return results;
}
