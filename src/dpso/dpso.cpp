/**
 * @file dpso.cpp
 * @brief Distributed DMS-PSO-HS via MPI (Zhao et al., 2011).
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
 * @brief Generates a random double within [min, max).
 *
 * @param min Minimum bound.
 * @param max Maximum bound.
 * @param gen PRNG instance.
 * @return Random double.
 */
double random_double(double min, double max, std::mt19937 &gen) {
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
int random_int(int min, int max, std::mt19937 &gen) {
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
double euclidean_dist_squared(const std::vector<double> &v1,
                              const std::vector<double> &v2) {
  double sum = 0.0;
  for (size_t i = 0; i < v1.size(); ++i) {
    double diff = v1[i] - v2[i];
    sum += diff * diff;
  }
  return sum;
}

/**
 * @brief Applies Harmony Search to a sub-swarm in SoA.
 *
 * @param pbest_pos Swarm personal best positions (SoA).
 * @param pbest_val Swarm personal best values (SoA).
 * @param start_idx Start index of the sub-swarm.
 * @param end_idx End index of the sub-swarm (exclusive).
 * @param dim Dimensions.
 * @param f The test function to optimize.
 * @param gen PRNG instance.
 * @param lower_bound Vector of lower bounds.
 * @param upper_bound Vector of upper bounds.
 * @param current_iter Current iteration number.
 * @param max_iter Maximum number of iterations.
 * @param params Algorithm parameters.
 * @param new_harmony Buffer for harmony search results.
 */
static void apply_harmony_search(const std::vector<double> &pbest_pos,
                                 std::vector<double> &pbest_val, int start_idx,
                                 int end_idx, int dim, const TestFunction &f,
                                 std::mt19937 &gen,
                                 const std::vector<double> &lower_bound,
                                 const std::vector<double> &upper_bound,
                                 int current_iter, int max_iter,
                                 const DPSOParameters &params,
                                 std::vector<double> &new_harmony) {
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

    if (random_double(0.0, 1.0, gen) < params.hmcr) {
      int p_idx = start_idx + random_int(0, sub_swarm_size - 1, gen);
      new_harmony[d] = pbest_pos[p_idx * dim + d];
      if (random_double(0.0, 1.0, gen) < PAR)
        new_harmony[d] += random_double(-1.0, 1.0, gen) * bw;
    } else {
      new_harmony[d] = random_double(lower_bound[d], upper_bound[d], gen);
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
      const_cast<std::vector<double> &>(pbest_pos)[nearest_idx * dim + d] =
          new_harmony[d];
    }
    pbest_val[nearest_idx] = new_val;
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
static void process_sub_swarm(std::vector<double> &pos, std::vector<double> &vel,
                              std::vector<double> &pbest_pos,
                              std::vector<double> &pbest_val,
                              std::vector<double> &current_val, int start,
                              int end, unsigned int dim,
                              const std::vector<double> &lb,
                              const std::vector<double> &ub,
                              const std::vector<double> &v_max,
                              const DPSOParameters &params,
                              const TestFunction &f, std::mt19937 &gen,
                              int iter, int max_iter,
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
      double r1 = random_double(0.0, 1.0, gen);
      double r2 = random_double(0.0, 1.0, gen);
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
  apply_harmony_search(pbest_pos, pbest_val, start, end, dim, f, gen, lb, ub,
                       iter, max_iter, params, hs_buffer);
}

/**
 * @brief Context structure to preserve regrouping buffers without static
 * locals.
 */
struct RegroupContext {
  std::vector<double> send_buffer;
  std::vector<int> recvcounts;
  std::vector<int> displs;
  std::vector<double> recv_buffer;
  std::vector<int> indices;
};

/**
 * @brief Shuffles the particles in the swarm to regroup sub-swarms in SoA.
 *
 * @param pos local swarm positions.
 * @param vel local swarm velocities.
 * @param pbest_pos local swarm personal best positions.
 * @param pbest_val local swarm personal best values.
 * @param current_val local swarm current values.
 * @param dim Number of dimensions.
 * @param rank MPI rank.
 * @param size Total MPI processes.
 * @param g PRNG instance.
 * @param total_particles Total particles in global swarm.
 * @param ctx Context for persistent buffers.
 */
void regroup_particles(std::vector<double> &pos, std::vector<double> &vel,
                       std::vector<double> &pbest_pos,
                       std::vector<double> &pbest_val,
                       std::vector<double> &current_val, int dim, int rank,
                       int size, std::mt19937 &g, int total_particles,
                       RegroupContext &ctx) {
  int local_n = current_val.size();
  int p_data_size = 3 * dim + 2;

  ctx.send_buffer.clear();
  ctx.send_buffer.reserve(local_n * p_data_size);
  for (int i = 0; i < local_n; ++i) {
    ctx.send_buffer.insert(ctx.send_buffer.end(), pos.begin() + i * dim,
                           pos.begin() + (i + 1) * dim);
    ctx.send_buffer.insert(ctx.send_buffer.end(), vel.begin() + i * dim,
                           vel.begin() + (i + 1) * dim);
    ctx.send_buffer.insert(ctx.send_buffer.end(), pbest_pos.begin() + i * dim,
                           pbest_pos.begin() + (i + 1) * dim);
    ctx.send_buffer.push_back(pbest_val[i]);
    ctx.send_buffer.push_back(current_val[i]);
  }

  ctx.recvcounts.resize(size);
  ctx.displs.resize(size);
  int current_displ = 0;
  for (int i = 0; i < size; ++i) {
    int count = total_particles / size;
    if (i < total_particles % size)
      count++;
    ctx.recvcounts[i] = count * p_data_size;
    ctx.displs[i] = current_displ;
    current_displ += ctx.recvcounts[i];
  }

  ctx.recv_buffer.resize(total_particles * p_data_size);
  MPI_Allgatherv(ctx.send_buffer.data(), ctx.send_buffer.size(), MPI_DOUBLE,
                 ctx.recv_buffer.data(), ctx.recvcounts.data(),
                 ctx.displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

  if (rank == 0) {
    ctx.indices.resize(total_particles);
    std::iota(ctx.indices.begin(), ctx.indices.end(), 0);
    std::shuffle(ctx.indices.begin(), ctx.indices.end(), g);
  } else {
    ctx.indices.resize(total_particles);
  }

  MPI_Bcast(ctx.indices.data(), total_particles, MPI_INT, 0, MPI_COMM_WORLD);

  int global_idx_start = 0;
  for (int i = 0; i < rank; ++i) {
    int count = total_particles / size;
    if (i < total_particles % size)
      count++;
    global_idx_start += count;
  }

  for (int i = 0; i < local_n; ++i) {
    int picked = ctx.indices[global_idx_start + i];
    int base = picked * p_data_size;
    int off = 0;
    std::copy_n(ctx.recv_buffer.begin() + base + off, dim, pos.begin() + i * dim);
    off += dim;
    std::copy_n(ctx.recv_buffer.begin() + base + off, dim, vel.begin() + i * dim);
    off += dim;
    std::copy_n(ctx.recv_buffer.begin() + base + off, dim,
                pbest_pos.begin() + i * dim);
    off += dim;
    pbest_val[i] = ctx.recv_buffer[base + off++];
    current_val[i] = ctx.recv_buffer[base + off++];
  }
}

/**
 * @brief Executes the Distributed DMS-PSO-HS algorithm.
 *
 * Coordinates the evaluation, particle updates, and communication
 * over MPI to find the global minimum of the provided test function.
 *
 * @param f The test function to optimize.
 * @param dim Number of dimensions of the search space.
 * @param n_points_total Total number of particles in the global swarm.
 * @param max_iter Base reference for maximum number of iterations.
 * @param params Algorithm parameters (w, c1, c2, hmcr, par_min/max, etc.).
 * @param convergence_tol The tolerance used by the stopping criteria manager.
 * @return OutputObject containing the global best position, value, and
 * performance metrics.
 */
OutputObject dpso(const TestFunction &f, unsigned int dim,
                  unsigned int n_points_total, int max_iter,
                  const DPSOParameters &params, double convergence_tol) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  StoppingCriteriaManager stop_manager(max_iter, 2000, convergence_tol, 1e-3);

  unsigned int n_points_per_rank = n_points_total / size;
  if ((unsigned int)rank < (n_points_total % size)) {
    n_points_per_rank++;
  }

  if (n_points_per_rank < (unsigned int)params.sub_swarm_size) {
    if (rank == 0) {
      std::cerr << "Error: particles per rank (" << n_points_per_rank
                << ") less than sub-swarm size (" << params.sub_swarm_size
                << ").\n";
    }
    return OutputObject(f.get_name(), dim, n_points_total, {},
                        f.get_true_solution(), 0.0, {}, size, 0.0, 0,
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
  std::mt19937 gen(rd() + rank);

  unsigned int global_seed = 0;
  if (rank == 0) {
    global_seed = rd();
  }
  MPI_Bcast(&global_seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  std::mt19937 global_gen(global_seed);

  std::vector<double> pos(n_points_per_rank * dim);
  std::vector<double> vel(n_points_per_rank * dim);
  std::vector<double> pbest_pos(n_points_per_rank * dim);
  std::vector<double> pbest_val(n_points_per_rank,
                                std::numeric_limits<double>::max());
  std::vector<double> current_val(n_points_per_rank);

  for (unsigned int i = 0; i < n_points_per_rank; ++i) {
    std::vector<double> p_pos(dim);
    for (unsigned int d = 0; d < dim; ++d) {
      pos[i * dim + d] = random_double(lb[d], ub[d], gen);
      vel[i * dim + d] = random_double(-v_max[d], v_max[d], gen);
      pbest_pos[i * dim + d] = pos[i * dim + d];
      p_pos[d] = pos[i * dim + d];
    }
    current_val[i] = f.value(p_pos);
    pbest_val[i] = current_val[i];
  }

  OutputObject results(f.get_name(), dim, n_points_total, {},
                       f.get_true_solution(), 0.0, {}, size, 0.0, 0,
                       stop_manager);
  results.x_best.resize(dim);
  double global_best_val = std::numeric_limits<double>::max();
  double start_time = MPI_Wtime();
  int iter = 0;
  RegroupContext regroup_ctx;
  std::vector<double> hs_buffer(dim);

  while (true) {
    if (iter > 0 && iter % params.regrouping_period == 0)
      regroup_particles(pos, vel, pbest_pos, pbest_val, current_val, dim, rank,
                        size, global_gen, n_points_total, regroup_ctx);

    int num_sub_swarms = n_points_per_rank / params.sub_swarm_size;
    int remainder = n_points_per_rank % params.sub_swarm_size;

    for (int s = 0; s < num_sub_swarms; ++s) {
      int start = s * params.sub_swarm_size;
      int end = start + params.sub_swarm_size;
      process_sub_swarm(pos, vel, pbest_pos, pbest_val, current_val, start, end,
                        dim, lb, ub, v_max, params, f, gen, iter, max_iter,
                        hs_buffer);
    }
    if (remainder > 0) {
      int start = num_sub_swarms * params.sub_swarm_size;
      process_sub_swarm(pos, vel, pbest_pos, pbest_val, current_val, start,
                        n_points_per_rank, dim, lb, ub, v_max, params, f, gen,
                        iter, max_iter, hs_buffer);
    }

    struct {
      double val;
      int rank;
    } local_min, global_min;
    local_min.val = std::numeric_limits<double>::max();
    local_min.rank = rank;
    int local_best_idx = -1;
    for (int i = 0; i < (int)n_points_per_rank; ++i) {
      if (pbest_val[i] < local_min.val) {
        local_min.val = pbest_val[i];
        local_best_idx = i;
      }
    }
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                  MPI_COMM_WORLD);

    std::vector<double> global_best_pos(dim);
    if (rank == global_min.rank && local_best_idx != -1) {
      for (unsigned int d = 0; d < dim; ++d)
        global_best_pos[d] = pbest_pos[local_best_idx * dim + d];
    }
    MPI_Bcast(global_best_pos.data(), dim, MPI_DOUBLE, global_min.rank,
              MPI_COMM_WORLD);

    double local_sum = 0.0;
    for (unsigned int i = 0; i < n_points_per_rank; ++i) {
      double d_sq = 0.0;
      for (unsigned int d = 0; d < dim; ++d) {
        double diff = pos[i * dim + d] - global_best_pos[d];
        d_sq += diff * diff;
      }
      local_sum += std::sqrt(d_sq);
    }

    double global_sum_dist = 0.0;
    MPI_Allreduce(&local_sum, &global_sum_dist, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    double global_avg_dist = global_sum_dist / n_points_total;

    if (rank == 0)
      results.conv_history.push_back(global_min.val);
    global_best_val = global_min.val;

    stop_manager.increment_iterations();
    iter++;

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
  for (int i = 0; i < (int)n_points_per_rank; ++i) {
    if (pbest_val[i] < loc.val) {
      loc.val = pbest_val[i];
      best_local_idx = i;
    }
  }
  MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  if (rank == glob.rank && best_local_idx != -1) {
    for (unsigned int d = 0; d < dim; ++d)
      results.x_best[d] = pbest_pos[best_local_idx * dim + d];
  }
  MPI_Bcast(results.x_best.data(), dim, MPI_DOUBLE, glob.rank, MPI_COMM_WORLD);
  results.f_val = glob.val;
  results.execution_time = MPI_Wtime() - start_time;
  results.iterations = iter;
  return results;
}
