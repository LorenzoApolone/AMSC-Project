/**
 * @file dpso.cpp
 * @brief Distributed DMS-PSO-HS via MPI (Zhao et al., 2011).
 */

#include "../interfaces/interfaces.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include "methods_dpso.hpp"
#include "../pso/particle.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
template <typename Fn> void time_mpi_call(double &accumulator, Fn &&fn) {
  const double start = MPI_Wtime();
  fn();
  accumulator += MPI_Wtime() - start;
}
} // namespace

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
static void apply_harmony_search(
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
static void process_sub_swarm(
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
  apply_harmony_search(pos, current_val, pbest_pos, pbest_val, start, end, dim,
                       f, gen, lb, ub, iter, max_iter, params, hs_buffer);
}

/**
 * @brief Context structure to preserve regrouping buffers without static
 * locals.
 */
struct RegroupContext {
  std::vector<double> send_buffer;
  std::vector<int> sendcounts;
  std::vector<int> sdispls;
  std::vector<int> recvcounts;
  std::vector<int> rdispls;
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
                       RegroupContext &ctx, double &t_bcast,
                       double &t_alltoall) {
  int local_n = current_val.size();
  int p_data_size = 3 * dim + 2;

  std::vector<int> counts(size);
  std::vector<int> starts(size);
  int current_start = 0;
  for (int i = 0; i < size; ++i) {
    int count = total_particles / size;
    if (i < total_particles % size)
      count++;
    counts[i] = count;
    starts[i] = current_start;
    current_start += count;
  }

  ctx.indices.resize(total_particles);
  std::iota(ctx.indices.begin(), ctx.indices.end(), 0);
  std::shuffle(ctx.indices.begin(), ctx.indices.end(), g);

  std::vector<int> inv_indices(total_particles);
  for (int k = 0; k < total_particles; ++k) {
    inv_indices[ctx.indices[k]] = k;
  }

  ctx.sendcounts.assign(size, 0);
  ctx.recvcounts.assign(size, 0);

  std::vector<int> dest_ranks(local_n);
  for (int j = 0; j < local_n; ++j) {
    int orig_g = starts[rank] + j;
    int k = inv_indices[orig_g];
    int d_rank = 0;
    while (d_rank < size - 1 && k >= starts[d_rank + 1]) {
      d_rank++;
    }
    dest_ranks[j] = d_rank;
    ctx.sendcounts[d_rank] += p_data_size;
  }

  std::vector<int> src_ranks(local_n);
  for (int j = 0; j < local_n; ++j) {
    int k = starts[rank] + j;
    int orig_g = ctx.indices[k];
    int s_rank = 0;
    while (s_rank < size - 1 && orig_g >= starts[s_rank + 1]) {
      s_rank++;
    }
    src_ranks[j] = s_rank;
    ctx.recvcounts[s_rank] += p_data_size;
  }

  ctx.sdispls.assign(size, 0);
  ctx.rdispls.assign(size, 0);
  for (int i = 1; i < size; ++i) {
    ctx.sdispls[i] = ctx.sdispls[i - 1] + ctx.sendcounts[i - 1];
    ctx.rdispls[i] = ctx.rdispls[i - 1] + ctx.recvcounts[i - 1];
  }

  ctx.send_buffer.resize(local_n * p_data_size);
  std::vector<int> current_sdispl = ctx.sdispls;
  for (int j = 0; j < local_n; ++j) {
    int d_rank = dest_ranks[j];
    int offset = current_sdispl[d_rank];

    std::copy_n(pos.begin() + j * dim, dim, ctx.send_buffer.begin() + offset);
    offset += dim;
    std::copy_n(vel.begin() + j * dim, dim, ctx.send_buffer.begin() + offset);
    offset += dim;
    std::copy_n(pbest_pos.begin() + j * dim, dim,
                ctx.send_buffer.begin() + offset);
    offset += dim;
    ctx.send_buffer[offset++] = pbest_val[j];
    ctx.send_buffer[offset++] = current_val[j];

    current_sdispl[d_rank] += p_data_size;
  }

  ctx.recv_buffer.resize(local_n * p_data_size);

  time_mpi_call(t_alltoall, [&]() {
    MPI_Alltoallv(ctx.send_buffer.data(), ctx.sendcounts.data(),
                  ctx.sdispls.data(), MPI_DOUBLE, ctx.recv_buffer.data(),
                  ctx.recvcounts.data(), ctx.rdispls.data(), MPI_DOUBLE,
                  MPI_COMM_WORLD);
  });

  std::vector<std::vector<int>> recv_j_by_src(size);
  for (int j = 0; j < local_n; ++j) {
    recv_j_by_src[src_ranks[j]].push_back(j);
  }

  for (int s = 0; s < size; ++s) {
    if (recv_j_by_src[s].empty())
      continue;

    std::sort(recv_j_by_src[s].begin(), recv_j_by_src[s].end(),
              [&](int j1, int j2) {
                return ctx.indices[starts[rank] + j1] <
                       ctx.indices[starts[rank] + j2];
              });

    int offset = ctx.rdispls[s];
    for (int j : recv_j_by_src[s]) {
      std::copy_n(ctx.recv_buffer.begin() + offset, dim, pos.begin() + j * dim);
      offset += dim;
      std::copy_n(ctx.recv_buffer.begin() + offset, dim, vel.begin() + j * dim);
      offset += dim;
      std::copy_n(ctx.recv_buffer.begin() + offset, dim,
                  pbest_pos.begin() + j * dim);
      offset += dim;
      pbest_val[j] = ctx.recv_buffer[offset++];
      current_val[j] = ctx.recv_buffer[offset++];
    }
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
                  const DPSOParameters &params, double convergence_tol,
                  unsigned int seed) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double t_bcast = 0.0;
  double t_allreduce = 0.0;
  double t_alltoall = 0.0;

  int scaled_stagnation = std::max(100, max_iter / 4);
  StoppingCriteriaManager stop_manager(max_iter, scaled_stagnation, 1e-8, 1e-6);

  unsigned int n_points_per_rank = n_points_total / size;
  if ((unsigned int)rank < (n_points_total % size)) {
    n_points_per_rank++;
  }

  if (n_points_per_rank < (unsigned int)params.sub_swarm_size) {
    throw std::invalid_argument("Error: particles per rank (" +
                                std::to_string(n_points_per_rank) +
                                ") less than sub-swarm size (" +
                                std::to_string(params.sub_swarm_size) + ").");
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
  std::mt19937 gen(seed == 0 ? rd() + rank : seed + rank);

  unsigned int global_seed = 0;
  if (rank == 0) {
    global_seed = seed == 0 ? rd() : seed;
  }
  time_mpi_call(t_bcast, [&]() {
    MPI_Bcast(&global_seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  });
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

  const int SYNC_PERIOD = std::max(1, std::min(50, max_iter / 200));
  std::vector<double> global_best_pos(dim, 0.0);
  double global_avg_dist = 0.0;

  while (true) {
    if (iter > 0 && iter % params.regrouping_period == 0)
      regroup_particles(pos, vel, pbest_pos, pbest_val, current_val, dim, rank,
                        size, global_gen, n_points_total, regroup_ctx, t_bcast,
                        t_alltoall);

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

    bool do_sync = (iter % SYNC_PERIOD == 0) || (iter == 0);

    if (do_sync) {
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
      time_mpi_call(t_allreduce, [&]() {
        MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC,
                      MPI_COMM_WORLD);
      });

      if (rank == global_min.rank && local_best_idx != -1) {
        for (unsigned int d = 0; d < dim; ++d)
          global_best_pos[d] = pbest_pos[local_best_idx * dim + d];
      }
      time_mpi_call(t_bcast, [&]() {
        MPI_Bcast(global_best_pos.data(), dim, MPI_DOUBLE, global_min.rank,
                  MPI_COMM_WORLD);
      });

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
      time_mpi_call(t_allreduce, [&]() {
        MPI_Allreduce(&local_sum, &global_sum_dist, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
      });
      global_avg_dist = global_sum_dist / n_points_total;

      if (rank == 0)
        results.conv_history.push_back(global_min.val);
      global_best_val = global_min.val;
    } else {
      double local_best = std::numeric_limits<double>::max();
      for (int i = 0; i < (int)n_points_per_rank; ++i) {
        if (pbest_val[i] < local_best)
          local_best = pbest_val[i];
      }
      if (local_best < global_best_val)
        global_best_val = local_best;
    }

    stop_manager.increment_iterations();
    iter++;

    if (do_sync && stop_manager.should_stop(global_best_val, global_avg_dist))
      break;
    if (iter >= max_iter)
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
  time_mpi_call(t_allreduce, [&]() {
    MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  });
  if (rank == glob.rank && best_local_idx != -1) {
    for (unsigned int d = 0; d < dim; ++d)
      results.x_best[d] = pbest_pos[best_local_idx * dim + d];
  }
  time_mpi_call(t_bcast, [&]() {
    MPI_Bcast(results.x_best.data(), dim, MPI_DOUBLE, glob.rank,
              MPI_COMM_WORLD);
  });
  results.f_val = glob.val;
  results.execution_time = MPI_Wtime() - start_time;
  results.iterations = iter;

  results.comm_bcast_s = t_bcast;
  results.comm_allreduce_s = t_allreduce;
  results.comm_allgather_s = t_alltoall;
  results.comm_barrier_s = 0.0;
  results.wait_total_s = 0.0;
  results.comm_total_s = t_bcast + t_allreduce + t_alltoall;
  results.compute_total_s = results.execution_time - results.comm_total_s;

  return results;
}
