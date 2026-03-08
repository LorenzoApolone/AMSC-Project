
#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "pso_topology.hpp"
#include <algorithm>
#include <limits>
#include <mpi.h>
#include <random>
#include <vector>
#include <iostream>

/**
 * @brief Parallel PSO with neighborhood topology and MPI communication.
 *
 */

// Struct containing a serie of parameter needed for the PSO
struct PSOHyperparameters {
  static constexpr double C1 = 1.49618;
  static constexpr double C2 = 1.49618;

  static constexpr double W_MAX = 0.9;
  static constexpr double W_MIN = 0.4;

  static constexpr double V_INIT_FACTOR = 0.1;
  
};

OutputObject pso_topology(const TestFunction &f,
                       int d,
                      StoppingCriteriaManager &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list,  double &t_allgatherv_tot) {

  // MPI Setup
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //  Chronometer start
  double start_time = MPI_Wtime();

  // Particle Distribution
  int local_n = n_points / size;
  int remainder = n_points % size;
  if (rank < remainder)
    local_n++;

  // counts/displs (global indexing)
  std::vector<int> counts(size), displs(size);
  // gather local_n from all ranks to compute displs
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  //  it's used for the gid calculation               
  displs[0] = 0;
  for (int r = 1; r < size; ++r)
    displs[r] = displs[r - 1] + counts[r - 1];

  
  std::vector<double> pos(local_n * d);
  std::vector<double> vel(local_n * d);
  std::vector<double> pbest_pos = pos;
  std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

  // Global best used for the stopping criterion
  std::vector<double> gbest_pos(d);
  double gbest_val = std::numeric_limits<double>::max();

  // History (Rank 0 only)
  std::vector<double> history;

  // Domain bounds
  auto bounds = f.get_domain();
  double LB = bounds.first;
  double UB = bounds.second;

  // Random generators to dinstribute particel on the domain
  std::mt19937 gen(rank + 42 + std::hash<std::string>{}(f.get_name()));
  std::uniform_real_distribution<> dis(LB, UB);
  std::uniform_real_distribution<> dis_01(0.0, 1.0);
  std::uniform_real_distribution<> vel_dis(-1.0, 1.0);
  double range = UB - LB;

  // Distribute particles randomly in the domain, initialize velocities and pbest
  for (int i = 0; i < local_n; ++i) {
    std::vector<double> current_pos(d);
    for (int j = 0; j < d; ++j) {
      int idx = i * d + j;
      pos[idx] = dis(gen);
      vel[idx] = vel_dis(gen) * range * PSOHyperparameters::V_INIT_FACTOR;
      pbest_pos[idx] = pos[idx];
      current_pos[j] = pos[idx];
    }
    double fitness = f.value(current_pos);
    pbest_val[i] = fitness;
    if (fitness < gbest_val) {
      gbest_val = fitness;
      gbest_pos = current_pos;
    }
  }

  // Struct to return results
  struct {
    double val;
    int rank;
  } loc_data, glob_data;

  // Find and share global best among all ranks
  loc_data.val = gbest_val;
  loc_data.rank = rank;
  MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  gbest_val = glob_data.val;
  MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

  // Buffers for Allgatherv of pbests (needed for lbest) ---
  std::vector<double> all_pbest_val(n_points);
  std::vector<double> all_pbest_pos(n_points * d);
  std::vector<int> counts_d(size), displs_d(size);

  for (int r = 0; r < size; ++r) {
    counts_d[r] = counts[r] * d;
    displs_d[r] = displs[r] * d;
  }
  // end of initialization

  // -----------------------------
  // Main Loop 
  // -----------------------------

  int iter = 0;
  bool must_stop = false;
  int max_iter_limit = stop.get_max_iters();

  // Main loop that stop when at least one of the stopping criterion is met 
  while (!must_stop) {

    // incrementing iteration counter of the stopping criteria
    if (rank == 0){
          stop.increment_iterations();

    }

    // Compute inertia weight for the current iteration
    double current_w =
        PSOHyperparameters::W_MAX -
        ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) * (double)iter /
         (double)max_iter_limit);

    // 1) Update particles locally (pos/vel) using the previous neighborhood knowledge

    double t_start = MPI_Wtime();
    MPI_Allgatherv(pbest_val.data(), local_n, MPI_DOUBLE,
                   all_pbest_val.data(), counts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
    
    MPI_Allgatherv(pbest_pos.data(), local_n * d, MPI_DOUBLE,
                   all_pbest_pos.data(), counts_d.data(), displs_d.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
    double t_end = MPI_Wtime();
    // Accumulate Allgatherv time, useful to monitor this bottleneck 
    t_allgatherv_tot += (t_end - t_start);


    // 2) For each local particle, compute lbest from adjacency_list using all_pbest_
    for (int i = 0; i < local_n; ++i) {
      int gid = displs[rank] + i; // global id of this particle

    if (gid < 0 || gid >= n_points) {
        std::cerr << "Invalid gid: " << gid << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
      // Find best among {gid} U neighbors(gid)
      int best_gid = gid;
      double best_val = all_pbest_val[gid];

      // include neighbors
      for (int neigh : adjacency_list[gid]) {
        if (neigh < 0 || neigh >= n_points) {
            std::cerr << "Invalid neighbor id: " << neigh
                      << " for particle " << gid << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (all_pbest_val[neigh] < best_val) {
          best_val = all_pbest_val[neigh];
          best_gid = neigh;
        }
      }

      // Update vel/pos using lbest 
      for (int j = 0; j < d; ++j) {
        double r1 = dis_01(gen);
        double r2 = dis_01(gen);
        double lbest_j = all_pbest_pos[best_gid * d + j];
        int idx = i * d + j;

        vel[idx] =
            current_w * vel[idx] +
            PSOHyperparameters::C1 * r1 * (pbest_pos[idx] - pos[idx]) +
            PSOHyperparameters::C2 * r2 * (lbest_j - pos[idx]);

        pos[idx] += vel[idx];

        // clamp
        if (pos[idx] < LB) pos[idx] = LB;
        if (pos[idx] > UB) pos[idx] = UB;
      }

      // Evaluate and update pbest locally
      std::vector<double> current_pos(d);
      for (int j = 0; j < d; ++j) {
        current_pos[j] = pos[i * d + j];
      }
      double current_fit = f.value(current_pos);

      if (current_fit < pbest_val[i]) {
        pbest_val[i] = current_fit;
        for (int j = 0; j < d; ++j) {
          pbest_pos[i * d + j] = pos[i * d + j];
        }
      }
    }

    // 3) Compute global best for stopping/output 
    
    loc_data.val = std::numeric_limits<double>::max();
    loc_data.rank = rank;

    // local best among this rank's pbests
    double local_best_val = std::numeric_limits<double>::max();
    int local_best_idx = -1;
    for (int i = 0; i < local_n; ++i) {
      if (pbest_val[i] < local_best_val) {
        local_best_val = pbest_val[i];
        local_best_idx = i;
      }
    }

    // reduce which rank has global best pbest value
    loc_data.val = local_best_val;
    loc_data.rank = rank;
    MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    // broadcast gbest position from winner rank
    gbest_val = glob_data.val;

    if (rank == glob_data.rank && local_best_idx >= 0) {
      for (int j = 0; j < d; ++j) {
        gbest_pos[j] = pbest_pos[local_best_idx * d + j];
      }
    }

    MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);
    double local_sum_dist = 0.0;

    for (int i = 0; i < local_n; ++i) {
      double dist_sq = 0.0;
      for (int j = 0; j < d; ++j) {
        double diff = pos[i * d + j] - gbest_pos[j];   
        dist_sq += diff * diff;
      }
      local_sum_dist += std::sqrt(dist_sq);
    }

    double global_sum_dist = 0.0;
    MPI_Allreduce(&local_sum_dist, &global_sum_dist, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double avg_distance = global_sum_dist / static_cast<double>(n_points);

    // 4) Stopping criterion, the check is done only on rank 0.
    int stop_signal = 0;
    if (rank == 0) {
      double err = f.error(gbest_pos);
      history.push_back(err);
      if (stop.should_stop(gbest_val, avg_distance)) {
         stop_signal = 1;
         
        
      }
    }
    // Broadcast stop signal to all ranks
    MPI_Bcast(&stop_signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    must_stop = (stop_signal != 0);

    iter++;
  }

  //------------------------------
  // End of main loop
  //------------------------------

  double end_time = MPI_Wtime();

  OutputObject out(f.get_name(),
                   d,
                   n_points,
                   gbest_pos,
                   f.get_true_solution(),
                   gbest_val,
                   history,
                   size,
                   end_time - start_time,
                   iter,
                   stop);

  return out;
}