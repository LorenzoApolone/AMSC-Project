// file con all gather e senza timer 


#include "../interfaces.hpp"
#include "create_network.hpp"
#include "pso_topology.hpp"
#include <algorithm>
#include <limits>
#include <mpi.h>
#include <random>
#include <vector>
#include <iostream>

struct PSOHyperparameters {
  static constexpr double C1 = 1.49618;
  static constexpr double C2 = 1.49618;

  static constexpr double W_MAX = 0.9;
  static constexpr double W_MIN = 0.4;

  static constexpr double V_INIT_FACTOR = 0.1;
  
};



OutputObject pso_normal(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list, bool &converged) {

  // --- MPI Setup ---
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

  // --- Data Structures ---
  std::vector<std::vector<double>> pos(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> vel(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> pbest_pos = pos;
  std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

  // Global best (used for stopping + output)
  std::vector<double> gbest_pos(d);
  double gbest_val = std::numeric_limits<double>::max();

  // History (Rank 0 only)
  std::vector<double> history;

  // Domain bounds
  auto bounds = f.get_domain();
  double LB = bounds.first;
  double UB = bounds.second;

  // Random generators to dinstribute particel on the domain
  std::mt19937 gen(rank + 42);
  std::uniform_real_distribution<> dis(LB, UB);
  std::uniform_real_distribution<> dis_01(0.0, 1.0);

  // --- Initialization ---
  for (int i = 0; i < local_n; ++i) {
    for (int j = 0; j < d; ++j) {
      pos[i][j] = dis(gen);
      vel[i][j] = (dis(gen) - dis(gen)) * PSOHyperparameters::V_INIT_FACTOR;
      pbest_pos[i][j] = pos[i][j];
    }

    double fitness = f.value(pos[i]);
    pbest_val[i] = fitness;

    if (fitness < gbest_val) {
      gbest_val = fitness;
      gbest_pos = pos[i];
    }
  }

  // --- Initial Sync of gbest (optional, just to have something consistent) ---
  struct {
    double val;
    int rank;
  } loc_data, glob_data;

  loc_data.val = gbest_val;
  loc_data.rank = rank;

  MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  gbest_val = glob_data.val;
  MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

  // --- Buffers for Allgatherv of pbests (needed for lbest) ---
  std::vector<double> all_pbest_val(n_points);
  std::vector<double> all_pbest_pos(n_points * d);

  std::vector<double> local_pbest_pos_flat(local_n * d);

  std::vector<int> counts_d(size), displs_d(size);
  for (int r = 0; r < size; ++r) {
    counts_d[r] = counts[r] * d;
    displs_d[r] = displs[r] * d;
  }
  // end of initialization


  // --- Main Loop ---
  int iter = 0;
  bool must_stop = false;
  int max_iter_limit = stop.get_max_iter();
  int counter =0;
  while (!must_stop) {
    counter++;
    double current_w =
        PSOHyperparameters::W_MAX -
        ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) * (double)iter /
         (double)max_iter_limit);

    // 1) Update particles locally (pos/vel) using the previous neighborhood knowledge
  
    // ---- current pbest (val + pos) to all ranks ----
    for (int i = 0; i < local_n; ++i)
      for (int j = 0; j < d; ++j)
        local_pbest_pos_flat[i * d + j] = pbest_pos[i][j];
    
    MPI_Allgatherv(pbest_val.data(), local_n, MPI_DOUBLE,
                   all_pbest_val.data(), counts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
   

    MPI_Allgatherv(local_pbest_pos_flat.data(), local_n * d, MPI_DOUBLE,
                   all_pbest_pos.data(), counts_d.data(), displs_d.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
   
    // 2) For each local particle, compute lbest from adjacency_list using all_pbest_
    for (int i = 0; i < local_n; ++i) {
      int gid = displs[rank] + i; // global id of this particle

      // Find best among {gid} U neighbors(gid)
      int best_gid = gid;
      double best_val = all_pbest_val[gid];

      // include neighbors
      for (int neigh : adjacency_list[gid]) {
        if (all_pbest_val[neigh] < best_val) {
          best_val = all_pbest_val[neigh];
          best_gid = neigh;
        }
      }

      // Update vel/pos using lbest (best_gid)
      for (int j = 0; j < d; ++j) {
        double r1 = dis_01(gen);
        double r2 = dis_01(gen);

        double lbest_j = all_pbest_pos[best_gid * d + j];

        vel[i][j] =
            current_w * vel[i][j] +
            PSOHyperparameters::C1 * r1 * (pbest_pos[i][j] - pos[i][j]) +
            PSOHyperparameters::C2 * r2 * (lbest_j - pos[i][j]);

        pos[i][j] += vel[i][j];

        // clamp
        if (pos[i][j] < LB) pos[i][j] = LB;
        if (pos[i][j] > UB) pos[i][j] = UB;
      }

      // Evaluate and update pbest locally
      double current_fit = f.value(pos[i]);

      if (current_fit < pbest_val[i]) {
        pbest_val[i] = current_fit;
        pbest_pos[i] = pos[i];
      }
    }

    // 3) Compute global best for stopping/output (rank 0 uses gathered arrays)
    
    loc_data.val = std::numeric_limits<double>::max();
    loc_data.rank = rank;

    // local best among this rank's pbests
    double local_best_val = std::numeric_limits<double>::max();
    int local_best_idx = 0;
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
    if (rank == glob_data.rank) {
      gbest_val = local_best_val;
      gbest_pos = pbest_pos[local_best_idx];
    }
    MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

    // 4) Stopping criterion (rank 0 decides)
    int stop_signal = 0;
    if (rank == 0) {
      double err = f.error(gbest_pos);
      history.push_back(err);
      if (stop.should_stop(iter, err)){
         stop_signal = 1;
         if (iter < stop.get_max_iter())
           converged = true;
        
      }
    }
    MPI_Bcast(&stop_signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    must_stop = (stop_signal != 0);

    iter++;
  }

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