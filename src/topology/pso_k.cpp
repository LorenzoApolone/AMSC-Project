#include "../interfaces.hpp"
#include "create_network.hpp"
#include <algorithm>
#include <limits>
#include <mpi.h>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct PSOHyperparameters {
  static constexpr double C1 = 1.49618;
  static constexpr double C2 = 1.49618;

  static constexpr double W_MAX = 0.9;
  static constexpr double W_MIN = 0.4;

  static constexpr double V_INIT_FACTOR = 0.1;
};

OutputObject pso_small_debugger_2(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list,
                       bool &converged)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double start_time = MPI_Wtime();

  // ---------------- Particle distribution ----------------
  int local_n = n_points / size;
  int remainder = n_points % size;
  if (rank < remainder) local_n++;

  std::vector<int> counts(size), displs(size);
  MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  displs[0] = 0;
  for (int r = 1; r < size; ++r)
    displs[r] = displs[r - 1] + counts[r - 1];

  auto owner_of = [&](int gid) -> int {
    int r = int(std::upper_bound(displs.begin(), displs.end(), gid) - displs.begin()) - 1;
    if (r < 0) r = 0;
    if (r >= size) r = size - 1;
    return r;
  };

  // ---------------- Data structures ----------------
  std::vector<std::vector<double>> pos(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> vel(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> pbest_pos(local_n, std::vector<double>(d));
  std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

  std::vector<double> gbest_pos(d);
  double gbest_val = std::numeric_limits<double>::max();

  std::vector<double> history; // rank 0 only

  auto bounds = f.get_domain();
  double LB = bounds.first;
  double UB = bounds.second;

  std::mt19937 gen(rank + 42);
  std::uniform_real_distribution<> dis(LB, UB);
  std::uniform_real_distribution<> dis01(0.0, 1.0);

  // ---------------- Initialization ----------------
  for (int i = 0; i < local_n; ++i) {
    for (int j = 0; j < d; ++j) {
      pos[i][j] = dis(gen);
      vel[i][j] = (dis(gen) - dis(gen)) * PSOHyperparameters::V_INIT_FACTOR;
      pbest_pos[i][j] = pos[i][j];
    }
    double fit = f.value(pos[i]);
    pbest_val[i] = fit;
    if (fit < gbest_val) {
      gbest_val = fit;
      gbest_pos = pos[i];
    }
  }

  // Sync initial gbest
  struct { double val; int rank; } loc_data, glob_data;
  loc_data.val = gbest_val;
  loc_data.rank = rank;
  MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  gbest_val = glob_data.val;
  MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

  // =====================================================================
  // PREPROCESSING (one-time): ghost exchange plan
  // =====================================================================

  // need_ids_set[p]: remote GIDs (owned by p) needed by THIS rank
  std::vector<std::unordered_set<int>> need_ids_set(size);

  for (int i = 0; i < local_n; ++i) {
    int gid = displs[rank] + i;
    for (int neigh : adjacency_list[gid]) {
      int ow = owner_of(neigh);
      if (ow != rank) need_ids_set[ow].insert(neigh);
    }
  }

  std::vector<std::vector<int>> need_ids_by_peer(size);
  for (int p = 0; p < size; ++p) {
    need_ids_by_peer[p].assign(need_ids_set[p].begin(), need_ids_set[p].end());
    std::sort(need_ids_by_peer[p].begin(), need_ids_by_peer[p].end());
  }

  // Build request send buffers (IDs we request FROM each peer)
  std::vector<int> req_sendcounts(size, 0), req_sdispls(size, 0);
  for (int p = 0; p < size; ++p) req_sendcounts[p] = (int)need_ids_by_peer[p].size();
  std::partial_sum(req_sendcounts.begin(), req_sendcounts.end() - 1, req_sdispls.begin() + 1);

  int req_send_total = std::accumulate(req_sendcounts.begin(), req_sendcounts.end(), 0);
  std::vector<int> req_sendbuf(req_send_total);
  for (int p = 0; p < size; ++p) {
    std::copy(need_ids_by_peer[p].begin(), need_ids_by_peer[p].end(),
              req_sendbuf.begin() + req_sdispls[p]);
  }

  // Receive counts of requests incoming to us
  std::vector<int> req_recvcounts(size, 0), req_rdispls(size, 0);
  MPI_Alltoall(req_sendcounts.data(), 1, MPI_INT,
               req_recvcounts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  std::partial_sum(req_recvcounts.begin(), req_recvcounts.end() - 1, req_rdispls.begin() + 1);

  int req_recv_total = std::accumulate(req_recvcounts.begin(), req_recvcounts.end(), 0);
  std::vector<int> req_recvbuf(req_recv_total); // IDs that other ranks request FROM ME

  // Exchange request IDs
  MPI_Alltoallv(req_sendbuf.data(), req_sendcounts.data(), req_sdispls.data(), MPI_INT,
                req_recvbuf.data(), req_recvcounts.data(), req_rdispls.data(), MPI_INT,
                MPI_COMM_WORLD);

  // Map: requested remote gid -> index in ghost arrays (order == req_sendbuf)
  std::unordered_map<int, int> ghost_index;
  ghost_index.reserve((size_t)req_send_total * 2 + 1);
  for (int idx = 0; idx < req_send_total; ++idx)
    ghost_index[req_sendbuf[idx]] = idx;

  // Ghost arrays we receive each iteration (for req_sendbuf ids)
  std::vector<double> ghost_pbest_val(req_send_total, std::numeric_limits<double>::max());
  std::vector<double> ghost_pbest_pos((size_t)req_send_total * (size_t)d, 0.0);

  // Precompute double-counts/displs for position exchange
  std::vector<int> sendcounts_d(size, 0), sdispls_d(size, 0), recvcounts_d(size, 0), rdispls_d(size, 0);
  for (int p = 0; p < size; ++p) {
    sendcounts_d[p] = req_recvcounts[p] * d; // we send positions for ids requested from us
    recvcounts_d[p] = req_sendcounts[p] * d; // we receive positions for ids we requested
    sdispls_d[p] = req_rdispls[p] * d;
    rdispls_d[p] = req_sdispls[p] * d;
  }

  // =====================================================================
  // MAIN LOOP
  // =====================================================================
  int iter = 0;
  bool must_stop = false;
  const int max_iter_limit = stop.get_max_iter();

  while (!must_stop) {

    double current_w =
        PSOHyperparameters::W_MAX -
        ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) * (double)iter /
         (double)max_iter_limit);

    // ------------------------------------------------------------
    // (A) Ghost exchange: send pbests that others asked from us
    // ------------------------------------------------------------

    // --- exchange pbest_val ---
    std::vector<double> send_vals(req_recv_total);
    for (int k = 0; k < req_recv_total; ++k) {
      int gid = req_recvbuf[k];        // should be owned by this rank
      int li  = gid - displs[rank];    // local index
      send_vals[k] = pbest_val[li];
    }

    std::vector<double> recv_vals(req_send_total);

    MPI_Alltoallv(send_vals.data(), req_recvcounts.data(), req_rdispls.data(), MPI_DOUBLE,
                  recv_vals.data(), req_sendcounts.data(), req_sdispls.data(), MPI_DOUBLE,
                  MPI_COMM_WORLD);

    ghost_pbest_val = std::move(recv_vals);

    // --- exchange pbest_pos (flattened) ---
    std::vector<double> send_pos((size_t)req_recv_total * (size_t)d);
    for (int k = 0; k < req_recv_total; ++k) {
      int gid = req_recvbuf[k];
      int li  = gid - displs[rank];
      for (int j = 0; j < d; ++j)
        send_pos[(size_t)k * (size_t)d + (size_t)j] = pbest_pos[li][j];
    }

    std::vector<double> recv_pos((size_t)req_send_total * (size_t)d);

    MPI_Alltoallv(send_pos.data(), sendcounts_d.data(), sdispls_d.data(), MPI_DOUBLE,
                  recv_pos.data(), recvcounts_d.data(), rdispls_d.data(), MPI_DOUBLE,
                  MPI_COMM_WORLD);

    ghost_pbest_pos = std::move(recv_pos);

    // ------------------------------------------------------------
    // (B) Update local particles using lbest from local + ghost
    // ------------------------------------------------------------
    for (int i = 0; i < local_n; ++i) {
      int gid = displs[rank] + i;

      int best_gid = gid;
      double best_val = pbest_val[i];

      // scan neighbors
      for (int neigh : adjacency_list[gid]) {
        int ow = owner_of(neigh);

        double neigh_val;
        if (ow == rank) {
          int li = neigh - displs[rank];
          neigh_val = pbest_val[li];
        } else {
          auto it = ghost_index.find(neigh);
          if (it == ghost_index.end()) continue; // safety (should not happen if plan is correct)
          neigh_val = ghost_pbest_val[it->second];
        }

        if (neigh_val < best_val) {
          best_val = neigh_val;
          best_gid = neigh;
        }
      }

      // update using lbest position
      for (int j = 0; j < d; ++j) {
        double lbest_j;
        int ow = owner_of(best_gid);
        if (ow == rank) {
          int li = best_gid - displs[rank];
          lbest_j = pbest_pos[li][j];
        } else {
          auto it = ghost_index.find(best_gid);
          if (it == ghost_index.end()) lbest_j = pbest_pos[i][j]; // fallback
          else lbest_j = ghost_pbest_pos[(size_t)it->second * (size_t)d + (size_t)j];
        }

        double r1 = dis01(gen);
        double r2 = dis01(gen);

        vel[i][j] =
            current_w * vel[i][j] +
            PSOHyperparameters::C1 * r1 * (pbest_pos[i][j] - pos[i][j]) +
            PSOHyperparameters::C2 * r2 * (lbest_j - pos[i][j]);

        pos[i][j] += vel[i][j];

        if (pos[i][j] < LB) pos[i][j] = LB;
        if (pos[i][j] > UB) pos[i][j] = UB;
      }

      double fit = f.value(pos[i]);
      if (fit < pbest_val[i]) {
        pbest_val[i] = fit;
        pbest_pos[i] = pos[i];
      }
    }

    // ------------------------------------------------------------
    // (C) Global best for stopping/output
    // ------------------------------------------------------------
    double local_best_val = std::numeric_limits<double>::max();
    int local_best_idx = 0;
    for (int i = 0; i < local_n; ++i) {
      if (pbest_val[i] < local_best_val) {
        local_best_val = pbest_val[i];
        local_best_idx = i;
      }
    }

    loc_data.val = local_best_val;
    loc_data.rank = rank;
    MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    if (rank == glob_data.rank) {
      gbest_val = local_best_val;
      gbest_pos = pbest_pos[local_best_idx];
    }
    MPI_Bcast(gbest_pos.data(), d, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);

    // ------------------------------------------------------------
    // (D) Stop (rank 0 decides)
    // ------------------------------------------------------------
    int stop_signal = 0;
    if (rank == 0) {
      double err = f.error(gbest_pos);
      history.push_back(err);
      if (stop.should_stop(iter, err)) {
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
