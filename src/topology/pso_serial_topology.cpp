#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "pso_topology.hpp"
#include <algorithm>
#include <limits>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

/**
 * @brief Serial PSO with neighborhood topology.
 */
struct PSOHyperparameters {
  static constexpr double C1 = 1.49618;
  static constexpr double C2 = 1.49618;

  static constexpr double W_MAX = 0.9;
  static constexpr double W_MIN = 0.4;

  static constexpr double V_INIT_FACTOR = 0.1;
};

OutputObject pso_serial_topology(const TestFunction &f,
                          int d,
                          StoppingCriteriaManager &stop,
                          int n_points,
                          const std::vector<std::vector<int>> &adjacency_list) {


  // Chronometer start
  auto start_time = std::chrono::steady_clock::now();

  // All particles are local in serial
  int local_n = n_points;

  // Data structures
  std::vector<std::vector<double>> pos(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> vel(local_n, std::vector<double>(d));
  std::vector<std::vector<double>> pbest_pos(local_n, std::vector<double>(d));
  std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

  // Global best
  std::vector<double> gbest_pos(d);
  double gbest_val = std::numeric_limits<double>::max();

  // History
  std::vector<double> history;

  // Domain bounds
  auto bounds = f.get_domain();
  double LB = bounds.first;
  double UB = bounds.second;
  double range = UB - LB;

  // Random generators
  std::mt19937 gen( 42 + std::hash<std::string>{}(f.get_name()));
  std::uniform_real_distribution<> dis(LB, UB);
  std::uniform_real_distribution<> dis_01(0.0, 1.0);
  std::uniform_real_distribution<> vel_dis(-1.0, 1.0);

  // Initialization
  for (int i = 0; i < local_n; ++i) {
    for (int j = 0; j < d; ++j) {
      pos[i][j] = dis(gen);
      vel[i][j] = vel_dis(gen) * range * PSOHyperparameters::V_INIT_FACTOR;
      pbest_pos[i][j] = pos[i][j];
    }

    double fitness = f.value(pos[i]);
    pbest_val[i] = fitness;

    if (fitness < gbest_val) {
      gbest_val = fitness;
      gbest_pos = pos[i];
    }
  }

  // Main loop
  int iter = 0;
  bool must_stop = false;
  int max_iter_limit = stop.get_max_iters();

  while (!must_stop) {
    stop.increment_iterations();

    double current_w =
        PSOHyperparameters::W_MAX -
        ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) *
         static_cast<double>(iter) / static_cast<double>(max_iter_limit));

    // For each particle compute lbest from adjacency list
    for (int i = 0; i < local_n; ++i) {
      int gid = i; // in serial local id == global id

      if (gid < 0 || gid >= n_points) {
        std::cerr << "Invalid gid: " << gid << std::endl;
        std::exit(1);
      }

      // Find best among {gid} U neighbors(gid)
      int best_gid = gid;
      double best_val = pbest_val[gid];

      for (int neigh : adjacency_list[gid]) {
        if (neigh < 0 || neigh >= n_points) {
          std::cerr << "Invalid neighbor id: " << neigh
                    << " for particle " << gid << std::endl;
          std::exit(1);
        }

        if (pbest_val[neigh] < best_val) {
          best_val = pbest_val[neigh];
          best_gid = neigh;
        }
      }

      // Update velocity and position using lbest
      for (int j = 0; j < d; ++j) {
        double r1 = dis_01(gen);
        double r2 = dis_01(gen);
        double lbest_j = pbest_pos[best_gid][j];

        vel[i][j] =
            current_w * vel[i][j] +
            PSOHyperparameters::C1 * r1 * (pbest_pos[i][j] - pos[i][j]) +
            PSOHyperparameters::C2 * r2 * (lbest_j - pos[i][j]);

        pos[i][j] += vel[i][j];

        // Clamp
        if (pos[i][j] < LB) pos[i][j] = LB;
        if (pos[i][j] > UB) pos[i][j] = UB;
      }

      // Evaluate and update pbest
      double current_fit = f.value(pos[i]);

      if (current_fit < pbest_val[i]) {
        pbest_val[i] = current_fit;
        pbest_pos[i] = pos[i];
      }
    }

    // Recompute global best
    gbest_val = std::numeric_limits<double>::max();
    for (int i = 0; i < local_n; ++i) {
      if (pbest_val[i] < gbest_val) {
        gbest_val = pbest_val[i];
        gbest_pos = pbest_pos[i];
      }
    }

    // Compute average distance from global best
    double sum_dist = 0.0;
    for (int i = 0; i < local_n; ++i) {
      double dist_sq = 0.0;
      for (int j = 0; j < d; ++j) {
        double diff = pos[i][j] - gbest_pos[j];
        dist_sq += diff * diff;
      }
      sum_dist += std::sqrt(dist_sq);
    }

    double avg_distance = sum_dist / static_cast<double>(n_points);

    // Stopping criterion
    double err = f.error(gbest_pos);
    history.push_back(err);

    if (stop.should_stop(gbest_val, avg_distance)) {
      must_stop = true;
    }

    iter++;
  }

  auto end_time = std::chrono::steady_clock::now();
  double elapsed_time = std::chrono::duration<double>(end_time - start_time).count();
  OutputObject out(f.get_name(),
                   d,
                   n_points,
                   gbest_pos,
                   f.get_true_solution(),
                   gbest_val,
                   history,
                   1, 
                   elapsed_time,
                   iter,
                   stop);

  return out;
}