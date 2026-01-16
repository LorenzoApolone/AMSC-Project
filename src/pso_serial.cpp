/**
 * @file pso_serial.cpp
 * @brief Serial implementation of the Particle Swarm Optimization (PSO) algorithm.
 *
 * This file defines the function `pso_serial`, which optimizes a given `TestFunction`
 * using a standard PSO scheme in a single-threaded (non-parallel) environment.
 *
 * The algorithm uses time-varying inertia weight, personal and global best updates,
 * and a termination criterion based on maximum iterations or convergence tolerance.
 *
 * @details
 * - Velocity update rule:
 *    v_ij ← w*v_ij + C1*r1*(pbest_ij - x_ij) + C2*r2*(gbest_j - x_ij)
 * - Position update rule:
 *    x_ij ← x_ij + v_ij
 * - Bound handling: positions are clamped within [LB, UB].
 *
 * @note The random seed is fixed (42) for reproducibility.
 */

#include "interfaces.hpp"
#include <vector>
#include <random>
#include <chrono>

/**
 * @struct PSOHyperparameters
 * @brief Constants and configuration parameters for the PSO algorithm
 */
struct PSOHyperparameters {
  /// @brief Cognitive coefficient (local best influence)
  static constexpr double C1 = 1.49618;
  
  /// @brief Social coefficient (global best influence)
  static constexpr double C2 = 1.49618;

  // Dynamic inertia configuration
  /// @brief Maximum inertia weight (start of exploration)
  static constexpr double W_MAX = 0.9;
  
  /// @brief Minimum inertia weight (end of exploitation)
  static constexpr double W_MIN = 0.4;

};

/**
 * @brief Serial Particle Swarm Optimization (PSO) algorithm.
 *
 * @param f        Objective function to minimize (derived from TestFunction).
 * @param d        Dimensionality of the search space.
 * @param stop     Stopping criterion (max iterations, tolerance, etc.).
 * @param n_points Number of particles in the swarm.
 * @return OutputObject containing optimization results, convergence history, etc.
 */

OutputObject pso_serial(const TestFunction& f, int d, const StopCriterion& stop, int n_points) {
    // ------------------------------
    // Initialization phase
    // ------------------------------
    
    auto start_time = std::chrono::high_resolution_clock::now();

    std::mt19937 gen(42);  
    std::pair<double, double> bounds = f.get_domain();
    double LB = bounds.first;
    double UB = bounds.second;
    std::uniform_real_distribution<double> dist(LB, UB);
    std::uniform_real_distribution<double> dist01(0.0, 1.0);

    // Particle data structures
    int local_n = n_points;  
    std::vector<std::vector<double>> pos(local_n, std::vector<double>(d));
    std::vector<std::vector<double>> vel(local_n, std::vector<double>(d));
    std::vector<std::vector<double>> pbest_pos = pos; 
    std::vector<double> pbest_val(local_n, std::numeric_limits<double>::max());

    std::vector<double> gbest_pos(d);
    double gbest_val = std::numeric_limits<double>::max();
    //double prev_gbest_val = std::numeric_limits<double>::max();

    // Best value history
    std::vector<double> history;
    history.reserve(stop.get_max_iter());

    // ------------------------------
    // Random initialization
    // ------------------------------
    for (int i = 0; i < local_n; ++i) {
        for (int j = 0; j < d; ++j) {
            // Initial position uniformly in [LB, UB]
            pos[i][j] = dist(gen);
            // Initial velocity as a fraction of the range (±10%)
            vel[i][j] = (dist(gen) - dist(gen)) * 0.1;
            pbest_pos[i][j] = pos[i][j];
        }
        // Initial fitness evaluation
        double fitness = f.value(pos[i]);
        pbest_val[i] = fitness;
        // Update initial global best
        if (fitness < gbest_val) {
            gbest_val = fitness;
            gbest_pos = pos[i];
        }
    }

    // ------------------------------
    // Main optimization loop
    // ------------------------------

    int iter = 0;
    bool must_stop = false;
    int max_iter_limit = stop.get_max_iter();

    while (!must_stop) {
        // Dynamic inertia weight calculation (linearly decreases)
        double current_w = PSOHyperparameters::W_MAX - ((PSOHyperparameters::W_MAX - PSOHyperparameters::W_MIN) * iter / max_iter_limit);

        // Update each particle
        for (int i = 0; i < local_n; ++i) {
            for (int j = 0; j < d; ++j) {
                double r1 = dist01(gen);
                double r2 = dist01(gen);
                // Update velocity
                vel[i][j] = current_w * vel[i][j]
                          + PSOHyperparameters::C1 * r1 * (pbest_pos[i][j] - pos[i][j])
                          + PSOHyperparameters::C2 * r2 * (gbest_pos[j] - pos[i][j]);
                // Update position
                pos[i][j] += vel[i][j];
                // Boundary check: avoid going out of bounds
                if (pos[i][j] < LB) pos[i][j] = LB;
                if (pos[i][j] > UB) pos[i][j] = UB;
            }
            // Current fitness evaluation
            double current_fit = f.value(pos[i]);
            // Update personal best
            if (current_fit < pbest_val[i]) {
                pbest_val[i] = current_fit;
                pbest_pos[i] = pos[i];
            }
            // Update global best
            if (pbest_val[i] < gbest_val) {
                gbest_val = pbest_val[i];
                gbest_pos = pbest_pos[i];
            }
        }

        

        // ------------------------------
        // Check stopping criteria
        // ------------------------------
        double current_normalized_error = f.error(gbest_pos);
        // Save error to history
        history.push_back(current_normalized_error);
        
        if (stop.should_stop(iter, current_normalized_error)) {
            must_stop = true;
        }
        // Prepare for next iteration
        //prev_gbest_val = gbest_val;
        iter++;
        if (iter >= max_iter_limit) {
            // Reached maximum number of iterations
            must_stop = true;
        }
    }

    // ------------------------------
    // Wrap-up phase
    // ------------------------------
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // Construct output structure
    OutputObject out(
        f.get_name(),             ///< Function name
        d,                        ///< Dimensionality
        n_points,                 ///< Number of particles
        gbest_pos,                ///< Best position found
        f.get_true_solution(),    ///< True (known) optimal solution
        gbest_val,                ///< Best fitness value
        history,                  ///< Convergence curve
        1,                        ///< Cores used (serial)
        elapsed.count(),          ///< Execution time [s]
        iter,                     ///< Number of iterations performed
        stop                      ///< Stop criterion used
    );
    return out;
}