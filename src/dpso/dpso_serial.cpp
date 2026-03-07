/**
 * @file dpso_serial.cpp
 * @brief Serial implementation of the Distributed Particle Swarm Optimization (DPSO) algorithm.
 *
 * This file implements a hybrid DPSO-Harmony Search metaheuristic for global optimization
 * in a single-threaded environment. The swarm is divided into fixed-size sub-swarms using
 * a local-best (lbest) topology. Particles are periodically regrouped (shuffled) to maintain
 * diversity, and Harmony Search is applied to each sub-swarm for local refinement.
 */

#include <vector>
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <chrono>
#include "particle.hpp"
#include "methods_dpso.hpp"
#include "interfaces.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"

/**
 * @brief Generates a random double in [min, max].
 * @param min Minimum value.
 * @param max Maximum value.
 * @return Random double in [min, max].
 */
static double random_double_serial(double min, double max) {
    static std::mt19937 gen(12345); 
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

/**
 * @brief Calculates the Euclidean distance between two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Euclidean distance between v1 and v2.
 */
static double euclidean_dist_serial(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

/**
 * @brief Applies Harmony Search to a subpopulation of particles.
 * @param swarm Vector of particles.
 * @param start_idx Start index of the subpopulation.
 * @param end_idx End index of the subpopulation.
 * @param f Test function to optimize.
 * @param lower_bound Lower bounds for each dimension.
 * @param upper_bound Upper bounds for each dimension.
 * @param current_iter Current iteration.
 * @param max_iter Maximum number of iterations.
 */
static void apply_harmony_search_serial(std::vector<Particle>& swarm, 
                          int start_idx, 
                          int end_idx, 
                          const TestFunction& f, 
                          const std::vector<double>& lower_bound,
                          const std::vector<double>& upper_bound,
                          int current_iter,
                          int max_iter) {
    int dim = lower_bound.size();
    int sub_swarm_size = end_idx - start_idx;
    if (sub_swarm_size <= 0) return; 
    const double HMCR = 0.98;
    const double PAR_min = 0.01;
    const double PAR_max = 0.99;
    double PAR = PAR_min + ((PAR_max - PAR_min) / max_iter) * current_iter;
    std::vector<double> new_harmony(dim);
    for (int d = 0; d < dim; ++d) {
        double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
        double bw_min = 0.0001;
        double bw = bw_max * std::exp((std::log(bw_min/bw_max) / max_iter) * current_iter);
        if (random_double_serial(0.0, 1.0) < HMCR) {
            if (sub_swarm_size > 0) {
                int random_member_idx = start_idx + (int)random_double_serial(0, sub_swarm_size - 0.001);
                if (random_member_idx < start_idx) random_member_idx = start_idx;
                if (random_member_idx >= end_idx) random_member_idx = end_idx - 1;
                new_harmony[d] = swarm[random_member_idx].best_position[d];
                if (random_double_serial(0.0, 1.0) < PAR) {
                    new_harmony[d] += random_double_serial(-1.0, 1.0) * bw;
                }
            } else {
                new_harmony[d] = random_double_serial(lower_bound[d], upper_bound[d]);
            }
        } else {
            new_harmony[d] = random_double_serial(lower_bound[d], upper_bound[d]);
        }
        if (new_harmony[d] < lower_bound[d]) new_harmony[d] = lower_bound[d];
        if (new_harmony[d] > upper_bound[d]) new_harmony[d] = upper_bound[d];
    }
    double new_harmony_val = f.value(new_harmony);
    int nearest_idx = -1;
    double min_dist = std::numeric_limits<double>::max();
    for (int i = start_idx; i < end_idx; ++i) {
        double d = euclidean_dist_serial(new_harmony, swarm[i].best_position);
        if (d < min_dist) {
            min_dist = d;
            nearest_idx = i;
        }
    }
    if (nearest_idx != -1 && new_harmony_val < swarm[nearest_idx].best_value) {
        swarm[nearest_idx].best_position = new_harmony;
        swarm[nearest_idx].best_value = new_harmony_val;
    }
}

/**
 * @brief Shuffles particles in-place to regroup sub-swarms.
 * @param swarm Vector of particles to shuffle.
 */
static void regroup_particles_serial(std::vector<Particle>& swarm) {
    static std::mt19937 g(12345);
    std::shuffle(swarm.begin(), swarm.end(), g);
}

/**
 * @brief Serial DPSO (Distributed Particle Swarm Optimization) algorithm.
 * @param f Test function to optimize.
 * @param dim Number of dimensions.
 * @param n_points_total Total number of particles.
 * @param max_iter Maximum number of iterations.
 * @return OutputObject with optimization results.
 */
OutputObject dpso_serial(const TestFunction& f, 
                 unsigned int dim, 
                 unsigned int n_points_total, 
                 int max_iter) {
    StoppingCriteriaManager stop_manager(max_iter, 500, 1e-8, 1e-3);

    const double w = 0.729; 
    const double c1 = 1.49445;
    const double c2 = 1.49445;
    const int regrouping_period = 5; 
    const int sub_swarm_size = 5; 

    if (n_points_total < (unsigned int)sub_swarm_size) {
        std::cerr << "Error: Total particles (" << n_points_total 
                  << ") less than sub-swarm size (" << sub_swarm_size << ")." << std::endl;
        return OutputObject(f.get_name(), dim, n_points_total, {}, f.get_true_solution(), 0.0, {}, 1, 0.0, 0, stop_manager);
    }
    if (n_points_total % sub_swarm_size != 0) {
        std::cerr << "Warning: Total particles (" << n_points_total 
                  << ") not divisible by sub-swarm size (" << sub_swarm_size << ")." << std::endl;
    }

    const auto& domain = f.get_domain();
    std::vector<double> lb(dim, domain.first);
    std::vector<double> ub(dim, domain.second);
    std::vector<double> v_max(dim);
    if (lb[0] > ub[0]) {
        for(unsigned int d=0; d<dim; ++d) {
            lb[d] = -ub[d];
        }
    }
    for(unsigned int d=0; d<dim; ++d) {
        v_max[d] = 0.2 * (ub[d] - lb[d]);
    }

    std::vector<Particle> swarm;
    swarm.reserve(n_points_total);
    for (unsigned int i = 0; i < n_points_total; ++i) {
        Particle p(dim);
        for (unsigned int d = 0; d < dim; ++d) {
            p.position[d] = random_double_serial(lb[d], ub[d]);
            p.velocity[d] = random_double_serial(-v_max[d], v_max[d]);
            p.best_position[d] = p.position[d];
        }
        p.current_value = f.value(p.position);
        p.best_value = p.current_value;
        swarm.push_back(p);
    }

    OutputObject results(f.get_name(), dim, n_points_total,
                         {}, f.get_true_solution(), 0.0, {}, 1, 0.0, 0, stop_manager);
    results.x_best.resize(dim);

    double global_best_val = std::numeric_limits<double>::max();
    auto start_time = std::chrono::high_resolution_clock::now();

    int iter = 0;
    while (true) {
        if (iter > 0 && iter % regrouping_period == 0) {
             regroup_particles_serial(swarm);
        }

        int num_sub_swarms = swarm.size() / sub_swarm_size;
        int remainder = swarm.size() % sub_swarm_size;

        for (int s = 0; s < num_sub_swarms; ++s) {
            int start = s * sub_swarm_size;
            int end = start + sub_swarm_size;
            if (end > (int)swarm.size()) end = swarm.size();

            int lbest_idx = -1;
            double lbest_val = std::numeric_limits<double>::max();
            for (int i = start; i < end; ++i) {
                if (swarm[i].best_value < lbest_val) {
                    lbest_val = swarm[i].best_value;
                    lbest_idx = i;
                }
            }
            if (lbest_idx == -1) continue; 
            std::vector<double> lbest_pos = swarm[lbest_idx].best_position;

            for (int i = start; i < end; ++i) {
                Particle& p = swarm[i];
                bool in_bounds = true;
                for (unsigned int d = 0; d < dim; ++d) {
                    double r1 = random_double_serial(0.0, 1.0);
                    double r2 = random_double_serial(0.0, 1.0);
                    p.velocity[d] = w * p.velocity[d] +
                                    c1 * r1 * (p.best_position[d] - p.position[d]) +
                                    c2 * r2 * (lbest_pos[d] - p.position[d]);
                    if (p.velocity[d] > v_max[d]) p.velocity[d] = v_max[d];
                    if (p.velocity[d] < -v_max[d]) p.velocity[d] = -v_max[d];
                    p.position[d] += p.velocity[d];
                    if (p.position[d] < lb[d] || p.position[d] > ub[d]) {
                        in_bounds = false;
                    }
                }
                if (in_bounds) {
                    p.current_value = f.value(p.position);
                    if (p.current_value < p.best_value) {
                        p.best_value = p.current_value;
                        p.best_position = p.position;
                    }
                }
            }
            apply_harmony_search_serial(swarm, start, end, f, lb, ub, iter, max_iter);
        }

        // Handle remainder particles (if any)
        if (remainder > 0) {
            int start = num_sub_swarms * sub_swarm_size;
            int end = swarm.size();
            if (end > start) {
                int lbest_idx = -1;
                double lbest_val = std::numeric_limits<double>::max();
                for (int i = start; i < end; ++i) {
                    if (swarm[i].best_value < lbest_val) {
                        lbest_val = swarm[i].best_value;
                        lbest_idx = i;
                    }
                }
                if (lbest_idx != -1) {
                    std::vector<double> lbest_pos = swarm[lbest_idx].best_position;
                    for (int i = start; i < end; ++i) {
                        Particle& p = swarm[i];
                        bool in_bounds = true;
                        for (unsigned int d = 0; d < dim; ++d) {
                            double r1 = random_double_serial(0.0, 1.0);
                            double r2 = random_double_serial(0.0, 1.0);
                            p.velocity[d] = w * p.velocity[d] +
                                            c1 * r1 * (p.best_position[d] - p.position[d]) +
                                            c2 * r2 * (lbest_pos[d] - p.position[d]);
                            if (p.velocity[d] > v_max[d]) p.velocity[d] = v_max[d];
                            if (p.velocity[d] < -v_max[d]) p.velocity[d] = -v_max[d];
                            p.position[d] += p.velocity[d];
                            if (p.position[d] < lb[d] || p.position[d] > ub[d]) {
                                in_bounds = false;
                            }
                        }
                        if (in_bounds) {
                            p.current_value = f.value(p.position);
                            if (p.current_value < p.best_value) {
                                p.best_value = p.current_value;
                                p.best_position = p.position;
                            }
                        }
                    }
                    apply_harmony_search_serial(swarm, start, end, f, lb, ub, iter, max_iter);
                }
            }
        }

        /// @name Global best detection across all particles
        /// @{
        double current_global_min = std::numeric_limits<double>::max();
        int best_idx = -1;
        for(int i = 0; i < (int)swarm.size(); ++i) {
            if (swarm[i].best_value < current_global_min) {
                current_global_min = swarm[i].best_value;
                best_idx = i;
            }
        }

        std::vector<double> global_best_position(dim);
        if (best_idx != -1) {
            global_best_position = swarm[best_idx].best_position;
        }
        /// @}

        /// @name Swarm diversity: average distance to global best
        /// @{
        double sum_dist = 0.0;
        for(const auto& p : swarm) {
            sum_dist += euclidean_dist_serial(p.position, global_best_position);
        }
        double avg_dist = swarm.size() > 0 ? sum_dist / swarm.size() : 0.0;
        /// @}

        results.conv_history.push_back(current_global_min);
        global_best_val = current_global_min;

        stop_manager.increment_iterations();
        iter++;

        /// Stopping criteria check (max iterations, stagnation, or convergence)
        if (stop_manager.should_stop(global_best_val, avg_dist)) {
            break;
        }
    }

    /// @name Extract final best solution
    /// @{
    double best_val = std::numeric_limits<double>::max();
    int best_idx_final = -1;
    for(int i = 0; i < (int)swarm.size(); ++i) {
        if(swarm[i].best_value < best_val) {
            best_val = swarm[i].best_value;
            best_idx_final = i;
        }
    }
    if (best_idx_final != -1) {
        results.x_best = swarm[best_idx_final].best_position;
        results.f_val = best_val;
    }
    /// @}

    auto end_time = std::chrono::high_resolution_clock::now();
    results.execution_time = std::chrono::duration<double>(end_time - start_time).count();
    results.iterations = iter;
    return results;
}