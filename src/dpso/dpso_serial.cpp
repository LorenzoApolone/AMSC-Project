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
 * @brief Generates a random double within [min, max) using a fixed seed.
 * 
 * @param min Minimum bound.
 * @param max Maximum bound.
 * @return A random double.
 */
static double random_double_serial(double min, double max) {
    static std::mt19937 gen(12345);
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

/**
 * @brief Generates a random integer within [min, max] using a fixed seed.
 * 
 * Avoids precision issues resulting from casting double variables.
 * 
 * @param min Minimum integer (inclusive).
 * @param max Maximum integer (inclusive).
 * @return A random integer.
 */
static int random_int_serial(int min, int max) {
    static std::mt19937 gen(54321);
    std::uniform_int_distribution<> dis(min, max);
    return dis(gen);
}

static double euclidean_dist_serial(const std::vector<double>& v1, const std::vector<double>& v2) {
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

    const double HMCR    = 0.98;
    const double PAR_min = 0.01;
    const double PAR_max = 0.99;
    double PAR = PAR_min + ((PAR_max - PAR_min) / max_iter) * current_iter;

    std::vector<double> new_harmony(dim);
    double iter_ratio = (double)current_iter / max_iter;

    for (int d = 0; d < dim; ++d) {
        double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
        double bw_min = 0.0001;
        // Optimization: minimize log calls by precomputing operands
        double bw = bw_max * std::exp(std::log(bw_min / bw_max) * iter_ratio);

        if (random_double_serial(0.0, 1.0) < HMCR) {
            int idx = start_idx + random_int_serial(0, sub_swarm_size - 1);
            new_harmony[d] = swarm[idx].best_position[d];
            if (random_double_serial(0.0, 1.0) < PAR)
                new_harmony[d] += random_double_serial(-1.0, 1.0) * bw;
        } else {
            new_harmony[d] = random_double_serial(lower_bound[d], upper_bound[d]);
        }
        new_harmony[d] = std::max(lower_bound[d], std::min(upper_bound[d], new_harmony[d]));
    }

    double new_val = f.value(new_harmony);
    int nearest_idx = -1;
    double min_dist = std::numeric_limits<double>::max();
    for (int i = start_idx; i < end_idx; ++i) {
        double d = euclidean_dist_serial(new_harmony, swarm[i].best_position);
        if (d < min_dist) { min_dist = d; nearest_idx = i; }
    }
    if (nearest_idx != -1 && new_val < swarm[nearest_idx].best_value) {
        swarm[nearest_idx].best_position = new_harmony;
        swarm[nearest_idx].best_value    = new_val;
    }
}

/**
 * @brief Executes one local best (lbest) PSO + HS iteration on a sub-swarm.
 * 
 * Used for both complete sub-swarms and the possible remainder group.
 * 
 * @param swarm Reference to the entire swarm.
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
 * @param iter Current iteration number.
 * @param max_iter Maximum number of iterations.
 */
static void process_sub_swarm_serial(std::vector<Particle>& swarm,
                                     int start, int end,
                                     unsigned int dim,
                                     const std::vector<double>& lb,
                                     const std::vector<double>& ub,
                                     const std::vector<double>& v_max,
                                     double w, double c1, double c2,
                                     const TestFunction& f,
                                     int iter, int max_iter) {
    int lbest_idx = -1;
    double lbest_val = std::numeric_limits<double>::max();
    for (int i = start; i < end; ++i) {
        if (swarm[i].best_value < lbest_val) {
            lbest_val = swarm[i].best_value;
            lbest_idx = i;
        }
    }
    if (lbest_idx == -1) return;
    std::vector<double> lbest_pos = swarm[lbest_idx].best_position;

    for (int i = start; i < end; ++i) {
        Particle& p = swarm[i];
        for (unsigned int d = 0; d < dim; ++d) {
            double r1 = random_double_serial(0.0, 1.0);
            double r2 = random_double_serial(0.0, 1.0);
            p.velocity[d] = w * p.velocity[d]
                          + c1 * r1 * (p.best_position[d] - p.position[d])
                          + c2 * r2 * (lbest_pos[d]        - p.position[d]);
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
            p.best_value    = p.current_value;
            p.best_position = p.position;
        }
    }
    apply_harmony_search_serial(swarm, start, end, f, lb, ub, iter, max_iter);
}

static void regroup_particles_serial(std::vector<Particle>& swarm) {
    static std::mt19937 g(12345);
    std::shuffle(swarm.begin(), swarm.end(), g);
}

OutputObject dpso_serial(const TestFunction& f,
                         unsigned int dim,
                         unsigned int n_points_total,
                         int max_iter) {
    StoppingCriteriaManager stop_manager(max_iter, 500, 1e-8, 1e-3);

    // Parameters from Zhao et al. 2011, Tables 2-3
    const double w  = 0.729;
    const double c1 = 1.49445;
    const double c2 = 1.49445;
    const int regrouping_period = 5;
    const int sub_swarm_size    = 5;

    if (n_points_total < (unsigned int)sub_swarm_size) {
        std::cerr << "Error: total particles (" << n_points_total
                  << ") less than sub-swarm size (" << sub_swarm_size << ").\n";
        return OutputObject(f.get_name(), dim, n_points_total, {},
                            f.get_true_solution(), 0.0, {}, 1, 0.0, 0, stop_manager);
    }
    if (n_points_total % sub_swarm_size != 0)
        std::cerr << "Warning: total particles (" << n_points_total
                  << ") not divisible by sub-swarm size (" << sub_swarm_size << ").\n";

    const auto& domain = f.get_domain();
    std::vector<double> lb(dim, domain.first);
    std::vector<double> ub(dim, domain.second);
    std::vector<double> v_max(dim);
    if (lb[0] > ub[0])
        for (unsigned int d = 0; d < dim; ++d) lb[d] = -ub[d];
    for (unsigned int d = 0; d < dim; ++d)
        v_max[d] = 0.2 * (ub[d] - lb[d]);

    std::vector<Particle> swarm;
    swarm.reserve(n_points_total);
    for (unsigned int i = 0; i < n_points_total; ++i) {
        Particle p(dim);
        for (unsigned int d = 0; d < dim; ++d) {
            p.position[d]      = random_double_serial(lb[d],     ub[d]);
            p.velocity[d]      = random_double_serial(-v_max[d], v_max[d]);
            p.best_position[d] = p.position[d];
        }
        p.current_value = f.value(p.position);
        p.best_value    = p.current_value;
        swarm.push_back(p);
    }

    OutputObject results(f.get_name(), dim, n_points_total,
                         {}, f.get_true_solution(), 0.0, {}, 1, 0.0, 0, stop_manager);
    results.x_best.resize(dim);
    double global_best_val = std::numeric_limits<double>::max();
    auto start_time = std::chrono::high_resolution_clock::now();
    int iter = 0;

    while (true) {
        if (iter > 0 && iter % regrouping_period == 0)
            regroup_particles_serial(swarm);

        int num_sub_swarms = swarm.size() / sub_swarm_size;
        int remainder      = swarm.size() % sub_swarm_size;

        for (int s = 0; s < num_sub_swarms; ++s) {
            int start = s * sub_swarm_size;
            int end   = std::min(start + sub_swarm_size, (int)swarm.size());
            process_sub_swarm_serial(swarm, start, end, dim, lb, ub, v_max,
                                     w, c1, c2, f, iter, max_iter);
        }
        if (remainder > 0) {
            int start = num_sub_swarms * sub_swarm_size;
            process_sub_swarm_serial(swarm, start, (int)swarm.size(), dim, lb, ub, v_max,
                                     w, c1, c2, f, iter, max_iter);
        }

        double current_global_min = std::numeric_limits<double>::max();
        int best_idx = -1;
        for (int i = 0; i < (int)swarm.size(); ++i) {
            if (swarm[i].best_value < current_global_min) {
                current_global_min = swarm[i].best_value;
                best_idx = i;
            }
        }

        std::vector<double> global_best_pos(dim);
        if (best_idx != -1)
            global_best_pos = swarm[best_idx].best_position;

        double sum_dist = 0.0;
        for (const auto& p : swarm)
            sum_dist += euclidean_dist_serial(p.position, global_best_pos);
        double avg_dist = swarm.empty() ? 0.0 : sum_dist / swarm.size();

        results.conv_history.push_back(current_global_min);
        global_best_val = current_global_min;

        stop_manager.increment_iterations();
        iter++;

        if (stop_manager.should_stop(global_best_val, avg_dist))
            break;
    }

    double best_val = std::numeric_limits<double>::max();
    int best_idx_final = -1;
    for (int i = 0; i < (int)swarm.size(); ++i) {
        if (swarm[i].best_value < best_val) {
            best_val       = swarm[i].best_value;
            best_idx_final = i;
        }
    }
    if (best_idx_final != -1) {
        results.x_best = swarm[best_idx_final].best_position;
        results.f_val  = best_val;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    results.execution_time = std::chrono::duration<double>(end_time - start_time).count();
    results.iterations     = iter;
    return results;
}
