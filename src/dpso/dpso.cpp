#include "particle.hpp"
#include "methods_dpso.hpp"
#include "interfaces.hpp"
#include <mpi.h>
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>


double random_double(double min, double max, int rank) {
    static std::mt19937 gen(rank * 10000 + 12345); // Rimosso omp_get_thread_num e thread_local
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

double euclidean_dist(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
        sum = std::abs(sum);
    }
    return std::sqrt(sum);
}

void apply_harmony_search(std::vector<Particle>& swarm, 
                          int start_idx, 
                          int end_idx, 
                          const TestFunction& f, 
                          int rank,
                          const std::vector<double>& lower_bound,
                          const std::vector<double>& upper_bound,
                          int current_iter,
                          int max_iter) {
    int dim = lower_bound.size();
    int sub_swarm_size = end_idx - start_idx;
    if (sub_swarm_size <= 0) return; // Prevent invalid array sizes
    const double HMCR = 0.98;
    const double PAR_min = 0.01;
    const double PAR_max = 0.99;
    double PAR = PAR_min + ((PAR_max - PAR_min) / max_iter) * current_iter;
    std::vector<double> new_harmony(dim);
    for (int d = 0; d < dim; ++d) {
        double bw_max = 0.05 * (upper_bound[d] - lower_bound[d]);
        double bw_min = 0.0001;
        double bw = bw_max * std::exp((std::log(bw_min/bw_max) / max_iter) * current_iter);
        if (random_double(0.0, 1.0, rank) < HMCR) {
            if (sub_swarm_size > 0) {
                int random_member_idx = start_idx + (int)random_double(0, sub_swarm_size - 0.001, rank);
                // Clamp index to valid range
                if (random_member_idx < start_idx) random_member_idx = start_idx;
                if (random_member_idx >= end_idx) random_member_idx = end_idx - 1;
                new_harmony[d] = swarm[random_member_idx].best_position[d];
                if (random_double(0.0, 1.0, rank) < PAR) {
                    new_harmony[d] += random_double(-1.0, 1.0, rank) * bw;
                }
            } else {
                new_harmony[d] = random_double(lower_bound[d], upper_bound[d], rank);
            }
        } else {
            new_harmony[d] = random_double(lower_bound[d], upper_bound[d], rank);
        }
        if (new_harmony[d] < lower_bound[d]) new_harmony[d] = lower_bound[d];
        if (new_harmony[d] > upper_bound[d]) new_harmony[d] = upper_bound[d];
    }
    double new_harmony_val = f.value(new_harmony);
    int nearest_idx = -1;
    double min_dist = std::numeric_limits<double>::max();
    for (int i = start_idx; i < end_idx; ++i) {
        double d = euclidean_dist(new_harmony, swarm[i].best_position);
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

void regroup_particles(std::vector<Particle>& local_swarm, int dim, int rank, int size) {
    int local_n = local_swarm.size();
    int p_data_size = 3 * dim + 2; 
    std::vector<double> send_buffer;
    send_buffer.reserve(local_n * p_data_size);
    for (const auto& p : local_swarm) {
        send_buffer.insert(send_buffer.end(), p.position.begin(), p.position.end());
        send_buffer.insert(send_buffer.end(), p.velocity.begin(), p.velocity.end());
        send_buffer.insert(send_buffer.end(), p.best_position.begin(), p.best_position.end());
        send_buffer.push_back(p.best_value);
        send_buffer.push_back(p.current_value);
    }
    std::vector<double> recv_buffer(local_n * size * p_data_size);
    MPI_Allgather(send_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                  recv_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                  MPI_COMM_WORLD);
    std::vector<int> indices(local_n * size);
    std::iota(indices.begin(), indices.end(), 0);
    static std::mt19937 g(12345);
    std::shuffle(indices.begin(), indices.end(), g);
    int global_idx_start = rank * local_n;
    for (int i = 0; i < local_n; ++i) {
        int picked_idx = indices[global_idx_start + i];
        int base = picked_idx * p_data_size;
        int offset = 0;
        std::copy(recv_buffer.begin() + base + offset, 
                  recv_buffer.begin() + base + offset + dim, 
                  local_swarm[i].position.begin()); offset += dim;
        std::copy(recv_buffer.begin() + base + offset, 
                  recv_buffer.begin() + base + offset + dim, 
                  local_swarm[i].velocity.begin()); offset += dim;
        std::copy(recv_buffer.begin() + base + offset, 
                  recv_buffer.begin() + base + offset + dim, 
                  local_swarm[i].best_position.begin()); offset += dim;
        local_swarm[i].best_value = recv_buffer[base + offset++];
        local_swarm[i].current_value = recv_buffer[base + offset++];
    }
}

OutputObject pso_mpi(const TestFunction& f, 
                     unsigned int dim, 
                     const StopCriterion& stop, 
                     unsigned int n_points_total) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Divide total particles by number of ranks
    unsigned int n_points_per_rank = n_points_total / size;
    if (rank == 0 && n_points_total % size != 0) {
        std::cerr << "Warning: total particles (" << n_points_total << ") not divisible by number of ranks (" << size << ")." << std::endl;
    }
    const double w = 0.729; 
    const double c1 = 1.49445;
    const double c2 = 1.49445;
    const int regrouping_period = 10; 
    const int sub_swarm_size = 5; 
    if (n_points_per_rank < sub_swarm_size && rank == 0) {
        std::cerr << "Error: Particles per rank (" << n_points_per_rank 
                  << ") less than sub-swarm size (" << sub_swarm_size << ")." << std::endl;
        return OutputObject(f.get_name(), dim, n_points_per_rank * size, {}, f.get_true_solution(), 0.0, {}, size, 0.0, 0, stop);
    }
    if (n_points_per_rank % sub_swarm_size != 0 && rank == 0) {
        std::cerr << "Warning: Particles per rank (" << n_points_per_rank 
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
    swarm.reserve(n_points_per_rank);
    for (unsigned int i = 0; i < n_points_per_rank; ++i) {
        Particle p(dim);
        for (unsigned int d = 0; d < dim; ++d) {
            p.position[d] = random_double(lb[d], ub[d], rank);
            p.velocity[d] = random_double(-v_max[d], v_max[d], rank); // Ensure velocity is initialized
            p.best_position[d] = p.position[d]; // Ensure best_position is initialized
        }
        p.current_value = f.value(p.position);
        p.best_value = p.current_value;
        swarm.push_back(p);
    }
    OutputObject results(f.get_name(), dim, n_points_per_rank * size,
                         {}, f.get_true_solution(), 0.0, {}, size, 0.0, 0, stop);
    results.x_best.resize(dim);
    double global_best_val = std::numeric_limits<double>::max();
    double start_time = MPI_Wtime();
    int iter = 0;
    int max_iter = stop.get_max_iter();
    while (!stop.should_stop(iter, global_best_val)) {
        if (iter > 0 && iter % regrouping_period == 0) {
             regroup_particles(swarm, dim, rank, size);
        }
        int num_sub_swarms = swarm.size() / sub_swarm_size;
        int remainder = swarm.size() % sub_swarm_size;
        // Rimosso #pragma omp parallel for
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
            if (lbest_idx == -1) continue; // skip if no valid lbest
            std::vector<double> lbest_pos = swarm[lbest_idx].best_position;
            for (int i = start; i < end; ++i) {
                Particle& p = swarm[i];
                bool in_bounds = true;
                for (unsigned int d = 0; d < dim; ++d) {
                    double r1 = random_double(0.0, 1.0, rank);
                    double r2 = random_double(0.0, 1.0, rank);
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
            apply_harmony_search(swarm, start, end, f, rank, lb, ub, iter, max_iter);
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
                            double r1 = random_double(0.0, 1.0, rank);
                            double r2 = random_double(0.0, 1.0, rank);
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
                    apply_harmony_search(swarm, start, end, f, rank, lb, ub, iter, max_iter);
                }
            }
        }
        double local_min_val = std::numeric_limits<double>::max();
        for(const auto& p : swarm) {
            if (p.best_value < local_min_val) local_min_val = p.best_value;
        }
        double current_global_min;
        MPI_Allreduce(&local_min_val, &current_global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (rank == 0) {
            results.conv_history.push_back(current_global_min);
            global_best_val = current_global_min;
        }
        MPI_Bcast(&global_best_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        iter++;
    }
    struct { double val; int rank; } loc_data, glob_data;
    loc_data.val = std::numeric_limits<double>::max();
    int best_idx_local = -1;
    for(int i=0; i < (int)swarm.size(); ++i) {
        if(swarm[i].best_value < loc_data.val) {
            loc_data.val = swarm[i].best_value;
            best_idx_local = i;
        }
    }
    loc_data.rank = rank;
    MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if (rank == glob_data.rank) {
        results.x_best = swarm[best_idx_local].best_position;
        results.f_val = glob_data.val;
    }
    MPI_Bcast(results.x_best.data(), dim, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);
    results.execution_time = MPI_Wtime() - start_time;
    results.iterations = iter;
    return results;
}