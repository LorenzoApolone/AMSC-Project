#include "methods_dpso.hpp"
#include "interfaces.hpp"
#include <mpi.h>
#include <omp.h>
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>

Particle::Particle(int dim) {
    position.resize(dim);
    velocity.resize(dim, 0.0);
    p_best_pos.resize(dim);
    p_best_val = std::numeric_limits<double>::max();
}

double random_double(double min, double max, int rank) {
    thread_local std::mt19937 gen(rank * 1000 + omp_get_thread_num() + 12345); 
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

void apply_harmony_search(std::vector<Particle>& swarm, 
                          int candidate_idx, 
                          const TestFunction& f, 
                          int rank,
                          const std::vector<double>& lower_bound,
                          const std::vector<double>& upper_bound) {
    const double HMCR = 0.95;  
    const double PAR = 0.1;    
    const double BW = 0.01;   
    int dim = swarm[candidate_idx].position.size();
    std::vector<double> new_pos = swarm[candidate_idx].position; 
    for (int d = 0; d < dim; ++d) {
        double r1 = random_double(0.0, 1.0, rank);
        if (r1 < HMCR) {
            double new_val = new_pos[d]; 
            double r2 = random_double(0.0, 1.0, rank);
            if (r2 < PAR) {
                double adjustment = random_double(-1.0, 1.0, rank) * BW;
                new_val += adjustment;
            }
            new_pos[d] = new_val;
        } else {
            new_pos[d] = random_double(lower_bound[d], upper_bound[d], rank);
        }
        if (new_pos[d] < lower_bound[d]) new_pos[d] = lower_bound[d];
        if (new_pos[d] > upper_bound[d]) new_pos[d] = upper_bound[d];
    }
    double new_val = f.value(new_pos);
    if (new_val < swarm[candidate_idx].p_best_val) {
        swarm[candidate_idx].position = new_pos;
        swarm[candidate_idx].current_val = new_val;
        swarm[candidate_idx].p_best_pos = new_pos;
        swarm[candidate_idx].p_best_val = new_val;
    }
}

void regroup_particles(std::vector<Particle>& local_swarm, int dim, int rank, int size) {
    int local_n = local_swarm.size();
    int p_size = 3 * dim + 2; 
    std::vector<double> send_buffer;
    send_buffer.reserve(local_n * p_size);
    for (const auto& p : local_swarm) {
        send_buffer.insert(send_buffer.end(), p.position.begin(), p.position.end());
        send_buffer.insert(send_buffer.end(), p.velocity.begin(), p.velocity.end());
        send_buffer.insert(send_buffer.end(), p.p_best_pos.begin(), p.p_best_pos.end());
        send_buffer.push_back(p.p_best_val);
        send_buffer.push_back(p.current_val);
    }
    std::vector<double> recv_buffer;
    if (rank == 0) {
        recv_buffer.resize(local_n * size * p_size);
    }
    MPI_Gather(send_buffer.data(), send_buffer.size(), MPI_DOUBLE,
               recv_buffer.data(), send_buffer.size(), MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::vector<int> indices(local_n * size);
        std::iota(indices.begin(), indices.end(), 0);
        static std::mt19937 g(42);
        std::shuffle(indices.begin(), indices.end(), g);
        std::vector<double> shuffled_buffer = recv_buffer;
        for (int i = 0; i < indices.size(); ++i) {
            int src_idx = indices[i];
            int dest_idx = i;
            std::copy(recv_buffer.begin() + src_idx * p_size, 
                      recv_buffer.begin() + (src_idx + 1) * p_size,
                      shuffled_buffer.begin() + dest_idx * p_size);
        }
        recv_buffer = shuffled_buffer;
    }
    MPI_Scatter(recv_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                send_buffer.data(), send_buffer.size(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    int idx = 0;
    for (auto& p : local_swarm) {
        std::copy(send_buffer.begin() + idx, send_buffer.begin() + idx + dim, p.position.begin()); idx += dim;
        std::copy(send_buffer.begin() + idx, send_buffer.begin() + idx + dim, p.velocity.begin()); idx += dim;
        std::copy(send_buffer.begin() + idx, send_buffer.begin() + idx + dim, p.p_best_pos.begin()); idx += dim;
        p.p_best_val = send_buffer[idx++];
        p.current_val = send_buffer[idx++];
    }
}

OutputObject pso_mpi(const TestFunction& f, 
                     unsigned int dim, 
                     const StopCriterion& stop, 
                     unsigned int n_points_per_rank) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const double w = 0.729; 
    const double c1 = 1.49;
    const double c2 = 1.49;
    const int regrouping_period = 10;
    const auto& domain = f.get_domain();
    std::vector<double> lb(dim, domain.first);
    std::vector<double> ub(dim, domain.second);
    std::vector<Particle> swarm;
    swarm.reserve(n_points_per_rank);
    for (unsigned int i = 0; i < n_points_per_rank; ++i) {
        Particle p(dim);
        for (unsigned int d = 0; d < dim; ++d) {
            p.position[d] = random_double(lb[d], ub[d], rank);
            p.velocity[d] = random_double(-(ub[d]-lb[d])*0.1, (ub[d]-lb[d])*0.1, rank);
            p.p_best_pos[d] = p.position[d];
        }
        p.current_val = f.value(p.position);
        p.p_best_val = p.current_val;
        swarm.push_back(p);
    }
    OutputObject results(f.get_name(), dim, n_points_per_rank * size,
                         {}, f.get_true_solution(), 0.0, {}, size, 0.0, 0, stop);
    results.x_best.resize(dim);
    double global_best_val = std::numeric_limits<double>::max();
    double local_best_val = std::numeric_limits<double>::max();
    std::vector<double> local_best_pos(dim);
    int local_best_idx = -1;
    for(int i=0; i < swarm.size(); ++i) {
        if(swarm[i].p_best_val < local_best_val) {
            local_best_val = swarm[i].p_best_val;
            local_best_pos = swarm[i].p_best_pos;
            local_best_idx = i;
        }
    }
    double start_time = MPI_Wtime();
    unsigned int iter = 0;
    while (iter < (unsigned int)stop.get_max_iter()) { 
        #pragma omp parallel for
        for (int i=0; i < swarm.size(); ++i) {
            Particle& p = swarm[i];
            for (unsigned int d = 0; d < dim; ++d) {
                double r1 = random_double(0.0, 1.0, rank);
                double r2 = random_double(0.0, 1.0, rank);
                p.velocity[d] = w * p.velocity[d] +
                                c1 * r1 * (p.p_best_pos[d] - p.position[d]) +
                                c2 * r2 * (local_best_pos[d] - p.position[d]);
                p.position[d] += p.velocity[d];
                if(p.position[d] < lb[d]) { p.position[d] = lb[d]; p.velocity[d] = 0; }
                if(p.position[d] > ub[d]) { p.position[d] = ub[d]; p.velocity[d] = 0; }
            }
            p.current_val = f.value(p.position);
            if (p.current_val < p.p_best_val) {
                p.p_best_val = p.current_val;
                p.p_best_pos = p.position;
            }
        }
        for (int i=0; i < swarm.size(); ++i) {
            if (swarm[i].current_val < local_best_val) {
                local_best_val = swarm[i].current_val;
                local_best_pos = swarm[i].position;
                local_best_idx = i;
            }
        }
        if (local_best_idx != -1) {
            apply_harmony_search(swarm, local_best_idx, f, rank, lb, ub);
            if (swarm[local_best_idx].p_best_val < local_best_val) {
                local_best_val = swarm[local_best_idx].p_best_val;
                local_best_pos = swarm[local_best_idx].p_best_pos;
            }
        }
        double current_global_min;
        MPI_Allreduce(&local_best_val, &current_global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        if (rank == 0) results.conv_history.push_back(current_global_min);
        if (iter > 0 && iter % regrouping_period == 0) {
             regroup_particles(swarm, dim, rank, size);
             local_best_val = std::numeric_limits<double>::max();
             local_best_idx = -1;
             for(int i=0; i < swarm.size(); ++i) {
                 if(swarm[i].p_best_val < local_best_val) {
                     local_best_val = swarm[i].p_best_val;
                     local_best_pos = swarm[i].p_best_pos;
                     local_best_idx = i;
                 }
             }
        }
        iter++;
    }
    struct { double val; int rank; } loc_data, glob_data;
    loc_data.val = local_best_val;
    loc_data.rank = rank;
    MPI_Allreduce(&loc_data, &glob_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if (rank == glob_data.rank) {
        results.x_best = local_best_pos;
    }
    MPI_Bcast(results.x_best.data(), dim, MPI_DOUBLE, glob_data.rank, MPI_COMM_WORLD);
    results.execution_time = MPI_Wtime() - start_time;
    return results;
}