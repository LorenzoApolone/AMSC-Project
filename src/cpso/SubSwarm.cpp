/**
 * @file SubSwarm.cpp
 * @brief Implementation of the SubSwarm class methods
 */
#include "SubSwarm.hpp"
#include <stdexcept>

// Constuctor of the SubSwarm class
SubSwarm::SubSwarm(int num_particles, const std::vector<int> &active_dimensions, double lower_bound, 
                   double upper_bound, const std::vector<std::vector<int>> &adj_list)
    : num_particles(num_particles), active_dims(active_dimensions),
      gbest_val(std::numeric_limits<double>::infinity()),
      gbest_particle_idx(-1), ignore_gbest_this_iter(false),
      adjacency_list(adj_list), bounds_lower(lower_bound),
      bounds_upper(upper_bound) {
  int dim = active_dims.size();
  gbest_pos.resize(dim, 0.0);
  
  positions.resize(num_particles * dim);
  velocities.resize(num_particles * dim);
  best_positions.resize(num_particles * dim);
  current_values.resize(num_particles, std::numeric_limits<double>::infinity());
  best_values.resize(num_particles, std::numeric_limits<double>::infinity());
}


void SubSwarm::initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f) {
  
  // Creating probability distributions for both positions and velocity
  std::uniform_real_distribution<double> dist_pos(bounds_lower, bounds_upper);
  double range = bounds_upper - bounds_lower;

  std::uniform_real_distribution<double> dist_vel(-0.1 * range, 0.1 * range);

  // We get the context vector to extract the values of the active dimensions
  const std::vector<double> &context_vec = ctx.get_full_vector();

  // For every particle we initialize its position and its velocity
  // We set the best position as the initial one
  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {
    std::vector<double> temp_pos(dim);
    for (int d = 0; d < dim; ++d) {
      int idx = i * dim + d;
      if (i == 0) {
        // Place the first particle at the context vector without energy
        positions[idx] = context_vec[active_dims[d]];
        velocities[idx] = 0.0;
      } else {
        positions[idx] = dist_pos(gen);
        velocities[idx] = dist_vel(gen);
      }
      best_positions[idx] = positions[idx];
      temp_pos[d] = positions[idx];
    }

    // Evaluate the particle position
    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);
    best_values[i] = current_values[i];

    // Update the best position and value if the current position is better
    if (current_values[i] < gbest_val) {
      gbest_val = current_values[i];
      gbest_pos = temp_pos;
      gbest_particle_idx = i;
    }
  }
}

void SubSwarm::recalculate_fitness(ContextVector &ctx, const TestFunction &f) {
  gbest_val = std::numeric_limits<double>::infinity();
  gbest_particle_idx = -1;

  int dim = active_dims.size();
  std::vector<double> temp_best_pos(dim);
  std::vector<double> temp_pos(dim);

  for (int i = 0; i < num_particles; ++i) {
    for (int d = 0; d < dim; ++d) {
      int idx = i * dim + d;
      temp_best_pos[d] = best_positions[idx];
      temp_pos[d] = positions[idx];
    }

    // Re-evaluate the historical best position of the particle using the new context
    best_values[i] = ctx.evaluate_particle(f, temp_best_pos, active_dims);
    
    // Current position can also change its relative fitness compared to other particles
    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);

    // If the particle's current position in this new context is strictly better than its historical memory
    if (current_values[i] < best_values[i]) {
      best_values[i] = current_values[i];
      for (int d = 0; d < dim; ++d) {
          best_positions[i * dim + d] = positions[i * dim + d];
      }
      temp_best_pos = temp_pos;
    }

    if (best_values[i] < gbest_val) {
      gbest_val = best_values[i];
      gbest_pos = temp_best_pos;
      gbest_particle_idx = i;
    }
  }
}

void SubSwarm::update_active_dims(const std::vector<int> &new_dims,
                                  const ContextVector &ctx, std::mt19937 &gen, bool is_owned) {
  
  if (new_dims.size() != active_dims.size()) {
    throw std::runtime_error("Size mismatch during dimension shuffling");
  }

  // We overwrite the active dimensions with the new ones
  active_dims = new_dims;

  // Since the dims changed, the previous gbest is not valid anymore, so we set it to infinity
  gbest_val = std::numeric_limits<double>::infinity();
  gbest_particle_idx = -1;
  gbest_pos.assign(active_dims.size(), 0.0);

  // We get the context vector to extract the values of the new active dimensions
  const std::vector<double> &context_vec = ctx.get_full_vector();
  int dim = active_dims.size();

  if (!is_owned) {
    for (int i = 0; i < num_particles; ++i) {
      current_values[i] = std::numeric_limits<double>::infinity();
      best_values[i] = std::numeric_limits<double>::infinity();
      for (int d = 0; d < dim; ++d) {
        int idx = i * dim + d;
        positions[idx] = context_vec[active_dims[d]];
        best_positions[idx] = positions[idx];
        velocities[idx] = 0.0;
      }
    }
    return;
  }

  double range = bounds_upper - bounds_lower;

  // The particles will start around the global best, but we provide a wider spread to maintain diversity
  std::uniform_real_distribution<double> dist_pos(-0.25 * range, 0.25 * range);
  
  // Create the velocity distribution once outside the loops
  std::uniform_real_distribution<double> dist_vel(-0.05 * range, 0.05 * range);

  // For every particle we update its position and velocity
  for (int i = 0; i < num_particles; ++i) {
    current_values[i] = std::numeric_limits<double>::infinity();
    best_values[i] = std::numeric_limits<double>::infinity();

    for (int d = 0; d < dim; ++d) {
      int idx = i * dim + d;
      if (i == 0) {
        // Place the first particle at the context vector without energy
        positions[idx] = context_vec[active_dims[d]];
        velocities[idx] = 0.0;
      } else {
        positions[idx] = context_vec[active_dims[d]] + dist_pos(gen);
        // Inject some initial energy on the new dimensions
        velocities[idx] = dist_vel(gen);
      }
      
      // Keep within bounds
      if (positions[idx] < bounds_lower) positions[idx] = bounds_lower;
      if (positions[idx] > bounds_upper) positions[idx] = bounds_upper;

      best_positions[idx] = positions[idx];
    }
  }
}


void SubSwarm::inject_velocities(std::mt19937 &gen, bool hard_reset) {

  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_vel(-0.15 * range, 0.15 * range);

  // For every particle, we check if it is the global best. If not, we will inject random velocities
  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {
    // We use the particle index instead of floating-point comparison for robustness
    bool is_gbest = (i == gbest_particle_idx);

    // If the particle is not the global best, we inject random velocities and clear cognitive memory
    if (!is_gbest) {
      for (int d = 0; d < dim; ++d) {
        int idx = i * dim + d;

        // Double the original range for explosive escape
        velocities[idx] = dist_vel(gen) * 2.0; 
        
        if (hard_reset) {
          // Reset the personal memory so cognitive component doesn't pull it back instantly
          best_positions[idx] = positions[idx]; 
        }
      }
      if (hard_reset) {
        best_values[i] = current_values[i]; // Reset personal best fitness
      }
    }
  }
}

void SubSwarm::reset_gbest_attraction() {
  ignore_gbest_this_iter = true;
}

void SubSwarm::update_velocities_and_positions(double w, double c1, double c2,
                                               std::mt19937 &gen, double progress_ratio) {
  std::uniform_real_distribution<double> dist01(0.0, 1.0);


  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {
    
    std::vector<double> lbest_pos(dim);
    for (int d = 0; d < dim; ++d) {
        lbest_pos[d] = best_positions[i * dim + d];
    }
    double lbest_val = best_values[i];

    // If the adjacency list is not empty, we look for the best neighbor
    if (!adjacency_list.empty()) {
      for (int neighbor_idx : adjacency_list[i]) {
        if (best_values[neighbor_idx] < lbest_val) {
          
          // If the neighbor is better, we update the local best
          lbest_val = best_values[neighbor_idx];
          for (int d = 0; d < dim; ++d) {
              lbest_pos[d] = best_positions[neighbor_idx * dim + d];
          }
        }
      }
    } else {
      // If the adjacency list is empty, we use the global best
      lbest_pos = gbest_pos;
      lbest_val = gbest_val;
    }

    double v_max = (bounds_upper - bounds_lower) * (0.2 - 0.15 * progress_ratio);

    // We update the velocity and the position of every particle using the PSO equations
    for (int d = 0; d < dim; ++d) {
      double r1 = dist01(gen);
      double r2 = dist01(gen);
      
      int idx = i * dim + d;

      double c1_act = (best_values[i] == std::numeric_limits<double>::infinity()) ? 0.0 : c1;
      double c2_act = (lbest_val == std::numeric_limits<double>::infinity() || ignore_gbest_this_iter) ? 0.0 : c2;

      velocities[idx] = w * velocities[idx] +
                      c1_act * r1 * (best_positions[idx] - positions[idx]) +
                      c2_act * r2 * (lbest_pos[d] - positions[idx]);

      if (velocities[idx] > v_max) velocities[idx] = v_max;
      else if (velocities[idx] < -v_max) velocities[idx] = -v_max;

      positions[idx] += velocities[idx];

      // We apply the bounds-checking to the particle
      if (positions[idx] < bounds_lower) {
        positions[idx] = bounds_lower; 
        velocities[idx] = -0.5 * velocities[idx];
      } else if (positions[idx] > bounds_upper) {
        positions[idx] = bounds_upper;
        velocities[idx] = -0.5 * velocities[idx];
      }
    }
  }

  // Reset the temporary flag after the velocity update
  ignore_gbest_this_iter = false;
}

void SubSwarm::evaluate_and_update(ContextVector &ctx, const TestFunction &f) {

  int dim = active_dims.size();
  std::vector<double> temp_pos(dim);

  for (int i = 0; i < num_particles; ++i) {
    for (int d = 0; d < dim; ++d) {
        temp_pos[d] = positions[i * dim + d];
    }
    
    // Evaluates for every particles its current fitness
    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);

    // If the current position is better than the best position, we update the best position
    if (current_values[i] < best_values[i]) {
      best_values[i] = current_values[i];
      for (int d = 0; d < dim; ++d) {
          best_positions[i * dim + d] = positions[i * dim + d];
      }
    }

    // If the best position is better than the global best position, we update the global best position
    if (best_values[i] < gbest_val) {
      gbest_val = best_values[i];
      gbest_pos = temp_pos;
      gbest_particle_idx = i;
    }
  }
}

const std::vector<double> &SubSwarm::get_gbest_pos() const { return gbest_pos; }
double SubSwarm::get_gbest_val() const { return gbest_val; }
const std::vector<int> &SubSwarm::get_active_dims() const { return active_dims; }
int SubSwarm::get_num_particles() const { return num_particles; }
const std::vector<double> &SubSwarm::get_positions() const { return positions; }
