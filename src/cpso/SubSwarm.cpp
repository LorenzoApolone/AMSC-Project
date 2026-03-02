/**
 * @file SubSwarm.cpp
 * @brief Implementation of the SubSwarm class methods.
 */
#include "SubSwarm.hpp"
#include <stdexcept>

// Constuctor of the SubSwarm class.
SubSwarm::SubSwarm(int num_particles, const std::vector<int> &active_dimensions, double lower_bound, 
                   double upper_bound, const std::vector<std::vector<int>> &adj_list)
    : active_dims(active_dimensions),
      gbest_val(std::numeric_limits<double>::infinity()),
      adjacency_list(adj_list), bounds_lower(lower_bound),
      bounds_upper(upper_bound) {
  int dim = active_dims.size();
  gbest_pos.resize(dim, 0.0);
  particles.reserve(num_particles);
  for (int i = 0; i < num_particles; ++i) {
    particles.emplace_back(dim);
  }
}


void SubSwarm::initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f) {
  
  // Creating Probability distributions for both positions and velocity
  std::uniform_real_distribution<double> dist_pos(bounds_lower, bounds_upper);
  double range = bounds_upper - bounds_lower;

  std::uniform_real_distribution<double> dist_vel(-0.1 * range, 0.1 * range);

  // For every particle we initialize its position and its velocity
  // We set the best position as the initial one
  for (auto &p : particles) {
    for (size_t d = 0; d < active_dims.size(); ++d) {
      p.position[d] = dist_pos(gen);
      p.velocity[d] = dist_vel(gen);
      p.best_position[d] = p.position[d];
    }

    // Evaluate the particle position
    p.current_value = ctx.evaluate_particle(f, p.position, active_dims);
    p.best_value = p.current_value;

    // Update the best position and value if the current position is better
    if (p.current_value < gbest_val) {
      gbest_val = p.current_value;
      gbest_pos = p.position;
    }
  }
}

void SubSwarm::update_active_dims(const std::vector<int> &new_dims,
                                  const ContextVector &ctx, std::mt19937 &gen) {
  
  // Safety Check
  if (new_dims.size() != active_dims.size()) {
    throw std::runtime_error("Size mismatch during dimension shuffling");
  }

  // We overwrite the active dimensions with the new ones
  active_dims = new_dims;

  // Since the dims changed, the previous gbest is not valid anymore, so we set it to infinity
  gbest_val = std::numeric_limits<double>::infinity();

  // We get the context vector to extract the values of the new active dimensions
  const std::vector<double> &context_vec = ctx.get_full_vector();
  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_vel(-0.01 * range, 0.01 * range);

  // For every particle we update its position and velocity
  for (auto &p : particles) {
    p.current_value = std::numeric_limits<double>::infinity();
    p.best_value = std::numeric_limits<double>::infinity();

    for (size_t d = 0; d < active_dims.size(); ++d) {
      
      // The particles will starts where the were left previously, and not from the random initialization
      p.position[d] = context_vec[active_dims[d]];
      p.best_position[d] = p.position[d];

      p.velocity[d] += dist_vel(gen);
    }
  }
}


void SubSwarm::inject_velocities(std::mt19937 &gen) {

  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_vel(-0.1 * range, 0.1 * range);

  // For every particle, we check if it is the global best. If not, we will inject random velocities
  for (auto &p : particles) {
    bool is_gbest = true;
    for (size_t d = 0; d < active_dims.size(); ++d) {
      if (p.best_position[d] != gbest_pos[d]) {
        is_gbest = false;
        break;
      }
    }

    // If the particle is not the global best, we inject random velocities
    if (!is_gbest) {
      for (size_t d = 0; d < active_dims.size(); ++d) {
        p.velocity[d] = dist_vel(gen);
      }
    }
  }
}

void SubSwarm::update_velocities_and_positions(double w, double c1, double c2,
                                               std::mt19937 &gen) {
  std::uniform_real_distribution<double> dist01(0.0, 1.0);


  for (size_t i = 0; i < particles.size(); ++i) {
    auto &p = particles[i];

    std::vector<double> lbest_pos = p.best_position;
    double lbest_val = p.best_value;

    // If the adjacency list is not empty, we look for the best neighbor
    if (!adjacency_list.empty()) {
      for (int neighbor_idx : adjacency_list[i]) {
        if (particles[neighbor_idx].best_value < lbest_val) {

          // If the neighbor is better, we update the local best
          lbest_val = particles[neighbor_idx].best_value;
          lbest_pos = particles[neighbor_idx].best_position;
        }
      }
    } else {
      // If the adjacency list is empty, we use the global best
      lbest_pos = gbest_pos;
    }

    // We update the velocity and the position of every particle using the PSO equations
    for (size_t d = 0; d < active_dims.size(); ++d) {
      double r1 = dist01(gen);
      double r2 = dist01(gen);
      p.velocity[d] = w * p.velocity[d] +
                      c1 * r1 * (p.best_position[d] - p.position[d]) +
                      c2 * r2 * (lbest_pos[d] - p.position[d]);

      p.position[d] += p.velocity[d];

      // We apply the bounds-checking to the particle
      if (p.position[d] < bounds_lower) {
        p.position[d] = bounds_lower; 
        p.velocity[d] = -0.5 * p.velocity[d];
      } else if (p.position[d] > bounds_upper) {
        p.position[d] = bounds_upper;
        p.velocity[d] = -0.5 * p.velocity[d];
      }
    }
  }
}

void SubSwarm::evaluate_and_update(ContextVector &ctx, const TestFunction &f) {

  for (auto &p : particles) {
    // Evaluates for every particles its current fitness
    p.current_value = ctx.evaluate_particle(f, p.position, active_dims);

    // If the current position is better than the best position, we update the best position
    if (p.current_value < p.best_value) {
      p.best_value = p.current_value;
      p.best_position = p.position;
    }

    // If the best position is better than the global best position, we update the global best position
    if (p.best_value < gbest_val) {
      gbest_val = p.best_value;
      gbest_pos = p.best_position;
    }
  }
}

// Getter methods
const std::vector<double> &SubSwarm::get_gbest_pos() const { return gbest_pos; }
double SubSwarm::get_gbest_val() const { return gbest_val; }
const std::vector<int> &SubSwarm::get_active_dims() const { return active_dims; }
const std::vector<CPSOParticle> &SubSwarm::get_particles() const { return particles; }
