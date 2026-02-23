#include "SubSwarm.hpp"

SubSwarm::SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
                   double lower_bound, double upper_bound,
                   const std::vector<std::vector<int>> &adj_list)
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

void SubSwarm::initialize(std::mt19937 &gen, ContextVector &ctx,
                          const TestFunction &f) {
  std::uniform_real_distribution<double> dist_pos(bounds_lower, bounds_upper);
  // Velocity initialization at ±10% of the range
  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_vel(-0.1 * range, 0.1 * range);

  for (auto &p : particles) {
    for (size_t d = 0; d < active_dims.size(); ++d) {
      p.position[d] = dist_pos(gen);
      p.velocity[d] = dist_vel(gen);
      p.best_position[d] = p.position[d];
    }

    // Initial evaluation by combining with the Context Vector
    p.current_value = ctx.evaluate_particle(f, p.position, active_dims);
    p.best_value = p.current_value;

    // Updates the global best local to the sub-swarm
    if (p.current_value < gbest_val) {
      gbest_val = p.current_value;
      gbest_pos = p.position;
    }
  }
}

void SubSwarm::update_velocities_and_positions(double w, double c1, double c2,
                                               std::mt19937 &gen) {
  std::uniform_real_distribution<double> dist01(0.0, 1.0);

  for (size_t i = 0; i < particles.size(); ++i) {
    auto &p = particles[i];

    // 1. Determine the network's lBest (Local Best)
    std::vector<double> lbest_pos = p.best_position; // Itself as starting point
    double lbest_val = p.best_value;

    if (!adjacency_list.empty()) {
      for (int neighbor_idx : adjacency_list[i]) {
        if (particles[neighbor_idx].best_value < lbest_val) {
          lbest_val = particles[neighbor_idx].best_value;
          lbest_pos = particles[neighbor_idx].best_position;
        }
      }
    } else {
      // No list = Fully Connected (Global Topology)
      lbest_pos = gbest_pos;
    }

    // 2. Coordinate update
    for (size_t d = 0; d < active_dims.size(); ++d) {
      double r1 = dist01(gen);
      double r2 = dist01(gen);

      // Standard PSO equation, but influenced by lBest (Scale-Free network)
      p.velocity[d] = w * p.velocity[d] +
                      c1 * r1 * (p.best_position[d] - p.position[d]) +
                      c2 * r2 * (lbest_pos[d] - p.position[d]);

      p.position[d] += p.velocity[d];

      // Clamping to avoid exceeding permitted bounds
      if (p.position[d] < bounds_lower)
        p.position[d] = bounds_lower;
      if (p.position[d] > bounds_upper)
        p.position[d] = bounds_upper;
    }
  }
}

int SubSwarm::evaluate_and_update(ContextVector &ctx, const TestFunction &f) {
  int evaluations_done = 0;

  for (auto &p : particles) {
    // New evaluation
    p.current_value = ctx.evaluate_particle(f, p.position, active_dims);
    evaluations_done++;

    // Personal best update
    if (p.current_value < p.best_value) {
      p.best_value = p.current_value;
      p.best_position = p.position;
    }

    // Global best (local to the sub-swarm) update
    if (p.best_value < gbest_val) {
      gbest_val = p.best_value;
      gbest_pos = p.best_position;
    }
  }

  return evaluations_done;
}

const std::vector<double> &SubSwarm::get_gbest_pos() const { return gbest_pos; }
double SubSwarm::get_gbest_val() const { return gbest_val; }
const std::vector<int> &SubSwarm::get_active_dims() const {
  return active_dims;
}
const std::vector<CPSOParticle> &SubSwarm::get_particles() const {
  return particles;
}
