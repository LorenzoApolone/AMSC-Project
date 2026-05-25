/**
 * @file SubSwarm.cpp
 * @brief Implementation of the SubSwarm class methods
 */
#include "SubSwarm.hpp"
#include "NumericValidation.hpp"
#include <algorithm>
#include <stdexcept>

namespace {

void validate_dimension_set(const std::vector<int> &dims) {
  if (dims.empty()) {
    throw std::invalid_argument("active dimensions must not be empty");
  }

  // Sort a copy so uniqueness can be checked without changing the order.
  std::vector<int> sorted_dims = dims;
  std::sort(sorted_dims.begin(), sorted_dims.end());
  auto duplicate_it = std::adjacent_find(sorted_dims.begin(), sorted_dims.end());
  if (duplicate_it != sorted_dims.end()) {
    throw std::invalid_argument("active dimensions must be unique");
  }

  if (sorted_dims.front() < 0) {
    throw std::invalid_argument("active dimensions must be non-negative");
  }
}

void validate_adjacency_list(const std::vector<std::vector<int>> &adj_list,
                             int num_particles) {
  // The neighborhood graph must describe exactly the local particles.
  if (!adj_list.empty() &&
      adj_list.size() != static_cast<size_t>(num_particles)) {
    throw std::invalid_argument(
        "adjacency list size must match the number of particles");
  }

  for (const auto &neighbors : adj_list) {
    for (int neighbor_idx : neighbors) {
      if (neighbor_idx < 0 || neighbor_idx >= num_particles) {
        throw std::out_of_range(
            "adjacency list contains an invalid particle index");
      }
    }
  }
}

} // namespace

SubSwarm::SubSwarm(int num_particles, const std::vector<int> &active_dimensions,
                   double lower_bound, double upper_bound,
                   const std::vector<std::vector<int>> &adj_list)
    : num_particles(num_particles), active_dims(active_dimensions),
      gbest_val(std::numeric_limits<double>::infinity()),
      gbest_particle_idx(-1), ignore_gbest_this_iter(false),
      adjacency_list(adj_list), bounds_lower(lower_bound),
      bounds_upper(upper_bound) {
      
  // Safety Checks      
  if (num_particles <= 0) {
    throw std::invalid_argument("subswarm must contain at least one particle");
  }
  ensure_finite_value(lower_bound, "subswarm lower bound");
  ensure_finite_value(upper_bound, "subswarm upper bound");
  if (bounds_lower > bounds_upper) {
    throw std::invalid_argument("lower bound must be <= upper bound");
  }
  validate_dimension_set(active_dimensions);
  validate_adjacency_list(adj_list, num_particles);

  int dim = active_dims.size();
  gbest_pos.resize(dim, 0.0);

  positions.resize(num_particles * dim);
  velocities.resize(num_particles * dim);
  best_positions.resize(num_particles * dim);
  current_values.resize(num_particles, std::numeric_limits<double>::infinity());
  best_values.resize(num_particles, std::numeric_limits<double>::infinity());
}

void SubSwarm::initialize(std::mt19937 &gen, ContextVector &ctx, const TestFunction &f) {
  // Initialize positions in the domain with small initial velocities around zero.
  std::uniform_real_distribution<double> dist_pos(bounds_lower, bounds_upper);
  double range = bounds_upper - bounds_lower;

  std::uniform_real_distribution<double> dist_vel(-0.1 * range, 0.1 * range);

  const std::vector<double> &context_vec = ctx.get_full_vector();
  ensure_finite_vector(context_vec, "context vector before subswarm initialization");

  // Initialize positions, velocities and personal-best memories.
  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {
    std::vector<double> temp_pos(dim);
    for (int d = 0; d < dim; ++d) {
      int idx = i * dim + d;
      if (i == 0) {
        // Use the current context as a stable anchor particle.
        positions[idx] = context_vec[active_dims[d]];
        velocities[idx] = 0.0;
      } else {
        positions[idx] = dist_pos(gen);
        velocities[idx] = dist_vel(gen);
      }
      best_positions[idx] = positions[idx];
      temp_pos[d] = positions[idx];
    }

    ensure_finite_vector(temp_pos, "initialized particle position");

    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);
    best_values[i] = current_values[i];

    // Track the best local candidate seen during initialization.
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

    ensure_finite_vector(temp_best_pos, "particle personal-best position");
    ensure_finite_vector(temp_pos, "particle current position");

    best_values[i] = ctx.evaluate_particle(f, temp_best_pos, active_dims);

    // The current position may now dominate the stored memory under the new context.
    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);

    // Promote the current position if the reshaped context makes it better.
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
  validate_dimension_set(new_dims);

  const std::vector<int> old_active_dims = active_dims;
  const std::vector<double> old_positions = positions;
  const std::vector<double> old_velocities = velocities;
  const std::vector<double> old_best_positions = best_positions;

  // Switch the subswarm to the new cooperative block.
  active_dims = new_dims;

  // Reshuffling invalidates the previous subswarm best.
  gbest_val = std::numeric_limits<double>::infinity();
  gbest_particle_idx = -1;
  gbest_pos.assign(active_dims.size(), 0.0);

  const std::vector<double> &context_vec = ctx.get_full_vector();
  ensure_finite_vector(context_vec, "context vector before subswarm initialization");
  int dim = active_dims.size();

    // Map each new local coordinate to the matching coordinate in the previous block, when it exists.
  std::vector<int> old_local_index(dim, -1);
  for (int new_local = 0; new_local < dim; ++new_local) {
    for (int old_local = 0; old_local < dim; ++old_local) {
      if (old_active_dims[old_local] == active_dims[new_local]) {
        old_local_index[new_local] = old_local;
        break;
      }
    }
  }

  // Placeholders keep only the reshuffled coordinates needed to mirror the shared context.
  if (!is_owned) {
    for (int i = 0; i < num_particles; ++i) {
      current_values[i] = std::numeric_limits<double>::infinity();
      best_values[i] = std::numeric_limits<double>::infinity();
      for (int d = 0; d < dim; ++d) {
        int idx = i * dim + d;
        positions[idx] = context_vec[active_dims[d]];
        best_positions[idx] = positions[idx];
        velocities[idx] = 0.0;
        ensure_finite_value(positions[idx], "non-owned reshuffle position");
      }
    }
    return;
  }

  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_pos(-0.25 * range, 0.25 * range);
  std::uniform_real_distribution<double> dist_vel(-0.05 * range, 0.05 * range);

  // Preserve coordinates and memories only for dimensions that remain in the same subswarm.
  for (int i = 0; i < num_particles; ++i) {
    current_values[i] = std::numeric_limits<double>::infinity();
    best_values[i] = std::numeric_limits<double>::infinity();

    for (int d = 0; d < dim; ++d) {
      int idx = i * dim + d;
      int old_local = old_local_index[d];

      if (i == 0) {
        positions[idx] = context_vec[active_dims[d]];
        velocities[idx] = 0.0;
        best_positions[idx] = positions[idx];
        continue;
      }

      if (old_local >= 0) {
        int old_idx = i * dim + old_local;
        positions[idx] = old_positions[old_idx];
        velocities[idx] = old_velocities[old_idx];
        best_positions[idx] = old_best_positions[old_idx];
      } else {
        positions[idx] = context_vec[active_dims[d]] + dist_pos(gen);
        velocities[idx] = dist_vel(gen);
      }

      positions[idx] = std::clamp(positions[idx], bounds_lower, bounds_upper);
      if (old_local < 0) {
        best_positions[idx] = positions[idx];
      }

      ensure_finite_value(positions[idx], "reshuffled particle position");
      ensure_finite_value(velocities[idx], "reshuffled particle velocity");
      ensure_finite_value(best_positions[idx], "reshuffled particle personal best");
    }
  }
}


void SubSwarm::inject_velocities(std::mt19937 &gen, bool hard_reset) {
  double range = bounds_upper - bounds_lower;
  std::uniform_real_distribution<double> dist_vel(-0.15 * range, 0.15 * range);

  // Kick every non-best particle away from its current trajectory.
  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {
    // Compare by particle index instead of floating-point equality.
    bool is_gbest = (i == gbest_particle_idx);

    // Optionally clear local memory so the particle can explore a new region.
    if (!is_gbest) {
      for (int d = 0; d < dim; ++d) {
        int idx = i * dim + d;

        // Use a wider perturbation to produce a visible escape step.
        velocities[idx] = dist_vel(gen) * 2.0; 
        ensure_finite_value(velocities[idx], "injected particle velocity");
        
        if (hard_reset) {
          // Reset the personal memory so cognitive component doesn't pull it back instantly
          best_positions[idx] = positions[idx]; 
        }
      }
      if (hard_reset) {
        best_values[i] = current_values[i]; // Keep the personal-best fitness aligned with the reset memory
      }
    }
  }
}

void SubSwarm::reset_gbest_attraction() {
  ignore_gbest_this_iter = true;
}

void SubSwarm::update_velocities_and_positions(double w, double c1, double c2,
                                               std::mt19937 &gen, double progress_ratio) {
  ensure_finite_value(progress_ratio, "progress ratio");
  ensure_finite_value(w, "inertia weight");
  ensure_finite_value(c1, "cognitive coefficient");
  ensure_finite_value(c2, "social coefficient");
  std::uniform_real_distribution<double> dist01(0.0, 1.0);


  int dim = active_dims.size();
  for (int i = 0; i < num_particles; ++i) {

    std::vector<double> lbest_pos(dim);
    int lbest_particle_idx = i;
    for (int d = 0; d < dim; ++d) {
        lbest_pos[d] = best_positions[i * dim + d];
    }
    double lbest_val = best_values[i];

    // Use the neighborhood graph when present, otherwise fall back to the subswarm best.
    if (!adjacency_list.empty()) {
      for (int neighbor_idx : adjacency_list[i]) {
        if (best_values[neighbor_idx] < lbest_val) {

          // Keep the best neighborhood guide available to this particle.
          lbest_val = best_values[neighbor_idx];
          lbest_particle_idx = neighbor_idx;
          for (int d = 0; d < dim; ++d) {
              lbest_pos[d] = best_positions[neighbor_idx * dim + d];
          }
        }
      }
    } else {
      // Without a neighborhood graph, the subswarm best acts as the local guide.
      lbest_pos = gbest_pos;
      lbest_val = gbest_val;
      lbest_particle_idx = gbest_particle_idx;
    }

    double v_max = (bounds_upper - bounds_lower) * (0.2 - 0.15 * progress_ratio);
    ensure_finite_value(v_max, "velocity clamp range");
    if (v_max < 0.0) {
      throw std::runtime_error("velocity clamp range must be >= 0");
    }

    // Apply the standard PSO update with the selected local guide.
    for (int d = 0; d < dim; ++d) {
      double r1 = dist01(gen);
      double r2 = dist01(gen);
      
      int idx = i * dim + d;

      double c1_act = (best_values[i] == std::numeric_limits<double>::infinity()) ? 0.0 : c1;
      // Temporarily remove social attraction when the caller requested one exploration-only step.
      const bool ignore_selected_global_best =
          ignore_gbest_this_iter && gbest_particle_idx >= 0 &&
          lbest_particle_idx == gbest_particle_idx;
      double c2_act =
          (lbest_val == std::numeric_limits<double>::infinity() ||
           ignore_selected_global_best)
              ? 0.0
              : c2;

      velocities[idx] = w * velocities[idx] +
                      c1_act * r1 * (best_positions[idx] - positions[idx]) +
                      c2_act * r2 * (lbest_pos[d] - positions[idx]);

      if (v_max == 0.0) {
        velocities[idx] = 0.0;
      } else if (velocities[idx] > v_max) {
        velocities[idx] = v_max;
      } else if (velocities[idx] < -v_max) {
        velocities[idx] = -v_max;
      }

      positions[idx] += velocities[idx];

      // Reflect particles back into the domain with a bounce.
      if (positions[idx] < bounds_lower) {
        positions[idx] = bounds_lower; 
        velocities[idx] = -0.5 * velocities[idx];
      } else if (positions[idx] > bounds_upper) {
        positions[idx] = bounds_upper;
        velocities[idx] = -0.5 * velocities[idx];
      }

      ensure_finite_value(velocities[idx], "updated particle velocity");
      ensure_finite_value(positions[idx], "updated particle position");
    }
  }

  // Reset the temporary flag after the velocity update
  ignore_gbest_this_iter = false;
}

void SubSwarm::evaluate_and_update(ContextVector &ctx, const TestFunction &f) {
  int dim = active_dims.size();
  std::vector<double> temp_pos(dim);
  std::vector<double> temp_best_pos(dim);

  for (int i = 0; i < num_particles; ++i) {
    for (int d = 0; d < dim; ++d) {
        temp_pos[d] = positions[i * dim + d];
        temp_best_pos[d] = best_positions[i * dim + d];
    }
    
    ensure_finite_vector(temp_pos, "particle position before evaluation");

    current_values[i] = ctx.evaluate_particle(f, temp_pos, active_dims);

    // Refresh the personal best if the new position improves it.
    if (current_values[i] < best_values[i]) {
      best_values[i] = current_values[i];
      for (int d = 0; d < dim; ++d) {
          best_positions[i * dim + d] = positions[i * dim + d];
          temp_best_pos[d] = positions[i * dim + d];
      }
    }

    // Update the subswarm best from the refreshed personal memories.
    if (best_values[i] < gbest_val) {
      gbest_val = best_values[i];
      gbest_pos = temp_best_pos;
      gbest_particle_idx = i;
    }
  }
}

const std::vector<double> &SubSwarm::get_gbest_pos() const { return gbest_pos; }
double SubSwarm::get_gbest_val() const { return gbest_val; }
const std::vector<int> &SubSwarm::get_active_dims() const { return active_dims; }
int SubSwarm::get_num_particles() const { return num_particles; }
const std::vector<double> &SubSwarm::get_positions() const { return positions; }
