/**
 * @file SubSwarmTopologyFactory.cpp
 * @brief Implements the topology builders used by CPSO sub-swarms.
 */
#include "SubSwarmTopologyFactory.hpp"
#include "../topology/create_network.hpp"
#include <algorithm>
#include <stdexcept>

void validate_subswarm_topology_config(const SubSwarmTopologyConfig &config,
                                       int particles_per_swarm) {
  if (particles_per_swarm <= 0) {
    throw std::invalid_argument("particles per swarm must be > 0");
  }

  switch (config.type) {
  case NetworkType::SMALL_WORLD:
    if (config.small_world_rewire_probability < 0.0 ||
        config.small_world_rewire_probability > 1.0) {
      throw std::invalid_argument(
          "small-world rewiring probability must be in [0, 1]");
    }
    break;
  case NetworkType::SCALE_FREE:
    if (config.scale_free_links_per_node <= 0) {
      throw std::invalid_argument("scale-free links per node must be > 0");
    }
    break;
  case NetworkType::RANDOM:
    if (config.random_edge_probability < 0.0 ||
        config.random_edge_probability > 1.0) {
      throw std::invalid_argument(
          "random graph edge probability must be in [0, 1]");
    }
    break;
  default:
    throw std::invalid_argument("unsupported NetworkType for CPSO subswarm");
  }
}

std::vector<std::vector<int>>
build_subswarm_topology(const SubSwarmTopologyConfig &config,
                        int particles_per_swarm, std::mt19937 &gen) {
  validate_subswarm_topology_config(config, particles_per_swarm);

  if (particles_per_swarm == 1) {
    return std::vector<std::vector<int>>(1);
  }

  // Delegate the concrete graph construction to the shared topology builders.
  std::vector<std::vector<int>> adjacency_list;
  switch (config.type) {
  case NetworkType::SMALL_WORLD:
    create_network(particles_per_swarm, config.small_world_rewire_probability,
                   adjacency_list, gen);
    break;
  case NetworkType::SCALE_FREE:
    create_scale_free_network(
        particles_per_swarm,
        std::min(config.scale_free_links_per_node, particles_per_swarm - 1),
        adjacency_list, gen);
    break;
  case NetworkType::RANDOM:
    create_random_network(particles_per_swarm, config.random_edge_probability,
                          adjacency_list, gen);
    break;
  default:
    throw std::invalid_argument("unsupported NetworkType for CPSO subswarm");
  }

  return adjacency_list;
}

std::vector<std::vector<int>> build_subswarm_topology(NetworkType type,
                                                      int particles_per_swarm,
                                                      std::mt19937 &gen) {
  return build_subswarm_topology(SubSwarmTopologyConfig::from_type(type),
                                 particles_per_swarm, gen);
}
