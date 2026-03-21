/**
 * @file SubSwarmTopologyFactory.hpp
 * @brief Declares topology helpers used to build CPSO sub-swarm neighborhoods.
 */
#pragma once

#include "../topology/NetworkType.hpp"
#include <random>
#include <vector>

/**
 * @struct SubSwarmTopologyConfig
 * @brief Describes the topology parameters used for one CPSO sub-swarm.
 */
struct SubSwarmTopologyConfig {
  NetworkType type = NetworkType::SCALE_FREE;
  double small_world_rewire_probability = 0.3;
  int scale_free_links_per_node = 2;
  double random_edge_probability = 0.5;

  /**
   * @brief Builds a configuration with defaults for the selected topology type.
   */
  static SubSwarmTopologyConfig from_type(NetworkType type) {
    SubSwarmTopologyConfig config;
    config.type = type;
    return config;
  }
};

/**
 * @brief Validates a topology configuration against the current sub-swarm size.
 */
void validate_subswarm_topology_config(const SubSwarmTopologyConfig &config,
                                       int particles_per_swarm);

/**
 * @brief Builds the adjacency list for a sub-swarm using an explicit configuration.
 */
std::vector<std::vector<int>>
build_subswarm_topology(const SubSwarmTopologyConfig &config,
                        int particles_per_swarm, std::mt19937 &gen);

/**
 * @brief Builds the adjacency list for a sub-swarm from a topology type only.
 */
std::vector<std::vector<int>> build_subswarm_topology(NetworkType type,
                                                      int particles_per_swarm,
                                                      std::mt19937 &gen);
