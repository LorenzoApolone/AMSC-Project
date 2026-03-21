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
  // Topology family used to build the neighborhood graph.
  NetworkType type = NetworkType::SCALE_FREE;
  
  // Rewiring probability used by the small-world builder.
  double small_world_rewire_probability = 0.3;    

  // Number of links attached by each new node in a scale-free graph.
  int scale_free_links_per_node = 2;
  
  // Edge probability used by the random graph builder.
  double random_edge_probability = 0.5;           

  /**
   * @brief Builds a configuration with defaults for the selected topology type.
   * @param type Topology family to embed in the returned configuration.
   * @return A configuration initialized with the default parameters for that topology family.
   */
  static SubSwarmTopologyConfig from_type(NetworkType type) {
    SubSwarmTopologyConfig config;
    config.type = type;
    return config;
  }
};

/**
 * @brief Validates a topology configuration against the current sub-swarm size.
 * @param config Topology configuration to validate.
 * @param particles_per_swarm Number of particles in the target sub-swarm.
 */
void validate_subswarm_topology_config(const SubSwarmTopologyConfig &config,
                                       int particles_per_swarm);

/**
 * @brief Builds the adjacency list for a sub-swarm using an explicit configuration.
 * @param config Topology configuration to instantiate.
 * @param particles_per_swarm Number of particles in the target sub-swarm.
 * @param gen Random generator used by the stochastic topology builders.
 * @return The adjacency list that describes the neighborhood graph of the sub-swarm.
 */
std::vector<std::vector<int>>
build_subswarm_topology(const SubSwarmTopologyConfig &config,
                        int particles_per_swarm, std::mt19937 &gen);

/**
 * @brief Builds the adjacency list for a sub-swarm from a topology type only.
 * @param type Topology family to instantiate.
 * @param particles_per_swarm Number of particles in the target sub-swarm.
 * @param gen Random generator used by the stochastic topology builders.
 * @return The adjacency list that describes the neighborhood graph of the sub-swarm.
 */
std::vector<std::vector<int>> build_subswarm_topology(NetworkType type,
                                                      int particles_per_swarm,
                                                      std::mt19937 &gen);
