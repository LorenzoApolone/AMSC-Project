/**
 * @file CPSOParallel.hpp
 * @brief Declares the MPI-based Cooperative PSO solver.
 */
#pragma once

#include "CPSOBase.hpp"

/**
 * @class CPSOParallel
 * @brief Executes CPSO with MPI-aware subswarm ownership and synchronization.
 */
class CPSOParallel : public CPSOBase {
public:
  /**
   * @brief Constructs a parallel solver where all sub-swarms use the same topology family with its default parameters.
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               NetworkType topology, int shuffle_freq = 50,
               int stagnation_patience = 50, double w_start = 0.9,
               double w_end = 0.4, double coeff1 = 1.49618,
               double coeff2 = 1.49618,
               unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a parallel solver where all sub-swarms share the same fully specified topology configuration.
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               const SubSwarmTopologyConfig &topology_config,
               int shuffle_freq = 50, int stagnation_patience = 50,
               double w_start = 0.9, double w_end = 0.4,
               double coeff1 = 1.49618, double coeff2 = 1.49618,
               unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a parallel solver with one topology family per sub-swarm, each using the default parameters of that family.
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               const std::vector<NetworkType> &topologies,
               int shuffle_freq = 50, int stagnation_patience = 50,
               double w_start = 0.9, double w_end = 0.4,
               double coeff1 = 1.49618, double coeff2 = 1.49618,
               unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a parallel solver with one fully specified topology configuration for each sub-swarm.
   */
  CPSOParallel(int k_subswarms, int num_particles_per_swarm,
               const std::vector<SubSwarmTopologyConfig> &topologies,
               int shuffle_freq = 50, int stagnation_patience = 50,
               double w_start = 0.9, double w_end = 0.4,
               double coeff1 = 1.49618, double coeff2 = 1.49618,
               unsigned int seed = DEFAULT_SEED);

protected:
  /** @copydoc CPSOBase::use_distributed_swarm_ownership */
  bool use_distributed_swarm_ownership() const override { return true; }

  /** @copydoc CPSOBase::run_optimization_loop */
  CpsoRunArtifacts run_optimization_loop(const TestFunction &f,
                                         StoppingCriteriaManager &stop_manager,
                                         std::vector<SubSwarm> &swarms,
                                         std::vector<std::mt19937> &gens,
                                         ContextVector &context,
                                         std::mt19937 &global_gen) override;
};