/**
 * @file CPSOSerial.hpp
 * @brief Declares the serial Cooperative PSO solver.
 */
#pragma once

#include "CPSOBase.hpp"

/**
 * @class CPSOSerial
 * @brief Executes the CPSO loop sequentially on a single process.
 */
class CPSOSerial : public CPSOBase {
public:
  /** @brief Constructs a serial solver with a shared topology type. */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm, NetworkType topology,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4,
             double coeff1 = 1.49618, double coeff2 = 1.49618,
             unsigned int seed = DEFAULT_SEED);

  /** @brief Constructs a serial solver with a shared topology configuration. */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const SubSwarmTopologyConfig &topology_config,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4,
             double coeff1 = 1.49618, double coeff2 = 1.49618,
             unsigned int seed = DEFAULT_SEED);

  /** @brief Constructs a serial solver with one topology type per sub-swarm. */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const std::vector<NetworkType> &topologies,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4,
             double coeff1 = 1.49618, double coeff2 = 1.49618,
             unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a serial solver with one explicit topology configuration per sub-swarm.
   */
  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const std::vector<SubSwarmTopologyConfig> &topologies,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4,
             double coeff1 = 1.49618, double coeff2 = 1.49618,
             unsigned int seed = DEFAULT_SEED);

protected:
  /** @copydoc CPSOBase::run_optimization_loop */
  CpsoRunArtifacts run_optimization_loop(const TestFunction &f,
                                         StoppingCriteriaManager &stop_manager,
                                         std::vector<SubSwarm> &swarms,
                                         std::vector<std::mt19937> &gens,
                                         ContextVector &context,
                                         std::mt19937 &global_gen) override;
};
