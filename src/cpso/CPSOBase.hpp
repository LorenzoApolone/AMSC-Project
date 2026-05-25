/**
 * @file CPSOBase.hpp
 * @brief Defines the shared base class used by the serial and MPI CPSO solvers.
 */
#pragma once

#include "../interfaces/interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "CPSORunArtifacts.hpp"
#include "ContextVector.hpp"
#include "SubSwarm.hpp"
#include "SubSwarmTopologyFactory.hpp"
#include <random>
#include <vector>

/**
 * @class CPSOBase
 * @brief Provides the common CPSO configuration.
 *
 * The base class owns the solver hyperparameters, builds the initial
 * sub-swarms and context vector, and delegates the iteration loop to the
 * concrete serial or parallel implementation.
 */
class CPSOBase {
private:
  int num_subswarms;
  int particles_per_swarm;
  int shuffle_freq;
  int stagnation_patience;
  double w_max, w_min;
  double c1, c2;
  std::vector<SubSwarmTopologyConfig> subswarm_topologies;
  unsigned int master_seed;

protected:
  int get_num_subswarms() const { return num_subswarms; }
  int get_particles_per_swarm() const { return particles_per_swarm; }
  int get_shuffle_freq() const { return shuffle_freq; }
  int get_stagnation_patience() const { return stagnation_patience; }
  double get_w_max() const { return w_max; }
  double get_w_min() const { return w_min; }
  double get_c1() const { return c1; }
  double get_c2() const { return c2; }
  const std::vector<SubSwarmTopologyConfig> &get_subswarm_topologies() const {
    return subswarm_topologies;
  }
  unsigned int get_master_seed() const { return master_seed; }

public:
  /** @brief Default deterministic seed used when no seed is provided. */
  static constexpr unsigned int DEFAULT_SEED = 5489u;

  /**
   * @brief Constructs a solver where all sub-swarms use the same topology family with its default parameters.
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           NetworkType topology, int shuffle_freq,
           int stagnation_patience, double w_start,
           double w_end, double coeff1, double coeff2,
           unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a solver where all sub-swarms share the same fully specified topology configuration.
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           const SubSwarmTopologyConfig &topology_config, int shuffle_freq,
           int stagnation_patience, double w_start, double w_end,
           double coeff1, double coeff2, unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a solver with one topology family per sub-swarm, each using the default parameters of that family.
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           const std::vector<NetworkType> &topologies,
           int shuffle_freq, int stagnation_patience,
           double w_start, double w_end, double coeff1, double coeff2,
           unsigned int seed = DEFAULT_SEED);

  /**
   * @brief Constructs a solver with one fully specified topology configuration for each sub-swarm.
   */
  CPSOBase(int k_subswarms, int num_particles_per_swarm,
           const std::vector<SubSwarmTopologyConfig> &topologies,
           int shuffle_freq, int stagnation_patience, double w_start,
           double w_end, double coeff1, double coeff2,
           unsigned int seed = DEFAULT_SEED);

  virtual ~CPSOBase() = default;

  /**
   * @brief Executes CPSO and returns solver artifacts.
   * @param f Benchmark function to optimize.
   * @param stop_manager Stopping manager.
   * @return Raw CPSO artifacts containing the best solution and runtime metadata.
   */
  CpsoRunArtifacts optimize_raw(const TestFunction &f,
                                StoppingCriteriaManager &stop_manager);

  /**
   * @brief Executes CPSO and converts the result into the project-wide output type.
   * @param f Benchmark function to optimize.
   * @param stop_manager Stopping manager.
   * @return Output object used by the surrounding benchmark pipeline.
   */
  OutputObject optimize(const TestFunction &f,
                        StoppingCriteriaManager &stop_manager);

protected:
  /**
   * @brief Runs the solver-specific optimization loop.
   * @param f Benchmark function to optimize.
   * @param stop_manager Stopping manager.
   * @param swarms Prepared sub-swarms for the run.
   * @param gens Per-subswarm random generators.
   * @param context Shared CPSO context vector.
   * @param global_gen Generator used for global operations such as reshuffling.
   * @return Raw CPSO artifacts for the completed run.
   */
  virtual CpsoRunArtifacts
  run_optimization_loop(const TestFunction &f,
                        StoppingCriteriaManager &stop_manager,
                        std::vector<SubSwarm> &swarms,
                        std::vector<std::mt19937> &gens, ContextVector &context,
                        std::mt19937 &global_gen) = 0;
};
