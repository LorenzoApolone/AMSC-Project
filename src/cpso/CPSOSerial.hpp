#pragma once

#include "CPSOBase.hpp"

class CPSOSerial : public CPSOBase {
public:
  CPSOSerial(int k_subswarms, int num_particles_per_swarm, NetworkType topology,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const std::vector<NetworkType> &topologies,
             int shuffle_freq = 50, int stagnation_patience = 50,
             double w_start = 0.9, double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

protected:
  OutputObject run_optimization_loop(const TestFunction &f,
                                     StoppingCriteriaManager &stop_manager,
                                     std::vector<SubSwarm>& swarms, std::vector<std::mt19937>& gens,
                                     ContextVector& context, std::mt19937& global_gen) override;
};
