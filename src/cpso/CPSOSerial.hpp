#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include "../topology/create_network.hpp"
#include "ContextVector.hpp"
#include "SubSwarm.hpp"
#include <random>
#include <vector>

class CPSOSerial {
private:
  int num_subswarms;
  int particles_per_swarm;
  double w_max, w_min, c1, c2;
  std::vector<NetworkType> subswarm_topologies;

public:
  CPSOSerial(int k_subswarms, int num_particles_per_swarm, NetworkType topology,
             double w_start = 0.9, double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

  CPSOSerial(int k_subswarms, int num_particles_per_swarm,
             const std::vector<NetworkType> &topologies, double w_start = 0.9,
             double w_end = 0.4, double coeff1 = 1.49618,
             double coeff2 = 1.49618);

  OutputObject optimize(const TestFunction &f,
                        StoppingCriteriaManager &stop_manager);
};
