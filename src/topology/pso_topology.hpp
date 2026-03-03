#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include <vector>



OutputObject pso_topology(const TestFunction &f,
                       int d,
                        StoppingCriteriaManager &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list,  double &t_allgatherv_tot);
