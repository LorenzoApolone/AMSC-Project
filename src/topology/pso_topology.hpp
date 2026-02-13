#pragma once

#include "../interfaces.hpp"
#include <vector>




OutputObject pso_normal(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list, bool &converged);

OutputObject pso_small_timerv(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list, bool &converged, double &t_allgatherv_tot);
