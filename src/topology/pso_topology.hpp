#pragma once

#include "../interfaces.hpp"
#include <vector>


OutputObject pso_small(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list,
                       bool &converged);

OutputObject pso_small_debugger(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list, bool &converged);

OutputObject pso_small_timer(const TestFunction &f,
                       int d,
                       const StopCriterion &stop,
                       int n_points,
                       const std::vector<std::vector<int>> &adjacency_list, bool &converged, double &t_allgatherv_tot_out);