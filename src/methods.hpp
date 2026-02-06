#ifndef METHODS_HPP
#define METHODS_HPP
#include "interfaces.hpp"
#include <vector>

OutputObject pso_serial(const TestFunction& f, int d, const StopCriterion& stop, int n_points);
OutputObject pso_mpi(const TestFunction &f, int d, const StopCriterion &stop, int n_points, bool &converged);

#endif