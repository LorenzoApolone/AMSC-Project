#ifndef METHODS_HPP
#define METHODS_HPP
#include "interfaces.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <vector>

OutputObject pso_serial(const TestFunction& f, int d,  StoppingCriteriaManager& stop, int n_points, unsigned int seed = 12345);
OutputObject pso_mpi(const TestFunction &f, int d,  StoppingCriteriaManager &stop, int n_points, unsigned int seed = 12345);

#endif