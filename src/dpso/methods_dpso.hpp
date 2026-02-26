#ifndef METHODS_HPP
#define METHODS_HPP

#include "interfaces.hpp"

OutputObject pso_mpi(const TestFunction& f, 
                     unsigned int dim, 
                     const StopCriterion& stop, 
                     unsigned int n_points_per_rank);

#endif