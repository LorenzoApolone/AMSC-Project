#ifndef METHODS_HPP
#define METHODS_HPP

#include "interfaces.hpp"

OutputObject dpso(const TestFunction& f, 
                     unsigned int dim, 
                     unsigned int n_points_total, 
                     int max_iter);

#endif