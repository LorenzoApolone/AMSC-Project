#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include "stop_criterion.hpp"

class test_function {
public:
    test_function(const std::function<double(const std::vector<double>&)> &f_,
                  const std::vector<double> &true_solution_,
                  std::vector<std::array<double, 2>> &domain_,
                  const std::string name_);

    double operator()(const std::vector<double> &x) {
        return f(x);
    }

    // It computes the error w.r.t. the true solution.
    double error(const std::vector<double> &x) const;

    // It returns the domain as a vector of dimension d with
    // [(lower bound_1, upper bound_1), ..., (lower_bound_d, upper_bound_d)]
    const std::vector<std::array<double, 2>> &get_domain() const;

    const std::vector<double> &get_true_solution() const {
        return true_solution;
    }

    const std::string &to_string() const {
        return name;
    }

private:
    const std::function<double(const std::vector<double>&)> f;
    const StopCriterion stop_criterion;
    const std::vector<double> true_solution;
    const std::vector<std::array<double, 2>> domain;
    const std::string name;
};

class output_object {
public:
    // Criterion_idx indicates which kind of stopping criterion
    // the chosen method used. 
    // - 0: maximum delta_x
    // - 1: maximum number of iterations
    // - 2: patience 
    // TODO:
    // what is patience?
    output_object(const std::string function_name_,
                  int n_points_,
                  unsigned int n_cores_,
                  const StopCriterion &stop_criterion_,
                  unsigned int criterion_idx_,
                  const std::vector<double> &found_x_,
                  double found_value_,
                  const std::vector<double> &true_solution_,
                  const std::vector<double> &conv_history_,
                  double execution_time_);

    void terminal_info();
    void output_to_file();

private:
    const std::string function_name;
    const int n_points;
    const unsigned int n_cores;
    const StopCriterion stop_criterion;
    const unsigned int criterion_idx;
    const std::vector<double> found_x;
    const double found_value;
    const std::vector<double> true_solution;
    const std::vector<double> conv_history;
    const double execution_time;
};

// Note: it is not crazy to believe that the number
// of requested points may differ at the end of the program.
// Thus, this input is not a const ref.
output_object method(const test_function &f,
                     const StopCriterion &stop_criterion,
                     unsigned int dim,
                     unsigned int n_points,
                     unsigned int n_cores);

#endif
