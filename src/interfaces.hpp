/**
 * @file interfaces.hpp
 * @brief Contains core PSO interfaces: TestFunction, StopCriterion, and OutputObject.
 */

#ifndef INTERFACES_HPP
#define INTERFACES_HPP

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/**
 * @class TestFunction
 * @brief Abstract base class for optimization problems.
 */
class TestFunction
{
public:
  /**
   * @brief Constructor for the TestFunction.
   * @param dim_ Dimensionality of the problem
   * @param function_name_ The name of the function
   * @param domain_ A pair defining the search space boundaries [lower, upper]
   * @param true_solution_ Global optimum vector
   */
  TestFunction(unsigned int dim_,
               string function_name_,
               pair<double, double> domain_,
               const vector<double> true_solution_) : dim(dim_),
                                                      function_name(function_name_),
                                                      domain(domain_),
                                                      true_solution(true_solution_) {};

  unsigned int dim;                       ///< Dimensionality
  const string function_name;             ///< Function name
  const pair<double, double> domain;      ///< Search domain
  const vector<double> true_solution;     ///< Global optimum

  /// @return Function name
  const string get_name() const { return function_name; }

  /// @return Search domain pair [min, max]
  const pair<double, double> &get_domain() const { return domain; };

  /// @brief Calculates fitness f(x), this is a pure virtual function
  virtual double value(const vector<double> &x) const = 0;

  /// @return Known true solution vector x*
  const vector<double> &get_true_solution() const { return true_solution; };

  /**
   * @brief Calculates normalized RMSE error between x and true solution
   */
  double error(const vector<double> &x) const
  {
    auto truth = get_true_solution();
    const pair<double, double> bounds = get_domain();
    double range = bounds.second - bounds.first;
    if (range == 0) range = 1.0;

    double D = (double)x.size();
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
      double normalized_diff = (x[i] - truth[i]) / range;
      sum += pow(normalized_diff, 2);
    }
    return sqrt(sum / D);
  }
};

/**
 * @class StopCriterion
 * @brief Logic for stopping the algorithm based on iterations or tolerance
 */
class StopCriterion
{
  int max_iterations;
  double tol;

public:
  StopCriterion(int max_iter, double tolerance)
      : max_iterations(max_iter), tol(tolerance) {}

  /// @brief Checks if stopping condition is met
  bool should_stop(int current_iter, double current_error_metric) const
  {
    if (current_iter >= max_iterations) return true;
    if (abs(current_error_metric) < tol) return true;
    return false;
  }

  int get_max_iter() const { return max_iterations; }
  double get_tolerance() const { return tol; }
};

/**
 * @class OutputObject
 * @brief Stores final results and metrics
 */
class OutputObject
{
  void output_to_file(string filename)
  {
    cout << "Writing results to " << filename << "..." << endl;
    cout << "Final Fitness: " << f_val << " | Time: " << execution_time << "s" << endl;
  }

public:
  string function_name;             ///< Problem name
  int d;                            ///< Dimensions
  unsigned int n_points;            ///< Particle count
  vector<double> x_best;            ///< Best position found
  vector<double> x_star;            ///< True solution
  double f_val;                     ///< Best fitness value
  vector<double> conv_history;      ///< Convergence history
  int n_cores;                      ///< MPI cores used
  double execution_time;            ///< Execution time (in seconds)
  int iterations;                   ///< Total iterations
  const StopCriterion &stopcriterion;

  OutputObject(string function_name_,
               int d_,
               unsigned int n_points_,
               vector<double> x_best_,
               vector<double> x_star_,
               double f_val_,
               vector<double> conv_history_,
               int n_cores_,
               double execution_time_,
               int iterations_,
               const StopCriterion &stopcriterion_)
      : function_name(function_name_),
        d(d_),
        n_points(n_points_),
        x_best(x_best_),
        x_star(x_star_),
        f_val(f_val_),
        conv_history(conv_history_),
        n_cores(n_cores_),
        execution_time(execution_time_),
        iterations(iterations_),
        stopcriterion(stopcriterion_) {};

  void terminal_info();
  void output_to_file();
};

#endif