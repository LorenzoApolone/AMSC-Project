/**
 * @file interfaces.hpp
 * @brief Contains core PSO interfaces: TestFunction, StopCriterion, and OutputObject.
 */

#ifndef INTERFACES_HPP
#define INTERFACES_HPP
#include "interfaces/StoppingCriteriaManager.hpp"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

enum class FunctionTypology {
  UNIMODAL,
  MULTIMODAL,
  SEPARABLE,
  NON_SEPARABLE,
  DIFFERENTIABLE,
  NON_DIFFERENTIABLE,
  FLAT,
  COUPLED
};

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
               const vector<double> true_solution_,
               const vector<FunctionTypology> typologies_ = {}) : dim(dim_),
                                                      function_name(function_name_),
                                                      domain(domain_),
                                                      true_solution(true_solution_),
                                                      typologies(typologies_) {};

  unsigned int dim;                       ///< Dimensionality
  const string function_name;             ///< Function name
  const pair<double, double> domain;      ///< Search domain
  const vector<double> true_solution;     ///< Global optimum
  const vector<FunctionTypology> typologies; ///< Typologies of the function

  /// @return Function name
  const string get_name() const { return function_name; }

  /// @return Function typologies
  const vector<FunctionTypology>& get_typologies() const { return typologies; }

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
  
  // MPI Timing tracking for Benchmarks
  double comm_total_s = 0.0;
  double compute_total_s = 0.0;
  double comm_bcast_s = 0.0;
  double comm_allreduce_s = 0.0;
  double comm_allgather_s = 0.0; // In DPSO, we map Alltoallv to this for parity
  double comm_barrier_s = 0.0;
  double wait_total_s = 0.0;

  const StoppingCriteriaManager &stopcriterion;

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
               StoppingCriteriaManager &stopcriterion_)
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
void append_summary_csv_by_method(const std::string& method_name, int rep) const;
double get_best_fitness() const { return f_val; }
const std::vector<double>& get_best_position() const { return x_best; }
int get_iterations() const { return iterations; }
};


#endif