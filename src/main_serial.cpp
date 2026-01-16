/**
 * @file main_serial.cpp
 * @brief Run Particle Swarm Optimization (serial) over a suite of benchmark
 * functions.
 *
 * This program parses problem/configuration parameters from the command line,
 * constructs a set of standard continuous optimization test functions, and runs
 * a serial PSO solver on each function in turn. For each run, it prints summary
 * information to the terminal (and can optionally write results to a file).
 *
 * @usage
 *   ./main_serial <dim> <n_points> <max_iter> <delta_x>
 *
 * @param dim       Problem dimensionality (unsigned int > 0).
 * @param n_points  Swarm size / number of particles (unsigned int > 0).
 * @param max_iter  Maximum PSO iterations (unsigned int > 0).
 * @param delta_x   Stopping tolerance, typically a step/position tolerance
 * (double > 0).
 *
 * @requirements
 * - The following headers must provide the corresponding
 * interfaces/definitions:
 *     - methods.hpp: declares pso_serial(...) returning OutputObject.
 *     - interfaces.hpp: declares TestFunction, StopCriterion, OutputObject.
 *     - functions.cpp: defines the concrete TestFunction subclasses used below.
 * - C++17 (or later) for <memory>, std::make_unique, etc.
 */

#include "functions.cpp"
#include "interfaces.hpp"
#include "methods.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0]
              << " <dim> <n_points> <max_iter> <delta_x>\n";
    return 1;
  }

  unsigned int dim = std::atoi(argv[1]);
  unsigned int n_points = std::atoi(argv[2]);
  unsigned int max_iter = std::atoi(argv[3]);
  double delta_x = std::atof(argv[4]);

  StopCriterion stop_criterion(max_iter, delta_x);

  // Factory Definition
  std::unordered_map<std::string,
                     std::function<std::unique_ptr<TestFunction>(unsigned int)>>
      factory;

  factory["Sphere"] = [](unsigned int dim) {
    return std::make_unique<Sphere>(dim);
  };
  factory["Ellipsoid"] = [](unsigned int dim) {
    return std::make_unique<Ellipsoid>(dim);
  };
  factory["SumOfDiffPowers"] = [](unsigned int dim) {
    return std::make_unique<SumOfDiffPowers>(dim);
  };
  factory["DropWave"] = [](unsigned int dim) {
    return std::make_unique<DropWave>(dim);
  };
  factory["Weierstrass"] = [](unsigned int dim) {
    return std::make_unique<Weierstrass>(dim);
  };
  factory["Alpine1"] = [](unsigned int dim) {
    return std::make_unique<Alpine1>(dim);
  };
  factory["Ackley"] = [](unsigned int dim) {
    return std::make_unique<Ackley>(dim);
  };
  factory["Griewank"] = [](unsigned int dim) {
    return std::make_unique<Griewank>(dim);
  };
  factory["Rastrigin"] = [](unsigned int dim) {
    return std::make_unique<Rastrigin>(dim);
  };
  factory["HappyCat"] = [](unsigned int dim) {
    return std::make_unique<HappyCat>(dim);
  };
  factory["HGBat"] = [](unsigned int dim) {
    return std::make_unique<HGBat>(dim);
  };
  factory["Rosenbrock"] = [](unsigned int dim) {
    return std::make_unique<Rosenbrock>(dim);
  };
  factory["HighCondElliptic"] = [](unsigned int dim) {
    return std::make_unique<HighCondElliptic>(dim);
  };
  factory["Discus"] = [](unsigned int dim) {
    return std::make_unique<Discus>(dim);
  };
  factory["BentCigar"] = [](unsigned int dim) {
    return std::make_unique<BentCigar>(dim);
  };
  factory["PermdbFunc"] = [](unsigned int dim) {
    return std::make_unique<PermdbFunc>(dim);
  };
  factory["Schafferf7Func"] = [](unsigned int dim) {
    return std::make_unique<Schafferf7Func>(dim);
  };
  factory["ExpSchafferF6"] = [](unsigned int dim) {
    return std::make_unique<ExpSchafferF6>(dim);
  };
  factory["RotatedHyper"] = [](unsigned int dim) {
    return std::make_unique<RotatedHyper>(dim);
  };
  factory["Schwefel"] = [](unsigned int dim) {
    return std::make_unique<Schwefel>(dim);
  };
  factory["SumOfDifferentPowers2"] = [](unsigned int dim) {
    return std::make_unique<SumOfDifferentPowers2>(dim);
  };
  factory["XinSheYang1"] = [](unsigned int dim) {
    return std::make_unique<XinSheYang1>(dim);
  };
  factory["Schwefel221"] = [](unsigned int dim) {
    return std::make_unique<Schwefel221>(dim);
  };
  factory["Schwefel222"] = [](unsigned int dim) {
    return std::make_unique<Schwefel222>(dim);
  };
  factory["Salomon"] = [](unsigned int dim) {
    return std::make_unique<Salomon>(dim);
  };
  factory["ModifiedRidge"] = [](unsigned int dim) {
    return std::make_unique<ModifiedRidge>(dim);
  };
  factory["Zakharov"] = [](unsigned int dim) {
    return std::make_unique<Zakharov>(dim);
  };
  factory["ModifiedXinSheYang3"] = [](unsigned int dim) {
    return std::make_unique<ModifiedXinSheYang3>(dim);
  };
  factory["ModifiedXinSheYang5"] = [](unsigned int dim) {
    return std::make_unique<ModifiedXinSheYang5>(dim);
  };

  // List of functions to execute (maintaining the previous order)
  std::vector<std::string> function_names = {"Sphere",
                                             "Ellipsoid",
                                             "SumOfDiffPowers",
                                             "DropWave",
                                             "Weierstrass",
                                             "Alpine1",
                                             "Ackley",
                                             "Griewank",
                                             "Rastrigin",
                                             "HappyCat",
                                             "HGBat",
                                             "Rosenbrock",
                                             "HighCondElliptic",
                                             "Discus",
                                             "BentCigar",
                                             "PermdbFunc",
                                             "Schafferf7Func",
                                             "ExpSchafferF6",
                                             "RotatedHyper",
                                             "Schwefel",
                                             "SumOfDifferentPowers2",
                                             "XinSheYang1",
                                             "Schwefel221",
                                             "Schwefel222",
                                             "Salomon",
                                             "ModifiedRidge",
                                             "Zakharov",
                                             "ModifiedXinSheYang3",
                                             "ModifiedXinSheYang5"};

  // Run the solver
  for (const auto &name : function_names) {
    auto f_ptr = factory[name](dim);
    OutputObject result = pso_serial(*f_ptr, dim, stop_criterion, n_points);
    result.terminal_info();
    // result.output_to_file();
  }

  std::cout << "\nAll tests completed successfully." << std::endl;
  return 0;
}
