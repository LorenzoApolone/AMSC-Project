/**
 * @file CPSOBenchmarkUtils.hpp
 * @brief Helper utilities shared by the CPSO benchmark drivers.
 */
#pragma once

#include "../interfaces/functions.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @namespace cpso_benchmark
 * @brief Collects registry, convergence and reporting helpers for CPSO tests.
 */
namespace cpso_benchmark {

// Maps function names to builders.
using FunctionFactory = std::unordered_map<
    std::string, std::function<std::unique_ptr<TestFunction>(unsigned int)>>;

// Signature used by the static benchmark registry.
using FunctionBuilder = std::unique_ptr<TestFunction> (*)(unsigned int);

/**
 * @struct FunctionRegistryEntry
 * @brief One entry in the static CPSO benchmark registry.
 */
struct FunctionRegistryEntry {
  // Name of the Benchmark function.
  const char *name;

  // Factory function used to instantiate the benchmark.
  FunctionBuilder builder;
};

/**
 * @brief Instantiates one benchmark function of the requested type.
 * @tparam FunctionType Concrete benchmark function type.
 * @param dim Dimension used to build the function.
 * @return A newly allocated benchmark function of type `FunctionType`.
 */
template <typename FunctionType>
std::unique_ptr<TestFunction> make_function(unsigned int dim) {
  return std::make_unique<FunctionType>(dim);
}

/**
 * @brief Returns the static registry of benchmark functions exercised by CPSO tests.
 * @return An immutable registry shared by the serial and parallel CPSO drivers.
 */
inline const std::vector<FunctionRegistryEntry> &get_function_registry() {
  static const std::vector<FunctionRegistryEntry> registry = {
      {"Sphere", &make_function<Sphere>},
      {"Ellipsoid", &make_function<Ellipsoid>},
      {"SumOfDiffPowers", &make_function<SumOfDiffPowers>},
      {"QuinticFunction", &make_function<QuinticFunction>},
      {"DropWave", &make_function<DropWave>},
      {"Weierstrass", &make_function<Weierstrass>},
      {"Alpine1", &make_function<Alpine1>},
      {"Ackley", &make_function<Ackley>},
      {"Griewank", &make_function<Griewank>},
      {"Rastrigin", &make_function<Rastrigin>},
      {"HappyCat", &make_function<HappyCat>},
      {"HGBat", &make_function<HGBat>},
      {"Rosenbrock", &make_function<Rosenbrock>},
      {"HighCondElliptic", &make_function<HighCondElliptic>},
      {"Discus", &make_function<Discus>},
      {"BentCigar", &make_function<BentCigar>},
      {"Schafferf7Func", &make_function<Schafferf7Func>},
      {"ExpSchafferF6", &make_function<ExpSchafferF6>},
      {"RotatedHyper", &make_function<RotatedHyper>},
      {"Schwefel", &make_function<Schwefel>},
      {"SumOfDifferentPowers2", &make_function<SumOfDifferentPowers2>},
      {"XinSheYang1", &make_function<XinSheYang1>},
      {"Schwefel221", &make_function<Schwefel221>},
      {"Schwefel222", &make_function<Schwefel222>},
      {"Salomon", &make_function<Salomon>},
      {"ModifiedRidge", &make_function<ModifiedRidge>},
      {"Zakharov", &make_function<Zakharov>},
      {"ModifiedXinSheYang3", &make_function<ModifiedXinSheYang3>},
      {"ModifiedXinSheYang5", &make_function<ModifiedXinSheYang5>},
      {"Levy", &make_function<Levy>},
      {"Michalewicz", &make_function<Michalewicz>},
      {"Bohachevsky", &make_function<Bohachevsky>},
      {"Powell", &make_function<Powell>},
      {"DixonPrice", &make_function<DixonPrice>},
      {"StyblinskiTang", &make_function<StyblinskiTang>},
      {"Step", &make_function<Step>},
      {"Qing", &make_function<Qing>},
      {"Trid", &make_function<Trid>},
      {"Shubert", &make_function<Shubert>},
      {"Alpine2", &make_function<Alpine2>},
      {"Eggholder", &make_function<Eggholder>},
      {"Easom", &make_function<Easom>},
      {"Brown", &make_function<Brown>},
      {"Csendes", &make_function<Csendes>},
      {"Vincent", &make_function<Vincent>},
  };
  return registry;
}

/**
 * @brief Returns the ordered list of benchmark names from the static registry.
 * @return A vector of benchmark names that mirrors the static registry order.
 */
inline const std::vector<std::string> &get_test_names() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> collected;
    collected.reserve(get_function_registry().size());
    for (const auto &entry : get_function_registry()) {
      collected.emplace_back(entry.name);
    }
    return collected;
  }();
  return names;
}

/**
 * @brief Builds a hash map factory from the static function registry.
 * @return A table from function name to builder.
 */
inline FunctionFactory build_factory() {
  FunctionFactory factory;

  // Convert the ordered registry into a name-based lookup table for the drivers.
  for (const auto &entry : get_function_registry()) {
    factory[entry.name] = entry.builder;
  }
  return factory;
}

/**
 * @brief Converts a function typology enum to label.
 * @param typology Function typology to convert.
 * @return Enum of the typology for benchmark logs.
 */
inline std::string typology_to_string(FunctionTypology typology) {
  switch (typology) {
  case FunctionTypology::UNIMODAL:
    return "Unimodal";
  case FunctionTypology::MULTIMODAL:
    return "Multimodal";
  case FunctionTypology::SEPARABLE:
    return "Separable";
  case FunctionTypology::NON_SEPARABLE:
    return "Non-separable";
  case FunctionTypology::DIFFERENTIABLE:
    return "Differentiable";
  case FunctionTypology::NON_DIFFERENTIABLE:
    return "Non-differentiable";
  case FunctionTypology::FLAT:
    return "Flat";
  case FunctionTypology::COUPLED:
    return "Coupled";
  default:
    return "Unknown";
  }
}

/**
 * @struct ConvergenceResult
 * @brief Captures the convergence outcome of CPSO benchmark runs.
 */
struct ConvergenceResult {
  bool converged = false;
  bool failed = false;
  std::string status;
  double fitness_error = std::numeric_limits<double>::infinity();
};

/**
 * @brief Infers the benchmark-layer label for a run that did not reach the target tolerance.
 * @param stop_manager Stopping manager used during the run.
 * @return The non-convergence label that best explains why the run stopped.
 */
inline std::string
infer_non_convergence_status(const StoppingCriteriaManager &stop_manager) {
  if (stop_manager.get_current_iters() >= stop_manager.get_max_iters()) {
    return "No Convergence (Max Iterations)";
  }

  if (stop_manager.get_current_stagnation_iters() >=
      stop_manager.get_max_stagnation_iters()) {
    return "False Convergence (Stagnation)";
  }

  return "False Convergence (Low Diversity)";
}

/**
 * @brief Evaluates whether a benchmark run reached the requested tolerance.
 * @param function Benchmark function associated with the run.
 * @param result Generic solver output produced by the benchmark driver.
 * @param stop_manager Stopping manager used during the run.
 * @param tolerance Absolute fitness error threshold required for true convergence.
 * @return A structured convergence result that distinguishes success, failure and non-convergence.
 */
inline ConvergenceResult
evaluate_convergence(const TestFunction &function, const OutputObject &result,
                     const StoppingCriteriaManager &stop_manager,
                     double tolerance = 1e-4) {
  const double best_fitness = result.get_best_fitness();

  if (!std::isfinite(best_fitness)) {
    return {false, true, "Failure (Non-finite Best Fitness)",
            std::numeric_limits<double>::infinity()};
  }

  const double true_min = function.value(function.get_true_solution());
  if (!std::isfinite(true_min)) {
    return {false, true, "Failure (Non-finite Reference Fitness)",
            std::numeric_limits<double>::infinity()};
  }

  const double fitness_error = std::abs(best_fitness - true_min);
  if (!std::isfinite(fitness_error)) {
    return {false, true, "Failure (Non-finite Fitness Error)",
            std::numeric_limits<double>::infinity()};
  }

  if (fitness_error < tolerance) {
    return {true, false, "True Convergence", fitness_error};
  }

  return {false, false, infer_non_convergence_status(stop_manager),
          fitness_error};
}

/**
 * @struct BenchmarkSummary
 * @brief Aggregates convergence statistics by function typology.
 */
struct BenchmarkSummary {
  // Number of tested functions per typology.
  std::map<FunctionTypology, int> total_by_typology;

  // Number of converged functions per typology.
  std::map<FunctionTypology, int> converged_by_typology;

  // Number of failed functions per typology.
  std::map<FunctionTypology, int> failed_by_typology;

  // Descriptions for runs that did not converge but remained valid.
  std::vector<std::string> non_converged_functions;

  //Descriptions for technically failed runs.
  std::vector<std::string> failed_functions;

  /**
   * @brief Records the outcome of one benchmark function.
   * @param function_name Name of the tested benchmark function.
   * @param function Benchmark function instance associated with the run.
   * @param result Convergence outcome computed for that run.
   */
  void record(const std::string &function_name, const TestFunction &function,
              const ConvergenceResult &result) {
    // A function can belong to multiple typologies, each one is counted independently.
    for (auto typology : function.get_typologies()) {
      total_by_typology[typology]++;
      if (result.converged) {
        converged_by_typology[typology]++;
      }
      if (result.failed) {
        failed_by_typology[typology]++;
      }
    }

    if (result.converged) {
      return;
    }

    std::ostringstream line;
    line << "  - " << function_name;

    if (!function.get_typologies().empty()) {
      line << " [";
      for (size_t i = 0; i < function.get_typologies().size(); ++i) {
        if (i > 0) {
          line << ", ";
        }
        line << typology_to_string(function.get_typologies()[i]);
      }
      line << "]";
    }

    line << ": " << result.status << ", error=" << result.fitness_error;
    if (result.failed) {
      failed_functions.push_back(line.str());
      return;
    }

    non_converged_functions.push_back(line.str());
  }

  /**
   * @brief Prints the collected summary in a human-readable form.
   * @param os Output stream that receives the formatted summary.
   * @param label Section title printed before the summary body.
   */
  void print(std::ostream &os, const std::string &label) const {
    // Print aggregate counts first, then the detailed lists of individual functions.
    os << label << "\n";
    os << "Converged by typology:\n";

    if (total_by_typology.empty()) {
      os << "  - No typology metadata available.\n";
    } else {
      for (const auto &[typology, total_count] : total_by_typology) {
        int converged_count = 0;
        auto converged_it = converged_by_typology.find(typology);
        if (converged_it != converged_by_typology.end()) {
          converged_count = converged_it->second;
        }
        os << "  - " << typology_to_string(typology) << ": "
           << converged_count << "/" << total_count << "\n";
      }
    }

    os << "Failed by typology:\n";
    if (total_by_typology.empty()) {
      os << "  - No typology metadata available.\n";
    } else {
      for (const auto &[typology, total_count] : total_by_typology) {
        int failed_count = 0;
        auto failed_it = failed_by_typology.find(typology);
        if (failed_it != failed_by_typology.end()) {
          failed_count = failed_it->second;
        }
        os << "  - " << typology_to_string(typology) << ": "
           << failed_count << "/" << total_count << "\n";
      }
    }

    os << "Non-converged functions:\n";
    if (non_converged_functions.empty()) {
      os << "  - None\n";
    } else {
      for (const auto &function_line : non_converged_functions) {
        os << function_line << "\n";
      }
    }

    os << "Failed functions:\n";
    if (failed_functions.empty()) {
      os << "  - None\n";
    } else {
      for (const auto &function_line : failed_functions) {
        os << function_line << "\n";
      }
    }

    os << "\n";
  }
};

} // namespace cpso_benchmark
