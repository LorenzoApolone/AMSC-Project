/**
 * @file CPSORunArtifacts.hpp
 * @brief Defines the native CPSO result type and helper conversion utilities.
 */
#pragma once

#include "../interfaces.hpp"
#include "../interfaces/StoppingCriteriaManager.hpp"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

/**
 * @enum CpsoStopReason
 * @brief Describes why a CPSO run terminated.
 */
enum class CpsoStopReason {
  // The stopping cause could not be established.
  UNKNOWN,

  // The maximum iteration budget was exhausted.
  MAX_ITERATIONS,

  // The run stopped because progress stalled.
  STAGNATION,

  // The swarm diversity collapsed below the configured limit.
  LOW_DIVERSITY,

  // The run aborted because a numeric guard failed.
  NUMERIC_FAILURE,
};

/**
 * @struct CpsoRunArtifacts
 * @brief Native result object returned by the CPSO solvers.
 */
struct CpsoRunArtifacts {
  // Best context vector found during the run.
  std::vector<double> best_position;

  // Best fitness value associated with best_position.
  double best_fitness = std::numeric_limits<double>::infinity();

  // Best-fitness history collected during the optimization loop.
  std::vector<double> best_fitness_history;

  // Effective seed used by the run.
  unsigned int seed_used = 0u;

  // Number of processes or cores involved in the run.
  int cores = 1;

  // Number of iterations completed by the solver.
  int iterations = 0;

  // Aggregated active MPI communication time.
  double comm_total_s = 0.0;

  // Time spent in MPI_Allgather calls.
  double comm_allgather_s = 0.0;

  // Time spent in MPI_Bcast calls.
  double comm_bcast_s = 0.0;

  // Time spent in MPI_Allreduc calls.
  double comm_allreduce_s = 0.0;

  // Time spent in MPI_Barrier calls.
  double comm_barrier_s = 0.0;

  // Aggregated synchronization wait time.
  double wait_total_s = 0.0;

  // Estimated wall time spent outside MPI communication.
  double compute_total_s = 0.0;

  // Native CPSO stop reason reported by the solver.
  CpsoStopReason stop_reason = CpsoStopReason::UNKNOWN;

  // Optional detail message used when the run ends in numeric failure.
  std::string failure_message;
};

/**
 * @brief Maps stopping-manager flags to a CPSO stop reason.
 * @param stop_for_max_iters true when the run exhausted the iteration budget.
 * @param stop_for_low_diversity true when the swarm diversity fell below the configured limit.
 * @param stop_for_stagnation true when stagnation triggered the stop condition.
 * @return The CPSO stop reason that best matches the provided flags.
 */
inline CpsoStopReason infer_cpso_stop_reason(bool stop_for_max_iters,
                                             bool stop_for_low_diversity,
                                             bool stop_for_stagnation) {
  // Preserve the same priority order used by the CPSO stopping logic.
  if (stop_for_max_iters) {
    return CpsoStopReason::MAX_ITERATIONS;
  }
  if (stop_for_low_diversity) {
    return CpsoStopReason::LOW_DIVERSITY;
  }
  if (stop_for_stagnation) {
    return CpsoStopReason::STAGNATION;
  }
  return CpsoStopReason::UNKNOWN;
}

/**
 * @brief Converts a best-fitness history into an absolute error history.
 * @param f Benchmark function used to compute the reference optimum.
 * @param best_fitness_history Best-fitness history produced by the solver.
 * @return The absolute error history with respect to the known optimum of f.
 */
inline std::vector<double>
build_cpso_error_history(const TestFunction &f,
                         const std::vector<double> &best_fitness_history) {
  std::vector<double> error_history;
  error_history.reserve(best_fitness_history.size());

  const double true_minimum = f.value(f.get_true_solution());

  if (!std::isfinite(true_minimum)) {
    error_history.assign(best_fitness_history.size(),
                         std::numeric_limits<double>::infinity());
    return error_history;
  }

  for (double best_fitness : best_fitness_history) {
    if (!std::isfinite(best_fitness)) {
      error_history.push_back(std::numeric_limits<double>::infinity());
      continue;
    }
    error_history.push_back(std::abs(best_fitness - true_minimum));
  }

  return error_history;
}

/**
 * @brief Builds the status string exposed by the benchmark layer.
 * @param artifacts Native CPSO artifacts produced by the solver.
 * @param fallback_status Status string inferred by the generic benchmark layer.
 * @return A label string consistent with the native CPSO stop reason.
 */
inline std::string describe_cpso_status(
    const CpsoRunArtifacts &artifacts, const std::string &fallback_status) {
  // Keep the generic benchmark status unless CPSO explicitly reported a numeric failure.
  if (artifacts.stop_reason != CpsoStopReason::NUMERIC_FAILURE) {
    return fallback_status;
  }

  if (artifacts.failure_message.empty()) {
    return "Failure (Numeric Failure)";
  }

  return "Failure (Numeric Failure: " + artifacts.failure_message + ")";
}

/**
 * @brief Converts native CPSO artifacts to the project-wide benchmark output object.
 * @param f Benchmark function associated with the run.
 * @param artifacts Native CPSO artifacts produced by the solver.
 * @param stop_manager Stopping manager used during the run.
 * @return An OutputObject compatible with the rest of the benchmark framework.
 */
inline OutputObject build_cpso_output(const TestFunction &f,
                                      const CpsoRunArtifacts &artifacts,
                                      StoppingCriteriaManager &stop_manager) {
  // Repackage the native CPSO result in the generic format expected by the benchmark drivers.
  return OutputObject(
      f.get_name(), f.dim, 0u, artifacts.best_position, f.get_true_solution(),
      artifacts.best_fitness,
      build_cpso_error_history(f, artifacts.best_fitness_history),
      artifacts.cores, 0.0, artifacts.iterations, stop_manager);
}