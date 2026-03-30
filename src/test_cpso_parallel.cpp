#include "cpso/CPSOBenchmarkUtils.hpp"
#include "cpso/CPSOParallel.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 4) {
    if (rank == 0) {
      std::cerr << "Usage: " << argv[0]
                << " <dim> <k_subswarms> <particles_per_swarm> [max_iters] "
                   "[shuffle_freq] [stagnation_patience] [seed]\n";
    }
    MPI_Finalize();
    return 1;
  }

  unsigned int dim = std::stoi(argv[1]);
  int k_subswarms = std::stoi(argv[2]);
  int particles_per_swarm = std::stoi(argv[3]);
  int max_iters = argc > 4 ? std::stoi(argv[4]) : 1000;
  int shuffle_freq = argc > 5 ? std::stoi(argv[5]) : 50;
  int stagnation_patience = argc > 6 ? std::stoi(argv[6]) : 50;
  unsigned int seed =
      argc > 7 ? static_cast<unsigned int>(std::stoul(argv[7]))
               : CPSOBase::DEFAULT_SEED;
  auto env_flag_enabled = [](const char *name) {
    const char *raw_value = std::getenv(name);
    if (raw_value == nullptr || *raw_value == '\0') {
      return false;
    }

    std::string value(raw_value);
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) {
                     return static_cast<char>(std::tolower(ch));
                   });
    return value == "1" || value == "true" || value == "yes" ||
           value == "on";
  };
  const bool disable_greedy_merge =
      env_flag_enabled("CPSO_MPI_DISABLE_GREEDY_MERGE") ||
      env_flag_enabled("CPSO_PARALLEL_DISABLE_GREEDY_MERGE");

  auto factory = cpso_benchmark::build_factory();
  cpso_benchmark::BenchmarkSummary summary;

  if (rank == 0) {
    std::cout << "Running MPI tests with dim=" << dim << ", k=" << k_subswarms
              << ", particles/swarm=" << particles_per_swarm
              << ", max_iters=" << max_iters
              << ", shuffle_freq=" << shuffle_freq
              << ", stagnation_patience=" << stagnation_patience
              << ", seed=" << seed
              << ", greedy_merge_fallback="
              << (disable_greedy_merge ? "disabled" : "enabled")
              << "\n\n";
  }

  for (const auto &name : cpso_benchmark::get_test_names()) {
    if (rank == 0) {
      std::cout << "=== Testing " << name << " ===\n";
    }

    auto factory_it = factory.find(name);
    if (factory_it == factory.end()) {
      if (rank == 0) {
        std::cerr << "Factory not found for function " << name << "\n";
      }
      continue;
    }

    auto f2 = factory_it->second(dim);

    int scaled_stagnation = std::max(100, max_iters / 4);
    StoppingCriteriaManager stop2(max_iters, scaled_stagnation, 1e-8);
    CPSOParallel cpso_p(k_subswarms, particles_per_swarm,
                        NetworkType::SCALE_FREE, shuffle_freq,
                        stagnation_patience, 0.9, 0.4, 1.49618, 1.49618, seed);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();
    CpsoRunArtifacts artifacts = cpso_p.optimize_raw(*f2, stop2);
    OutputObject out_p = build_cpso_output(*f2, artifacts, stop2);
    out_p.n_points = static_cast<unsigned int>(particles_per_swarm * k_subswarms);
    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();
    out_p.execution_time = t_end - t_start;

    if (rank == 0) {
      auto result = cpso_benchmark::evaluate_convergence(*f2, out_p, stop2);
      result.status = describe_cpso_status(artifacts, result.status);
      summary.record(name, *f2, result);

      std::cout << "[CPSO-P MPI] Best Fitness: " << out_p.f_val
                << " | Error: " << result.fitness_error
                << " | Iters: " << stop2.get_current_iters()
                << " | Time: " << (t_end - t_start) << "s\n"
                << "             Status: " << result.status << "\n"
                << "             Comm: total=" << artifacts.comm_total_s
                << "s | compute=" << artifacts.compute_total_s
                << "s | allgather=" << artifacts.comm_allgather_s
                << "s | bcast=" << artifacts.comm_bcast_s
                << "s | allreduce=" << artifacts.comm_allreduce_s
                << "s | barrier=" << artifacts.comm_barrier_s
                << "s | wait=" << artifacts.wait_total_s << "s\n\n";
    }
  }

  if (rank == 0) {
    summary.print(std::cout, "=== Final CPSO-P Summary ===");
  }

  MPI_Finalize();
  return 0;
}
