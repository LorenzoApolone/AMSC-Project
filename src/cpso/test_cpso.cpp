#include "CPSOBenchmarkUtils.hpp"
#include "CPSOSerial.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <dim> <k_subswarms> <particles_per_swarm> [max_iters] [seed]\n";
    return 1;
  }

  unsigned int dim = std::stoi(argv[1]);
  int k_subswarms = std::stoi(argv[2]);
  int particles_per_swarm = std::stoi(argv[3]);
  int max_iters = argc > 4 ? std::stoi(argv[4]) : 1000;
  unsigned int seed =
      argc > 5 ? static_cast<unsigned int>(std::stoul(argv[5]))
               : CPSOBase::DEFAULT_SEED;

  auto factory = cpso_benchmark::build_factory();
  cpso_benchmark::BenchmarkSummary summary;

  std::cout << "Running tests with dim=" << dim << ", k=" << k_subswarms
            << ", particles/swarm=" << particles_per_swarm
            << ", max_iters=" << max_iters << ", seed=" << seed
            << "\n\n";

  for (const auto &name : cpso_benchmark::get_test_names()) {
    std::cout << "=== Testing " << name << " ===\n";

    auto factory_it = factory.find(name);
    if (factory_it == factory.end()) {
      std::cerr << "Factory not found for function " << name << "\n";
      continue;
    }

    auto f1 = factory_it->second(dim);
    int scaled_stagnation = std::max(100, max_iters / 4);
    StoppingCriteriaManager stop1(max_iters, scaled_stagnation, 1e-8);
    CPSOSerial cpso_s(k_subswarms, particles_per_swarm,
                      NetworkType::SCALE_FREE, 50, 50, 0.9, 0.4,
                      1.49618, 1.49618, seed);
    auto t_start = std::chrono::high_resolution_clock::now();
    CpsoRunArtifacts artifacts = cpso_s.optimize_raw(*f1, stop1);
    auto t_end = std::chrono::high_resolution_clock::now();
    OutputObject out_s = build_cpso_output(*f1, artifacts, stop1);
    out_s.n_points = static_cast<unsigned int>(particles_per_swarm * k_subswarms);
    out_s.execution_time =
        std::chrono::duration<double>(t_end - t_start).count();

    auto result = cpso_benchmark::evaluate_convergence(*f1, out_s, stop1);
    result.status = describe_cpso_status(artifacts, result.status);
    summary.record(name, *f1, result);

    std::cout << "[CPSO-S] Best Fitness: " << out_s.f_val
              << " | Error: " << result.fitness_error
              << " | Iters: " << stop1.get_current_iters()
              << " | Time: " << out_s.execution_time << "s\n"
              << "             Status: " << result.status << "\n\n";
  }

  summary.print(std::cout, "=== Final CPSO-S Summary ===");

  return 0;
}
