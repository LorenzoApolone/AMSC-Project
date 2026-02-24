#include "cpso/CPSOSerial.hpp"
#include "functions.cpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <dim> <k_subswarms> <particles_per_swarm> [max_fevals]\n";
    return 1;
  }

  unsigned int dim = std::stoi(argv[1]);
  int k_subswarms = std::stoi(argv[2]);
  int particles_per_swarm = std::stoi(argv[3]);
  int max_fevals = argc > 4 ? std::stoi(argv[4]) : 20000;

  std::vector<std::string> test_names = {"ModifiedXinSheYang3",
                                         "ModifiedXinSheYang5", "Michalewicz",
                                         "Powell", "DixonPrice"};

  auto get_function = [](const std::string &name,
                         unsigned int d) -> std::unique_ptr<TestFunction> {
    if (name == "ModifiedXinSheYang3")
      return std::make_unique<ModifiedXinSheYang3>(d);
    if (name == "ModifiedXinSheYang5")
      return std::make_unique<ModifiedXinSheYang5>(d);
    if (name == "Michalewicz")
      return std::make_unique<Michalewicz>(d);
    if (name == "Powell")
      return std::make_unique<Powell>(d);
    if (name == "DixonPrice")
      return std::make_unique<DixonPrice>(d);
    return nullptr;
  };

  std::cout << "Running tests with dim=" << dim << ", k=" << k_subswarms
            << ", particles/swarm=" << particles_per_swarm
            << ", max_fevals=" << max_fevals << "\n\n";

  for (const auto &name : test_names) {
    std::cout << "=== Testing " << name << " ===\n";

    // Test CPSO-S
    auto f1 = get_function(name, dim);
    StoppingCriteriaManager stop1(max_fevals, 100, 1e-8);
    CPSOSerial cpso_s(k_subswarms, particles_per_swarm);
    OutputObject out_s = cpso_s.optimize(*f1, stop1);

    std::cout << "[CPSO-S] Best Fitness: " << out_s.f_val
              << " | FEvals: " << stop1.get_current_fevals()
              << " | Iters: " << out_s.iterations
              << " | Time: " << out_s.execution_time << "s\n";
  }

  return 0;
}
