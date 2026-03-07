#include "cpso/CPSOSerial.hpp"
#include "functions.cpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <dim> <k_subswarms> <particles_per_swarm> [max_iters]\n";
    return 1;
  }

  unsigned int dim = std::stoi(argv[1]);
  int k_subswarms = std::stoi(argv[2]);
  int particles_per_swarm = std::stoi(argv[3]);
  int max_iters = argc > 4 ? std::stoi(argv[4]) : 1000;

  std::vector<std::string> test_names = {"Sphere",
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
                                         "ModifiedXinSheYang5",
                                         "Levy",
                                         "Michalewicz",
                                         "Bohachevsky",
                                         "Powell",
                                         "DixonPrice",
                                         "StyblinskiTang"};

  auto get_function = [](const std::string &name,
                         unsigned int d) -> std::unique_ptr<TestFunction> {
    if (name == "Sphere")
      return std::make_unique<Sphere>(d);
    if (name == "Ellipsoid")
      return std::make_unique<Ellipsoid>(d);
    if (name == "SumOfDiffPowers")
      return std::make_unique<SumOfDiffPowers>(d);
    if (name == "DropWave")
      return std::make_unique<DropWave>(d);
    if (name == "Weierstrass")
      return std::make_unique<Weierstrass>(d);
    if (name == "Alpine1")
      return std::make_unique<Alpine1>(d);
    if (name == "Ackley")
      return std::make_unique<Ackley>(d);
    if (name == "Griewank")
      return std::make_unique<Griewank>(d);
    if (name == "Rastrigin")
      return std::make_unique<Rastrigin>(d);
    if (name == "HappyCat")
      return std::make_unique<HappyCat>(d);
    if (name == "HGBat")
      return std::make_unique<HGBat>(d);
    if (name == "Rosenbrock")
      return std::make_unique<Rosenbrock>(d);
    if (name == "HighCondElliptic")
      return std::make_unique<HighCondElliptic>(d);
    if (name == "Discus")
      return std::make_unique<Discus>(d);
    if (name == "BentCigar")
      return std::make_unique<BentCigar>(d);
    if (name == "PermdbFunc")
      return std::make_unique<PermdbFunc>(d);
    if (name == "Schafferf7Func")
      return std::make_unique<Schafferf7Func>(d);
    if (name == "ExpSchafferF6")
      return std::make_unique<ExpSchafferF6>(d);
    if (name == "RotatedHyper")
      return std::make_unique<RotatedHyper>(d);
    if (name == "Schwefel")
      return std::make_unique<Schwefel>(d);
    if (name == "SumOfDifferentPowers2")
      return std::make_unique<SumOfDifferentPowers2>(d);
    if (name == "XinSheYang1")
      return std::make_unique<XinSheYang1>(d);
    if (name == "Schwefel221")
      return std::make_unique<Schwefel221>(d);
    if (name == "Schwefel222")
      return std::make_unique<Schwefel222>(d);
    if (name == "Salomon")
      return std::make_unique<Salomon>(d);
    if (name == "ModifiedRidge")
      return std::make_unique<ModifiedRidge>(d);
    if (name == "Zakharov")
      return std::make_unique<Zakharov>(d);
    if (name == "ModifiedXinSheYang3")
      return std::make_unique<ModifiedXinSheYang3>(d);
    if (name == "ModifiedXinSheYang5")
      return std::make_unique<ModifiedXinSheYang5>(d);
    if (name == "Levy")
      return std::make_unique<Levy>(d);
    if (name == "Michalewicz")
      return std::make_unique<Michalewicz>(d);
    if (name == "Bohachevsky")
      return std::make_unique<Bohachevsky>(d);
    if (name == "Powell")
      return std::make_unique<Powell>(d);
    if (name == "DixonPrice")
      return std::make_unique<DixonPrice>(d);
    if (name == "StyblinskiTang")
      return std::make_unique<StyblinskiTang>(d);
    return nullptr;
  };

  std::cout << "Running tests with dim=" << dim << ", k=" << k_subswarms
            << ", particles/swarm=" << particles_per_swarm
            << ", max_iters=" << max_iters << "\n\n";

  for (const auto &name : test_names) {
    std::cout << "=== Testing " << name << " ===\n";

    // Test CPSO-S
    auto f1 = get_function(name, dim);
    int scaled_stagnation = std::max(100, max_iters / 4);
    StoppingCriteriaManager stop1(max_iters, scaled_stagnation, 1e-8);
    CPSOSerial cpso_s(k_subswarms, particles_per_swarm,
                      NetworkType::SCALE_FREE);
    OutputObject out_s = cpso_s.optimize(*f1, stop1);

    double true_min = f1->value(f1->get_true_solution());
    double fitness_error = std::abs(out_s.f_val - true_min);
    double tol = 1e-4; // threshold for convergence

    std::string conv_status;
    if (fitness_error < tol) {
        conv_status = "True Convergence";
    } else if (stop1.get_current_iters() < stop1.get_max_iters()) {
        conv_status = "False Convergence (Stagnation)";
    } else {
        conv_status = "No Convergence (Max Iterations)";
    }

    std::cout << "[CPSO-S] Best Fitness: " << out_s.f_val
              << " | Error: " << fitness_error
              << " | Iters: " << stop1.get_current_iters()
              << " | Time: " << out_s.execution_time << "s\n"
              << "             Status: " << conv_status << "\n\n";
  }

  return 0;
}
