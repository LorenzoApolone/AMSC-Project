#include "cpso/CPSOParallel.hpp"
#include "functions.cpp"
#include "interfaces/StoppingCriteriaManager.hpp"
#include <iostream>
#include <memory>
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
                << " <dim> <k_subswarms> <particles_per_swarm> [max_fevals]\n";
    }
    MPI_Finalize();
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

  if (rank == 0) {
    std::cout << "Running MPI tests with dim=" << dim << ", k=" << k_subswarms
              << ", particles/swarm=" << particles_per_swarm
              << ", max_fevals=" << max_fevals << "\n\n";
  }

  for (const auto &name : test_names) {
    if (rank == 0) {
      std::cout << "=== Testing " << name << " ===\n";
    }

    auto f2 = get_function(name, dim);

    std::vector<NetworkType> topologies = {NetworkType::SCALE_FREE,
                                           NetworkType::SMALL_WORLD,
                                           NetworkType::FULLY_CONNECTED};

    StoppingCriteriaManager stop2(max_fevals, 100, 1e-8);
    CPSOParallel cpso_p(k_subswarms, particles_per_swarm, topologies);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();
    OutputObject out_p = cpso_p.optimize(*f2, stop2);
    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    if (rank == 0) {
      std::cout << "[CPSO-P MPI] Best Fitness: " << out_p.f_val
                << " | FEvals: " << stop2.get_current_fevals()
                << " | Iters: " << out_p.iterations
                << " | MPI Time: " << (t_end - t_start) << "s\n\n";
    }
  }

  MPI_Finalize();
  return 0;
}
