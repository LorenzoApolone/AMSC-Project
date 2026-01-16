#include "functions.cpp"
#include "methods.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

int main(int argc, char **argv)
{
  // Gathering input parameters

  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0]
              << " <dim> <n_points> <max_iter> <delta_x>\n";
    return 1;
  }

  MPI_Init(&argc, &argv);
  unsigned int dim = atoi(argv[1]);
  unsigned int n_points = atoi(argv[2]);
  unsigned int max_iter = atoi(argv[3]);
  double delta_x = atof(argv[4]);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Factory Definition
  std::unordered_map<std::string,
                     std::function<std::unique_ptr<TestFunction>(unsigned int)>>
      factory;

  factory["Sphere"] = [](unsigned int dim)
  {
    return std::make_unique<Sphere>(dim);
  };
  factory["Ellipsoid"] = [](unsigned int dim)
  {
    return std::make_unique<Ellipsoid>(dim);
  };
  factory["SumOfDiffPowers"] = [](unsigned int dim)
  {
    return std::make_unique<SumOfDiffPowers>(dim);
  };
  factory["DropWave"] = [](unsigned int dim)
  {
    return std::make_unique<DropWave>(dim);
  };
  factory["Weierstrass"] = [](unsigned int dim)
  {
    return std::make_unique<Weierstrass>(dim);
  };
  factory["Alpine1"] = [](unsigned int dim)
  {
    return std::make_unique<Alpine1>(dim);
  };
  factory["Ackley"] = [](unsigned int dim)
  {
    return std::make_unique<Ackley>(dim);
  };
  factory["Griewank"] = [](unsigned int dim)
  {
    return std::make_unique<Griewank>(dim);
  };
  factory["Rastrigin"] = [](unsigned int dim)
  {
    return std::make_unique<Rastrigin>(dim);
  };
  factory["HappyCat"] = [](unsigned int dim)
  {
    return std::make_unique<HappyCat>(dim);
  };
  factory["HGBat"] = [](unsigned int dim)
  {
    return std::make_unique<HGBat>(dim);
  };
  factory["Rosenbrock"] = [](unsigned int dim)
  {
    return std::make_unique<Rosenbrock>(dim);
  };
  factory["HighCondElliptic"] = [](unsigned int dim)
  {
    return std::make_unique<HighCondElliptic>(dim);
  };
  factory["Discus"] = [](unsigned int dim)
  {
    return std::make_unique<Discus>(dim);
  };
  factory["BentCigar"] = [](unsigned int dim)
  {
    return std::make_unique<BentCigar>(dim);
  };
  factory["PermdbFunc"] = [](unsigned int dim)
  {
    return std::make_unique<PermdbFunc>(dim);
  };
  factory["Schafferf7Func"] = [](unsigned int dim)
  {
    return std::make_unique<Schafferf7Func>(dim);
  };
  factory["ExpSchafferF6"] = [](unsigned int dim)
  {
    return std::make_unique<ExpSchafferF6>(dim);
  };
  factory["RotatedHyper"] = [](unsigned int dim)
  {
    return std::make_unique<RotatedHyper>(dim);
  };
  factory["Schwefel"] = [](unsigned int dim)
  {
    return std::make_unique<Schwefel>(dim);
  };
  factory["SumOfDifferentPowers2"] = [](unsigned int dim)
  {
    return std::make_unique<SumOfDifferentPowers2>(dim);
  };
  factory["XinSheYang1"] = [](unsigned int dim)
  {
    return std::make_unique<XinSheYang1>(dim);
  };
  factory["Schwefel221"] = [](unsigned int dim)
  {
    return std::make_unique<Schwefel221>(dim);
  };
  factory["Schwefel222"] = [](unsigned int dim)
  {
    return std::make_unique<Schwefel222>(dim);
  };
  factory["Salomon"] = [](unsigned int dim)
  {
    return std::make_unique<Salomon>(dim);
  };
  factory["ModifiedRidge"] = [](unsigned int dim)
  {
    return std::make_unique<ModifiedRidge>(dim);
  };
  factory["Zakharov"] = [](unsigned int dim)
  {
    return std::make_unique<Zakharov>(dim);
  };
  factory["ModifiedXinSheYang3"] = [](unsigned int dim)
  {
    return std::make_unique<ModifiedXinSheYang3>(dim);
  };
  factory["ModifiedXinSheYang5"] = [](unsigned int dim)
  {
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

  StopCriterion stop(max_iter, delta_x);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Run the solver
  for (const auto &name : function_names)
  {
    auto f_ptr = factory[name](dim);
    OutputObject result = pso_mpi(*f_ptr, dim, stop, n_points);
    if (rank == 0)
    {
      result.terminal_info();
      result.output_to_file();
    }
  }

  MPI_Finalize();
  return 0;
}