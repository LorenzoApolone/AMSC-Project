#include "dpso/methods_dpso.hpp"
#include "functions.hpp"
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

std::unordered_map<std::string, DPSOParameters> parse_params(const std::string &filename, int rank) {
  std::unordered_map<std::string, DPSOParameters> config;
  std::ifstream file(filename);
  if (!file.is_open()) {
    if (rank == 0)
      std::cerr << "Warning: Could not open parameter file " << filename
                << ". Using defaults.\n";
    return config;
  }
  std::string line;
  std::string current_section = "GLOBAL";
  config[current_section] = DPSOParameters();

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#')
      continue;
      
    if (line[0] == '[' && line.back() == ']') {
        current_section = line.substr(1, line.size() - 2);
        if (config.find(current_section) == config.end()) {
            config[current_section] = config.count("GLOBAL") ? config["GLOBAL"] : DPSOParameters();
        }
        continue;
    }

    std::istringstream iss(line);
    std::string key;
    double val;
    if (iss >> key >> val) {
      auto& params = config[current_section];
      if (key == "w")
        params.w = val;
      else if (key == "c1")
        params.c1 = val;
      else if (key == "c2")
        params.c2 = val;
      else if (key == "regrouping_period")
        params.regrouping_period = (int)val;
      else if (key == "sub_swarm_size")
        params.sub_swarm_size = (int)val;
      else if (key == "hmcr")
        params.hmcr = val;
      else if (key == "par_min")
        params.par_min = val;
      else if (key == "par_max")
        params.par_max = val;
    }
  }
  return config;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double start_time = MPI_Wtime();

  if (argc < 4) {
    if (rank == 0) {
      std::cout << "Usage: " << argv[0]
                << " <dim> <total_particles> <max_iter> [convergence_tol] "
                   "[params_file]\n";
    }
    MPI_Finalize();
    return 1;
  }

  unsigned int dim = std::atoi(argv[1]);
  unsigned int total_particles = std::atoi(argv[2]);
  int base_iter = std::atoi(argv[3]);

  // Factory for all benchmark functions
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
  factory["QuinticFunction"] = [](unsigned int dim) {
    return std::make_unique<QuinticFunction>(dim);
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
  factory["Levy"] = [](unsigned int dim) {
    return std::make_unique<Levy>(dim);
  };
  factory["Michalewicz"] = [](unsigned int dim) {
    return std::make_unique<Michalewicz>(dim);
  };
  factory["Bohachevsky"] = [](unsigned int dim) {
    return std::make_unique<Bohachevsky>(dim);
  };
  factory["Powell"] = [](unsigned int dim) {
    return std::make_unique<Powell>(dim);
  };
  factory["DixonPrice"] = [](unsigned int dim) {
    return std::make_unique<DixonPrice>(dim);
  };
  factory["StyblinskiTang"] = [](unsigned int dim) {
    return std::make_unique<StyblinskiTang>(dim);
  };

  std::vector<std::string> function_names = {"Sphere",
                                             "Ellipsoid",
                                             "SumOfDiffPowers",
                                             "QuinticFunction",
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

  int number_of_converged = 0;
  int number_of_functions = function_names.size();
  std::vector<std::string> functions_converged;
  double t_start = MPI_Wtime();

  int iters = base_iter;
  unsigned int ppr = total_particles;
  double true_convergence_tol = 1e-4;
  std::string param_file = "";
  if (argc >= 5) {
    true_convergence_tol = std::atof(argv[4]);
  }
  if (argc >= 6) {
    param_file = argv[5];
  }
  std::unordered_map<std::string, DPSOParameters> all_params;
  if (!param_file.empty()) {
    all_params = parse_params(param_file, rank);
  }
  int stopped_by_max_iter = 0;
  int incorrect_early_stop = 0;
  int correct_total = 0;
  for (const auto &name : function_names) {
    bool converged = false;
    auto f_ptr = factory[name](dim);
    
    DPSOParameters current_params;
    if (all_params.count(name)) {
        current_params = all_params[name];
    } else if (all_params.count("GLOBAL")) {
        current_params = all_params["GLOBAL"];
    }

    OutputObject res = dpso(*f_ptr, dim, ppr, iters, current_params);
    double fval = f_ptr->value(res.x_best);
    double err = f_ptr->error(res.x_best);
    bool true_converged = err < true_convergence_tol;
    if (rank == 0) {
      std::cout << std::left << std::setw(22) << name << std::right
                << std::setw(6) << dim << std::setw(8) << ppr // total particles
                << std::setw(8) << iters << "   " << std::scientific
                << std::setprecision(4) << std::setw(13) << fval << "   "
                << std::setw(13) << err << "   " << std::fixed
                << std::setprecision(2) << std::setw(8) << res.execution_time
                << "s" << std::endl;
      if (res.iterations >= iters)
        stopped_by_max_iter++;
      else if (!true_converged)
        incorrect_early_stop++;
      if (true_converged)
        correct_total++;
    }
  }
  if (rank == 0)
    std::cout << std::string(90, '-') << std::endl;

  double t_end = MPI_Wtime();
  if (rank == 0) {
    std::cout << "\nRESULT,"
              << "version=dpso,"
              << "time=" << (t_end - t_start) << ","
              << "conv=" << correct_total << "," // <-- usa correct_total!
              << "total=" << number_of_functions << ","
              << "\n";
    std::cout << "\n=== ALL BENCHMARKS COMPLETED ===" << std::endl;
    std::cout << "Tempo totale esecuzione: " << (t_end - t_start) << " s"
              << std::endl;
    std::cout << "Convergence rate: " << correct_total << "/"
              << number_of_functions << std::endl;
    std::cout << "Stopped by max iter: " << stopped_by_max_iter << std::endl;
    std::cout << "Incorrect when early stop: " << incorrect_early_stop
              << std::endl;
    std::cout << "Correct total: " << correct_total << std::endl;
  }

  MPI_Finalize();
  return 0;
}
