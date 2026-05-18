#include "dpso/methods_dpso.hpp"
#include "functions.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

std::unordered_map<std::string, DPSOParameters>
parse_params(const std::string &filename) {
  std::unordered_map<std::string, DPSOParameters> config;
  std::ifstream file(filename);
  if (!file.is_open()) {
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
        config[current_section] =
            config.count("GLOBAL") ? config["GLOBAL"] : DPSOParameters();
      }
      continue;
    }

    std::istringstream iss(line);
    std::string key;
    double val;
    if (iss >> key >> val) {
      auto &params = config[current_section];
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
  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
              << " <dim> <total_particles> <max_iter> [convergence_tol] "
                 "[params_file] [seed]\n";
    return 1;
  }

  unsigned int dim = std::atoi(argv[1]);
  unsigned int total_particles = std::atoi(argv[2]);
  int base_iter = std::atoi(argv[3]);

  std::unordered_map<std::string,
                     std::function<std::unique_ptr<TestFunction>(unsigned int)>>
      factory;
  // clang-format off
  factory["Sphere"] = [](unsigned int dim) { return std::make_unique<Sphere>(dim); };
  factory["Ellipsoid"] = [](unsigned int dim) { return std::make_unique<Ellipsoid>(dim); };
  factory["SumOfDiffPowers"] = [](unsigned int dim) { return std::make_unique<SumOfDiffPowers>(dim); };
  factory["QuinticFunction"] = [](unsigned int dim) { return std::make_unique<QuinticFunction>(dim); };
  factory["DropWave"] = [](unsigned int dim) { return std::make_unique<DropWave>(dim); };
  factory["Weierstrass"] = [](unsigned int dim) { return std::make_unique<Weierstrass>(dim); };
  factory["Alpine1"] = [](unsigned int dim) { return std::make_unique<Alpine1>(dim); };
  factory["Ackley"] = [](unsigned int dim) { return std::make_unique<Ackley>(dim); };
  factory["Griewank"] = [](unsigned int dim) { return std::make_unique<Griewank>(dim); };
  factory["Rastrigin"] = [](unsigned int dim) { return std::make_unique<Rastrigin>(dim); };
  factory["HappyCat"] = [](unsigned int dim) { return std::make_unique<HappyCat>(dim); };
  factory["HGBat"] = [](unsigned int dim) { return std::make_unique<HGBat>(dim); };
  factory["Rosenbrock"] = [](unsigned int dim) { return std::make_unique<Rosenbrock>(dim); };
  factory["HighCondElliptic"] = [](unsigned int dim) { return std::make_unique<HighCondElliptic>(dim); };
  factory["Discus"] = [](unsigned int dim) { return std::make_unique<Discus>(dim); };
  factory["BentCigar"] = [](unsigned int dim) { return std::make_unique<BentCigar>(dim); };
  factory["Schafferf7Func"] = [](unsigned int dim) { return std::make_unique<Schafferf7Func>(dim); };
  factory["ExpSchafferF6"] = [](unsigned int dim) { return std::make_unique<ExpSchafferF6>(dim); };
  factory["RotatedHyper"] = [](unsigned int dim) { return std::make_unique<RotatedHyper>(dim); };
  factory["Schwefel"] = [](unsigned int dim) { return std::make_unique<Schwefel>(dim); };
  factory["SumOfDifferentPowers2"] = [](unsigned int dim) { return std::make_unique<SumOfDifferentPowers2>(dim); };
  factory["XinSheYang1"] = [](unsigned int dim) { return std::make_unique<XinSheYang1>(dim); };
  factory["Schwefel221"] = [](unsigned int dim) { return std::make_unique<Schwefel221>(dim); };
  factory["Schwefel222"] = [](unsigned int dim) { return std::make_unique<Schwefel222>(dim); };
  factory["Salomon"] = [](unsigned int dim) { return std::make_unique<Salomon>(dim); };
  factory["ModifiedRidge"] = [](unsigned int dim) { return std::make_unique<ModifiedRidge>(dim); };
  factory["Zakharov"] = [](unsigned int dim) { return std::make_unique<Zakharov>(dim); };
  factory["ModifiedXinSheYang3"] = [](unsigned int dim) { return std::make_unique<ModifiedXinSheYang3>(dim); };
  factory["ModifiedXinSheYang5"] = [](unsigned int dim) { return std::make_unique<ModifiedXinSheYang5>(dim); };
  factory["Levy"] = [](unsigned int dim) { return std::make_unique<Levy>(dim); };
  factory["Michalewicz"] = [](unsigned int dim) { return std::make_unique<Michalewicz>(dim); };
  factory["Bohachevsky"] = [](unsigned int dim) { return std::make_unique<Bohachevsky>(dim); };
  factory["Powell"] = [](unsigned int dim) { return std::make_unique<Powell>(dim); };
  factory["DixonPrice"] = [](unsigned int dim) { return std::make_unique<DixonPrice>(dim); };
  factory["StyblinskiTang"] = [](unsigned int dim) { return std::make_unique<StyblinskiTang>(dim); };
  factory["Step"]           = [](unsigned int dim) { return std::make_unique<Step>(dim); };
  factory["Qing"]           = [](unsigned int dim) { return std::make_unique<Qing>(dim); };
  factory["Trid"]           = [](unsigned int dim) { return std::make_unique<Trid>(dim); };
  factory["Shubert"]        = [](unsigned int dim) { return std::make_unique<Shubert>(dim); };
  factory["Alpine2"]        = [](unsigned int dim) { return std::make_unique<Alpine2>(dim); };
  factory["Eggholder"]      = [](unsigned int dim) { return std::make_unique<Eggholder>(dim); };
  factory["Easom"]          = [](unsigned int dim) { return std::make_unique<Easom>(dim); };
  factory["Brown"]          = [](unsigned int dim) { return std::make_unique<Brown>(dim); };
  factory["Csendes"]        = [](unsigned int dim) { return std::make_unique<Csendes>(dim); };
  factory["Vincent"]        = [](unsigned int dim) { return std::make_unique<Vincent>(dim); };
  // clang-format on

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
                                             "StyblinskiTang",
                                             "Step",
                                             "Qing",
                                             "Trid",
                                             "Shubert",
                                             "Alpine2",
                                             "Eggholder",
                                             "Easom",
                                             "Brown",
                                             "Csendes",
                                             "Vincent"};

  int number_of_functions = function_names.size();
  int iters = base_iter;
  unsigned int ppr = total_particles;
  double true_convergence_tol = 1e-4;
  std::string param_file = "";
  unsigned int seed = 0;
  if (argc >= 5) {
    true_convergence_tol = std::atof(argv[4]);
  }
  if (argc >= 6) {
    param_file = argv[5];
  }
  if (argc >= 7) {
    seed = static_cast<unsigned int>(std::stoul(argv[6]));
  }

  std::unordered_map<std::string, DPSOParameters> all_params;
  if (!param_file.empty()) {
    all_params = parse_params(param_file);
  }

  using StrSet = std::unordered_map<std::string, std::vector<std::string>>;
  StrSet categories;
  categories["Unimodal"] = {"Sphere",
                            "Ellipsoid",
                            "SumOfDiffPowers",
                            "QuinticFunction",
                            "HighCondElliptic",
                            "Discus",
                            "BentCigar",
                            "RotatedHyper",
                            "Zakharov",
                            "Powell",
                            "DixonPrice",
                            "Qing",
                            "Trid",
                            "Brown"};
  categories["Multimodal"] = {"DropWave",
                              "Weierstrass",
                              "Alpine1",
                              "Ackley",
                              "Griewank",
                              "Rastrigin",
                              "HappyCat",
                              "HGBat",
                              "Rosenbrock",
                              "Schafferf7Func",
                              "ExpSchafferF6",
                              "Schwefel",
                              "SumOfDifferentPowers2",
                              "XinSheYang1",
                              "Schwefel221",
                              "Schwefel222",
                              "Salomon",
                              "ModifiedRidge",
                              "ModifiedXinSheYang3",
                              "ModifiedXinSheYang5",
                              "Levy",
                              "Michalewicz",
                              "Bohachevsky",
                              "StyblinskiTang",
                              "Shubert",
                              "Alpine2",
                              "Eggholder",
                              "Easom",
                              "Csendes",
                              "Vincent"};
  categories["Separable"] = {
      "Sphere",      "SumOfDiffPowers", "QuinticFunction",
      "Weierstrass", "Alpine1",         "Ackley",
      "Rastrigin",   "Schwefel",        "XinSheYang1",
      "Schwefel221", "Schwefel222",     "Salomon",
      "Michalewicz", "Bohachevsky",     "StyblinskiTang",
      "Levy",        "DixonPrice",      "Step",
      "Qing",        "Alpine2",         "Easom",
      "Csendes",     "Vincent"};
  categories["Non-separable"] = {"Ellipsoid",
                                 "DropWave",
                                 "Griewank",
                                 "HappyCat",
                                 "HGBat",
                                 "Rosenbrock",
                                 "HighCondElliptic",
                                 "Discus",
                                 "BentCigar",
                                 "Schafferf7Func",
                                 "ExpSchafferF6",
                                 "RotatedHyper",
                                 "SumOfDifferentPowers2",
                                 "ModifiedRidge",
                                 "Zakharov",
                                 "ModifiedXinSheYang3",
                                 "ModifiedXinSheYang5",
                                 "Powell",
                                 "Trid",
                                 "Shubert",
                                 "Eggholder",
                                 "Brown"};
  categories["Differentiable"] = {"Sphere",
                                  "Ellipsoid",
                                  "SumOfDiffPowers",
                                  "QuinticFunction",
                                  "Griewank",
                                  "HappyCat",
                                  "HGBat",
                                  "Rosenbrock",
                                  "HighCondElliptic",
                                  "Discus",
                                  "BentCigar",
                                  "RotatedHyper",
                                  "Schwefel",
                                  "SumOfDifferentPowers2",
                                  "Schwefel221",
                                  "Schwefel222",
                                  "Salomon",
                                  "ModifiedRidge",
                                  "Zakharov",
                                  "ModifiedXinSheYang3",
                                  "Levy",
                                  "Michalewicz",
                                  "Bohachevsky",
                                  "Powell",
                                  "DixonPrice",
                                  "StyblinskiTang",
                                  "Rastrigin",
                                  "Ackley",
                                  "Qing",
                                  "Trid",
                                  "Shubert",
                                  "Alpine2",
                                  "Eggholder",
                                  "Easom",
                                  "Brown",
                                  "Csendes",
                                  "Vincent"};
  categories["Non-differentiable"] = {
      "Weierstrass",   "Alpine1",     "Schafferf7Func",
      "ExpSchafferF6", "XinSheYang1", "ModifiedXinSheYang5",
      "Step"};
  categories["Flat"] = {"DropWave", "Step"};
  categories["Coupled"] = {"HappyCat", "HGBat", "Rosenbrock"};

  // Stopping counters (same 4-bucket style as ExperimentStats in topology)
  int stopped_by_maxiter_and_incorrect = 0;
  int stopped_by_maxiter_and_correct = 0;
  int incorrect_when_early_stop = 0;
  int correct_when_early_stop = 0;
  int correct_total = 0;
  std::vector<std::string> functions_converged;

  auto t_start = std::chrono::high_resolution_clock::now();

  for (const auto &name : function_names) {
    auto f_ptr = factory[name](dim);

    DPSOParameters current_params;
    if (all_params.count(name)) {
      current_params = all_params[name];
    } else if (all_params.count("GLOBAL")) {
      current_params = all_params["GLOBAL"];
    }

    try {
      OutputObject res =
          dpso_serial(*f_ptr, dim, ppr, iters, current_params, 1e-6, seed);
      double fval = f_ptr->value(res.x_best);
      double err = f_ptr->error(res.x_best);
      bool true_converged = err < true_convergence_tol;

      std::cout << std::left << std::setw(22) << name << std::right
                << std::setw(6) << dim << std::setw(8) << ppr << std::setw(8)
                << iters << "   " << std::scientific << std::setprecision(4)
                << std::setw(13) << fval << "   " << std::setw(13) << err
                << "   " << std::fixed << std::setprecision(6) << std::setw(10)
                << res.execution_time << "s" << std::endl;

      std::cout << "Comm times: Total 0.0 s | Compute " << std::fixed << std::setprecision(6) << res.execution_time 
                << " s | Bcast 0.0 s | Allreduce 0.0 s | Allgather 0.0 s | Wait 0.0 s" << std::endl;

      bool stopped_by_maxiter = (res.iterations >= iters);
      if (stopped_by_maxiter && !true_converged)
        stopped_by_maxiter_and_incorrect++;
      if (stopped_by_maxiter && true_converged)
        stopped_by_maxiter_and_correct++;
      if (!stopped_by_maxiter && !true_converged)
        incorrect_when_early_stop++;
      if (!stopped_by_maxiter && true_converged)
        correct_when_early_stop++;
      if (true_converged) {
        correct_total++;
        functions_converged.push_back(name);
      }
    } catch (const std::invalid_argument &e) {
      std::cerr << "Configuration Error for " << name << ": " << e.what()
                << "\n";
      return 1;
    }
  }
  std::cout << std::string(90, '-') << std::endl;

  auto t_end = std::chrono::high_resolution_clock::now();
  double total_time = std::chrono::duration<double>(t_end - t_start).count();

  // --- machine-readable summary line (same format as topology) ---
  std::cout << "\nRESULT,"
            << "version=dpso_serial,"
            << "time=" << total_time << ","
            << "conv=" << correct_total << ","
            << "total=" << number_of_functions << ","
            << "\n";

  // --- human-readable summary (same style as ExperimentStats in topology) ---
  std::cout << "\n=== ALL BENCHMARKS COMPLETED ===" << std::endl;
  std::cout << "Total execution time: " << total_time << " s" << std::endl;
  std::cout << "Stopped by max iter and incorrect: "
            << stopped_by_maxiter_and_incorrect << std::endl;
  std::cout << "Stopped by max iter and correct: "
            << stopped_by_maxiter_and_correct << std::endl;
  std::cout << "Incorrect when early stop: " << incorrect_when_early_stop
            << std::endl;
  std::cout << "Correct when early stop: " << correct_when_early_stop
            << std::endl;
  std::cout << "Correct total: " << correct_total << "/" << number_of_functions
            << std::endl;

  // --- per-category convergence breakdown ---
  std::unordered_set<std::string> converged_set(functions_converged.begin(),
                                                functions_converged.end());
  std::cout << "Converged by typology:" << std::endl;
  for (const auto &cat :
       {std::string("Unimodal"), std::string("Multimodal"),
        std::string("Separable"), std::string("Non-separable"),
        std::string("Differentiable"), std::string("Non-differentiable"),
        std::string("Flat"), std::string("Coupled")}) {
    const auto &cat_fns = categories.at(cat);
    int cat_total = 0;
    int cat_conv = 0;
    for (const auto &fn : cat_fns) {
      // only count functions that were actually tested
      bool in_test = (std::find(function_names.begin(), function_names.end(),
                                fn) != function_names.end());
      if (in_test) {
        cat_total++;
        if (converged_set.count(fn))
          cat_conv++;
      }
    }
    std::cout << "  - " << cat << ": " << cat_conv << "/" << cat_total
              << std::endl;
  }

  // --- functions that did not converge ---
  std::cout << std::endl;
  for (const auto &fn : function_names) {
    if (!converged_set.count(fn))
      std::cout << "Function " << fn << " did not converge." << std::endl;
  }

  return 0;
}
