/**
 * @file main_topology.cpp
 * @brief Main file for topology-based PSO experiments using MPI.
 *
 * This file runs and compares multiple PSO variants on a set of benchmark
 * functions:
 * - topology-based PSO with small-world communication graph
 * - topology-based PSO with scale-free communication graph
 * - topology-based PSO with random communication graph
 * - classic MPI PSO without explicit communication topology
 * 
 * Using: 
 * - the function in pso_topology.hpp, which implements the PSO algorithm with a neighborhood topology. 
 * - the function in create_network.hpp, which creates the communication topologies
 * - the stopping criteria defined in StoppingCriteriaManager.hpp
 * 
 * The program:
 * 1. initializes MPI,
 * 2. parses command-line parameters,
 * 3. initializes the functions,
 * 4. runs the selected PSO variants,
 * 5. collects convergence annd time,
 * 6. prints machine-readable and human-readable summaries.
 */

#include "create_network.hpp"
#include "confront.hpp"
#include "pso_topology.hpp"
#include "../methods.hpp"
#include "../functions.hpp" 
#include "../interfaces/StoppingCriteriaManager.hpp"
#include <mpi.h>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

enum class TopologyMode {
    SMALL_WORLD,
    SCALE_FREE,
    RANDOM
};
/**
 * @brief Broadcasts an adjacency list containing the communication graph from rank 0 to all MPI ranks.
 *
 *
 * @param[in,out] adjacency_list Adjacency list to broadcast. On rank 0 it must
 * contain the valid graph.
 * @param[in] n_points Number of graph nodes / particles.
 * @param[in] rank Rank of the calling MPI process.
 *
 */
static void bcast_adjacency_list(std::vector<std::vector<int>>& adjacency_list,
                                int n_points,
                                int rank)
{
  std::vector<int> degrees(n_points);
  std::vector<int> flat;
 
  if (rank == 0) {
    
    for (int i = 0; i < n_points; ++i) {
      degrees[i] = static_cast<int>(adjacency_list[i].size());
      flat.insert(flat.end(), adjacency_list[i].begin(), adjacency_list[i].end());
    }
  }

  MPI_Bcast(degrees.data(), n_points, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> offsets(n_points);
  offsets[0] = 0;
  for (int i = 1; i < n_points; ++i)
    offsets[i] = offsets[i - 1] + degrees[i - 1];

  int total = offsets[n_points - 1] + degrees[n_points - 1];

  if (rank != 0) flat.resize(total);

  MPI_Bcast(flat.data(), total, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) adjacency_list.resize(n_points);
  for (int i = 0; i < n_points; ++i) {
    adjacency_list[i].assign(flat.begin() + offsets[i],
                             flat.begin() + offsets[i] + degrees[i]);
  }
}

/**
 * @struct ExperimentStats
 * @brief Stores convergence and timing statistics for one experiment.
 *
 * This structure collects:
 * - stopping-condition statistics,
 * - number of converged benchmark functions,
 * - communication time spent in MPI_Allgatherv,
 * - total execution time of the experiment,
 * - names of functions for which convergence was achieved.
 */
struct ExperimentStats {
    int stopped_by_maxiter_and_incorrect = 0;
    int stopped_by_maxiter_and_correct = 0;
    int incorrect_when_early_stop = 0;
    int correct_when_early_stop = 0;
    int correct_total = 0;
    int number_of_converged = 0;
    double t_allgatherv = 0.0;
    double total_time = 0.0;
    std::vector<std::string> functions_converged;
};

/**
 * @brief Updates the statistics of one experiment.
 *
 * The function compares the final fitness returned by the optimizer with the
 * true optimum value of the function, if the module of the difference is less than 
 * delta_x, the result is considered correct. 
 *
 * @param[in] name Name of the function.
 * @param[in] result Output object returned by the optimizer.
 * @param[in] function  Function associated with the run.
 * @param[in] delta_x Tolerance used to decide whether the final result is correct.
 * @param[in] max_iter Maximum allowed number of iterations.
 * @param[in,out] stats Statistics object to update.
 */
static void update_experiment_stats(const std::string& name,
                                    const OutputObject& result,
                                    const TestFunction& function,
                                    double delta_x,
                                    unsigned int max_iter,
                                    ExperimentStats& stats)
{
    const double final_fitness = result.get_best_fitness();
    const double f_star = function.value(function.get_true_solution());

    const bool is_correct = (std::abs(final_fitness - f_star) <= delta_x);
    const bool stopped_by_maxiter = (result.iterations >= static_cast<int>(max_iter));

    if (is_correct) {
        stats.correct_total++;
    }

    if (stopped_by_maxiter && !is_correct) {
        stats.stopped_by_maxiter_and_incorrect++;
    }

    if (stopped_by_maxiter && is_correct) {
        stats.stopped_by_maxiter_and_correct++;
        stats.number_of_converged++;
        stats.functions_converged.push_back(name);
    }

    if (!stopped_by_maxiter && !is_correct) {
        stats.incorrect_when_early_stop++;
    }

    if (!stopped_by_maxiter && is_correct) {
        stats.correct_when_early_stop++;
        stats.number_of_converged++;
        stats.functions_converged.push_back(name);
    }
}
/**
 * @brief Prints a human-readable summary of experiment statistics.
 *
 * @param[in] label Name of the experiment to display.
 * @param[in] stats Statistics collected for that experiment.
 */
static void print_experiment_stats(const std::string& label,
                                   const ExperimentStats& stats)
{
    std::cout << label << ":\n";
    std::cout << "Stopped by max iter and incorrect: "
              << stats.stopped_by_maxiter_and_incorrect << std::endl;
    std::cout << "Stopped by max iter and correct: "
              << stats.stopped_by_maxiter_and_correct << std::endl;
    std::cout << "Incorrect when early stop: "
              << stats.incorrect_when_early_stop << std::endl;
    std::cout << "Correct when early stop: "
              << stats.correct_when_early_stop << std::endl;
    std::cout << "Correct total: "
              << stats.correct_total << std::endl;
}
/**
 * @brief Inizializing all the function.
 *
 *
 * @return A map from function names to construction lambdas.
 */
using FunctionFactory =
    std::unordered_map<std::string,
    std::function<std::unique_ptr<TestFunction>(unsigned int)>>;
static FunctionFactory build_factory()
{
    FunctionFactory factory;
    factory["Sphere"] = [](unsigned int dim){ return std::make_unique<Sphere>(dim); };
    factory["Ellipsoid"] = [](unsigned int dim){ return std::make_unique<Ellipsoid>(dim); };
    factory["SumOfDiffPowers"] = [](unsigned int dim){ return std::make_unique<SumOfDiffPowers>(dim); };
    factory["DropWave"] = [](unsigned int dim){ return std::make_unique<DropWave>(dim); };
    factory["Weierstrass"] = [](unsigned int dim){ return std::make_unique<Weierstrass>(dim); };
    factory["Alpine1"] = [](unsigned int dim){ return std::make_unique<Alpine1>(dim); };
    factory["Ackley"] = [](unsigned int dim){ return std::make_unique<Ackley>(dim); };
    factory["Griewank"] = [](unsigned int dim){ return std::make_unique<Griewank>(dim); };
    factory["Rastrigin"] = [](unsigned int dim){ return std::make_unique<Rastrigin>(dim); };
    factory["HappyCat"] = [](unsigned int dim){ return std::make_unique<HappyCat>(dim); };
    factory["HGBat"] = [](unsigned int dim){ return std::make_unique<HGBat>(dim); };
    factory["Rosenbrock"] = [](unsigned int dim){ return std::make_unique<Rosenbrock>(dim); };
    factory["HighCondElliptic"] = [](unsigned int dim){ return std::make_unique<HighCondElliptic>(dim); };
    factory["Discus"] = [](unsigned int dim){ return std::make_unique<Discus>(dim); };
    factory["BentCigar"] = [](unsigned int dim){ return std::make_unique<BentCigar>(dim); };
    factory["PermdbFunc"] = [](unsigned int dim){ return std::make_unique<PermdbFunc>(dim); };
    factory["Schafferf7Func"] = [](unsigned int dim){ return std::make_unique<Schafferf7Func>(dim); };
    factory["ExpSchafferF6"] = [](unsigned int dim){ return std::make_unique<ExpSchafferF6>(dim); };
    factory["RotatedHyper"] = [](unsigned int dim){ return std::make_unique<RotatedHyper>(dim); };
    factory["Schwefel"] = [](unsigned int dim){ return std::make_unique<Schwefel>(dim); };
    factory["SumOfDifferentPowers2"] = [](unsigned int dim){ return std::make_unique<SumOfDifferentPowers2>(dim); };
    factory["XinSheYang1"] = [](unsigned int dim){ return std::make_unique<XinSheYang1>(dim); };
    factory["Schwefel221"] = [](unsigned int dim){ return std::make_unique<Schwefel221>(dim); };
    factory["Schwefel222"] = [](unsigned int dim){ return std::make_unique<Schwefel222>(dim); };
    factory["Salomon"] = [](unsigned int dim){ return std::make_unique<Salomon>(dim); };
    factory["ModifiedRidge"] = [](unsigned int dim){ return std::make_unique<ModifiedRidge>(dim); };
    factory["Zakharov"] = [](unsigned int dim){ return std::make_unique<Zakharov>(dim); };
    factory["ModifiedXinSheYang3"] = [](unsigned int dim){ return std::make_unique<ModifiedXinSheYang3>(dim); };
    factory["ModifiedXinSheYang5"] = [](unsigned int dim){ return std::make_unique<ModifiedXinSheYang5>(dim); };
    factory["Levy"] = [](unsigned int dim){ return std::make_unique<Levy>(dim); };
    factory["Michalewicz"] = [](unsigned int dim){ return std::make_unique<Michalewicz>(dim); };
    factory["Bohachevsky"] = [](unsigned int dim){ return std::make_unique<Bohachevsky>(dim); };
    factory["Powell"] = [](unsigned int dim){ return std::make_unique<Powell>(dim); };
    factory["DixonPrice"] = [](unsigned int dim){ return std::make_unique<DixonPrice>(dim); };
    factory["StyblinskiTang"] = [](unsigned int dim){ return std::make_unique<StyblinskiTang>(dim); };
    return factory;
}
/**
 * @brief Returns the list of function names used in the experiments.
 *
 * @return Vector containing the names of all benchmark functions to test.
 */
static std::vector<std::string> build_function_names()
{
    return {
        "Sphere","Ellipsoid","SumOfDiffPowers","DropWave","Weierstrass","Alpine1","Ackley",
        "Griewank","Rastrigin","HappyCat","HGBat","Rosenbrock","HighCondElliptic","Discus",
        "BentCigar","PermdbFunc","Schafferf7Func","ExpSchafferF6","RotatedHyper","Schwefel",
        "SumOfDifferentPowers2","XinSheYang1","Schwefel221","Schwefel222","Salomon",
        "ModifiedRidge","Zakharov","ModifiedXinSheYang3","ModifiedXinSheYang5",
        "Levy","Michalewicz","Bohachevsky","Powell","DixonPrice","StyblinskiTang"
    };
}
/**
 * @brief Runs the topology PSO on all functions.
 *
 * For each benchmark function, the selected topology is generated,
 * broadcast to all MPI ranks, and then used by `pso_topology(...)`.
 *
 * @param[in] mode Topology type to use (small-world, scale-free, or random).
 * @param[in] rank Rank of the calling MPI process.
 * @param[in] dim Dimension of the optimization problem.
 * @param[in] n_points Number of particles in the swarm.
 * @param[in] max_iter Maximum number of iterations.
 * @param[in] delta_x Acceptance tolerance used to determine convergence.
 * @param[in] iterations_stagnation Maximum number of stagnation iterations.
 * @param[in] stagnation_tol Tolerance used for stagnation-based stopping.
 * @param[in] diversity_tol Tolerance used for diversity-based stopping.
 * @param[in] p_rewiring Rewiring probability for small-world topology.
 * @param[in] p_random Edge probability for random topology.
 * @param[in] m Number of edges added by each new node in the scale-free topology.
 * @param[in] function_names List of benchmark-function names.
 * @param[in] factory Factory used to instantiate benchmark functions.
 *
 * @return Statistics by the topology-based experiment.
 */
static ExperimentStats run_topology_experiment(
    TopologyMode mode,
    int rank,
    unsigned int dim,
    unsigned int n_points,
    unsigned int max_iter,
    double delta_x,
    int iterations_stagnation,
    double stagnation_tol,
    double diversity_tol,
    double p_rewiring,
    double p_random,
    int m,
    const std::vector<std::string>& function_names,
    const FunctionFactory& factory)
{
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    ExperimentStats stats;

    for (const auto& name : function_names) {
        auto function = factory.at(name)(dim);
        std::vector<std::vector<int>> adjacency_list;

        if (rank == 0) {
            switch (mode) {
                case TopologyMode::SMALL_WORLD:
                    create_network(static_cast<int>(n_points), p_rewiring, adjacency_list);
                    break;
                case TopologyMode::SCALE_FREE:
                    create_scale_free_network(static_cast<int>(n_points), m, adjacency_list);
                    break;
                case TopologyMode::RANDOM:
                    create_random_network(static_cast<int>(n_points), p_random, adjacency_list);
                    break;
            }
        }

        StoppingCriteriaManager stop(max_iter,
                                     iterations_stagnation,
                                     stagnation_tol,
                                     diversity_tol);

        bcast_adjacency_list(adjacency_list, static_cast<int>(n_points), rank);

        OutputObject result = pso_topology(*function,
                                           dim,
                                           stop,
                                           n_points,
                                           adjacency_list,
                                           stats.t_allgatherv);

        if (rank == 0) {
            update_experiment_stats(name, result, *function, delta_x, max_iter, stats);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    stats.total_time = MPI_Wtime() - t_start;

    return stats;
}
/**
 * @brief Runs the MPI PSO experiment on all functions.
 *
 * @param[in] rank Rank of the calling MPI process.
 * @param[in] dim Dimension of the optimization problem.
 * @param[in] n_points Number of particles in the swarm.
 * @param[in] max_iter Maximum number of iterations.
 * @param[in] delta_x Acceptance tolerance used to determine convergence.
 * @param[in] iterations_stagnation Maximum number of stagnation iterations.
 * @param[in] stagnation_tol Tolerance used for stagnation-based stopping.
 * @param[in] diversity_tol Tolerance used for diversity-based stopping.
 * @param[in] function_names List of benchmark-function names.
 * @param[in] factory Factory used to instantiate benchmark functions.
 *
 * @return Statistics collected over the PSO experiment.
 */
static ExperimentStats run_classic_experiment(
    int rank, 
    unsigned int dim,
    unsigned int n_points,
    unsigned int max_iter,
    double delta_x,
    int iterations_stagnation,
    double stagnation_tol,
    double diversity_tol,
    const std::vector<std::string>& function_names,
    const FunctionFactory& factory)
{
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    ExperimentStats stats;

    for (const auto& name : function_names) {
        auto function = factory.at(name)(dim);

        StoppingCriteriaManager stop(max_iter,
                                     iterations_stagnation,
                                     stagnation_tol,
                                     diversity_tol);

        OutputObject result = pso_mpi(*function, dim, stop, n_points);

        if (rank == 0) {
            update_experiment_stats(name, result, *function, delta_x, max_iter, stats);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    stats.total_time = MPI_Wtime() - t_start;

    return stats;
}

/**
 * @brief Runs the serial topology-based PSO experiment on all functions.
 *
 *
 * @param[in] mode Topology type to use.
 * @param[in] dim Dimension of the optimization problem.
 * @param[in] n_points Number of particles in the swarm.
 * @param[in] max_iter Maximum number of iterations.
 * @param[in] delta_x Acceptance tolerance used to determine convergence.
 * @param[in] iterations_stagnation Maximum number of stagnation iterations.
 * @param[in] stagnation_tol Tolerance used for stagnation-based stopping.
 * @param[in] diversity_tol Tolerance used for diversity-based stopping.
 * @param[in] p_rewiring Rewiring probability for small-world topology.
 * @param[in] p_random Edge probability for random topology.
 * @param[in] m Number of edges added by each new node in the scale-free topology.
 * @param[in] function_names List of benchmark-function names.
 * @param[in] factory Factory used to instantiate benchmark functions.
 *
 * @return Statistics collected over the serial topology-based experiment.
 */
static ExperimentStats run_serial_topology_experiment(
    TopologyMode mode,
    unsigned int dim,
    unsigned int n_points,
    unsigned int max_iter,
    double delta_x,
    int iterations_stagnation,
    double stagnation_tol,
    double diversity_tol,
    double p_rewiring,
    double p_random,
    int m,
    const std::vector<std::string>& function_names,
    const FunctionFactory& factory)
{
    double t_start = MPI_Wtime();

    ExperimentStats stats;

    for (const auto& name : function_names) {
        auto function = factory.at(name)(dim);
        std::vector<std::vector<int>> adjacency_list;

        switch (mode) {
            case TopologyMode::SMALL_WORLD:
                create_network(static_cast<int>(n_points), p_rewiring, adjacency_list);
                break;
            case TopologyMode::SCALE_FREE:
                create_scale_free_network(static_cast<int>(n_points), m, adjacency_list);
                break;
            case TopologyMode::RANDOM:
                create_random_network(static_cast<int>(n_points), p_random, adjacency_list);
                break;
        }

        StoppingCriteriaManager stop(max_iter,
                                     iterations_stagnation,
                                     stagnation_tol,
                                     diversity_tol);

        OutputObject result = pso_serial_topology(*function,
                                                  dim,
                                                  stop,
                                                  n_points,
                                                  adjacency_list);

        update_experiment_stats(name, result, *function, delta_x, max_iter, stats);
    }

    stats.total_time = MPI_Wtime() - t_start;
    return stats;
}
/**
 * @brief Main of the benchmark program.
 *
 * The program expects the following command-line arguments:
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 *
 * @return 0 on success, 1 if the input arguments are invalid.
 */
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Argouments: <dim> <n_points> <max_iter> <delta_x>
  if (argc < 5) {
    if (rank == 0) {
      std::cerr << "Usage: " << argv[0]
                << " <dim> <n_points> <max_iter> <delta_x>\n";
    }
    MPI_Finalize();
    return 1;
  }

  unsigned int dim     = std::atoi(argv[1]);
  unsigned int n_points= std::atoi(argv[2]);
  unsigned int max_iter= std::atoi(argv[3]);
  double delta_x       = std::atof(argv[4]);
  int m = 3;                                // for scale-free network, number of edges of each new node
  int iterations_stagnation = max_iter/2.2; // number of iterations for stagnation control
  double stagnation_tol = 1e-10;            // tolerance for stagnation-based stopping
  double diversity_tol = 1e-8;              // diversity tolerance for stopping criteria
  double p_rewiring = 0.05;                 // rewiring probability for small-world network
  double p_random = 0.08;                   // edge probability for random network

  // Factory Definition 
  FunctionFactory factory = build_factory();
  std::vector<std::string> function_names = build_function_names();
  const int number_of_functions = static_cast<int>(function_names.size());
  ExperimentStats small_stats = run_topology_experiment(
      TopologyMode::SMALL_WORLD,
      rank,
      dim,
      n_points,
      max_iter,
      delta_x,
      iterations_stagnation,
      stagnation_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory
  );

  ExperimentStats scale_stats = run_topology_experiment(
      TopologyMode::SCALE_FREE,
      rank,
      dim,
      n_points,
      max_iter,
      delta_x,
      iterations_stagnation,
      stagnation_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory
  );

  ExperimentStats random_stats = run_topology_experiment(
      TopologyMode::RANDOM,
      rank,
      dim,
      n_points,
      max_iter,
      delta_x,
      iterations_stagnation,
      stagnation_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory
  );

  ExperimentStats classic_stats = run_classic_experiment(
      rank,
      dim,
      n_points,
      max_iter,
      delta_x,
      iterations_stagnation,
      stagnation_tol,
      diversity_tol,
      function_names,
      factory
  );

  ExperimentStats serial_small_stats;
  ExperimentStats serial_scale_stats;
  ExperimentStats serial_random_stats;

  if (rank == 0) {
    serial_small_stats = run_serial_topology_experiment(
        TopologyMode::SMALL_WORLD,
        dim,
        n_points,
        max_iter,
        delta_x,
        iterations_stagnation,
        stagnation_tol,
        diversity_tol,
        p_rewiring,
        p_random,
        m,
        function_names,
        factory
    );

    serial_scale_stats = run_serial_topology_experiment(
        TopologyMode::SCALE_FREE,
        dim,
        n_points,
        max_iter,
        delta_x,
        iterations_stagnation,
        stagnation_tol,
        diversity_tol,
        p_rewiring,
        p_random,
        m,
        function_names,
        factory
    );

    serial_random_stats = run_serial_topology_experiment(
        TopologyMode::RANDOM,
        dim,
        n_points,
        max_iter,
        delta_x,
        iterations_stagnation,
        stagnation_tol,
        diversity_tol,
        p_rewiring,
        p_random,
        m,
        function_names,
        factory
    );
  }

  if (rank == 0) {

    std::array<std::vector<std::string>, 5> all = {
        small_stats.functions_converged,
        scale_stats.functions_converged,
        random_stats.functions_converged,
        classic_stats.functions_converged,
        function_names
    };
  
    print_experiment_stats("Small world", small_stats);
    print_experiment_stats("Scale free", scale_stats);
    print_experiment_stats("Random", random_stats);
    print_experiment_stats("Classic", classic_stats);
    print_experiment_stats("Serial small world", serial_small_stats);
    print_experiment_stats("Serial scale free", serial_scale_stats);
    print_experiment_stats("Serial random", serial_random_stats);

    int n = not_converged(all);

  //uniform output format for benchmarking
    std::cout << "\n";
    
    std::cout << "RESULT,"
              << "version=classic,"
              << "time=" << classic_stats.total_time << ","
              << "conv=" << classic_stats.number_of_converged << ","
              << "total=" << number_of_functions << ","
              << "\n";

    
    std::cout << "RESULT,"
              << "version=scale_free,"
              << "time=" << scale_stats.total_time << ","
              << "conv=" << scale_stats.number_of_converged << ","
              << "total=" << number_of_functions << ","
              << "\n";

  
    std::cout << "RESULT," << "version=small_world," << "time=" << small_stats.total_time << ","
              << "conv=" << small_stats.number_of_converged << "," << "total=" << number_of_functions << "," << "\n";

    std::cout << "RESULT," << "version=random," << "time=" << random_stats.total_time << ","
              << "conv=" << random_stats.number_of_converged << "," << "total=" << number_of_functions << "," << "\n";
              
    std::cout << "RESULT," << "not_converged=" << n << "," << "total=" << number_of_functions << ","<< "\n";

    std::cout << "\n";

  // output in human friendly format

    std::cout << "Total time classic PSO: " << classic_stats.total_time << " s\n";
    std::cout << "Convergence rate classic PSO: " << classic_stats.number_of_converged << "/" << number_of_functions  << std::endl;

    std::cout << "Total time small-world network timer version: " << small_stats.t_allgatherv << "/"  << small_stats.total_time << " s\n";
    std::cout << "Convergence rate small-world network: " << small_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time scale-free network timer version: " << scale_stats.t_allgatherv << "/"  << scale_stats.total_time << " s\n";
    std::cout << "Convergence rate scale-free network: " << scale_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time random network timer version: " << random_stats.t_allgatherv << "/"  << random_stats.total_time << " s\n";
    std::cout << "Convergence rate random network: " << random_stats.number_of_converged << "/" << number_of_functions << std::endl << std::endl;
    
    std::cout << "Total time serial small-world topology: " << serial_small_stats.total_time << " s\n";
    std::cout << "Convergence rate serial small-world topology: " << serial_small_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time serial scale-free topology: " << serial_scale_stats.total_time << " s\n";
    std::cout << "Convergence rate serial scale-free topology: " << serial_scale_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time serial random topology: " << serial_random_stats.total_time << " s\n";
    std::cout << "Convergence rate serial random topology: " << serial_random_stats.number_of_converged << "/" << number_of_functions << std::endl;

}
  MPI_Finalize();
  return 0;
}
