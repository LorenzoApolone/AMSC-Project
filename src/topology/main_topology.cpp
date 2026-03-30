/**
 * @file main_topology.cpp
 * @brief Main file for topology-based PSO experiments using MPI.
 *
 * This file execute and compare multiple PSO variants 
 * (Small-World, Scale-Free, Random, and Classic MPI) 
 * over a test suite of benchmark functions available in the "functions.hpp" file.
 * All the PSO variants are executed in parallel and in serial. 
 *
 * All setups, including graph generation, MPI broadcasting of the adjacency lists, 
 * and statistics collections, are taken from the "benchmark_utils.hpp" module.
 * 
 * The program:
 * 1. initializes MPI,
 * 2. parses command-line parameters,
 * 3. imports the function factory from benchmark utilities,
 * 4. runs the selected PSO variants sequentially,
 * 5. collects convergence and timing metrics,
 * 6. prints machine-readable and human-readable summaries.
 */


#include "benchmark_utils.hpp"
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

  // Argouments: <dim> <n_points> <max_iter> <delta_x> [seed]
  if (argc < 5) {
    if (rank == 0) {
      std::cerr << "Usage: " << argv[0]
                << " <dim> <n_points> <max_iter> <delta_x> [seed]\n";
    }
    MPI_Finalize();
    return 1;
  }

  unsigned int dim     = std::atoi(argv[1]);
  unsigned int n_points= std::atoi(argv[2]);
  unsigned int max_iter= std::atoi(argv[3]);
  double delta_x       = std::atof(argv[4]);
  unsigned int seed    = (argc > 5) ? static_cast<unsigned int>(std::stoul(argv[5])) : 12345;
  int m = 3;                                // for scale-free network, number of edges of each new node
  int iterations_stagnation = std::max(100, static_cast<int>(max_iter / 4)); // number of iterations for stagnation control
  double stagnation_tol = 1e-8;             // tolerance for stagnation-based stopping
  double stagnation_rel_tol = 1e-4;         // relative tolerance for stagnation-based stopping
  double diversity_tol = 1e-4;              // diversity tolerance for stopping criteria
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
      stagnation_rel_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory,
      seed
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
      stagnation_rel_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory,
      seed
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
      stagnation_rel_tol,
      diversity_tol,
      p_rewiring,
      p_random,
      m,
      function_names,
      factory,
      seed
  );

  ExperimentStats classic_stats = run_classic_experiment(
      rank,
      dim,
      n_points,
      max_iter,
      delta_x,
      iterations_stagnation,
      stagnation_tol,
      stagnation_rel_tol,
      diversity_tol,
      function_names,
      factory,
      seed
  );

  if (rank == 0) {
    
    std::array<std::vector<std::string>, 5> all = {
        small_stats.functions_converged,
        scale_stats.functions_converged,
        random_stats.functions_converged,
        classic_stats.functions_converged,
        function_names
    };
  
    /*
    print_experiment_stats("Small world", small_stats);
    print_experiment_stats("Scale free", scale_stats);
    print_experiment_stats("Random", random_stats);
    print_experiment_stats("Classic", classic_stats);
    */
    int n = not_converged(all, false);


  // output in human friendly format
    /*
    std::cout << "Total time classic PSO: " << classic_stats.total_time << " s\n";
    std::cout << "Convergence rate classic PSO: " << classic_stats.number_of_converged << "/" << number_of_functions  << std::endl;

    std::cout << "Total time small-world network timer version: " << small_stats.t_allgatherv << "/"  << small_stats.total_time << " s\n";
    std::cout << "Convergence rate small-world network: " << small_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time scale-free network timer version: " << scale_stats.t_allgatherv << "/"  << scale_stats.total_time << " s\n";
    std::cout << "Convergence rate scale-free network: " << scale_stats.number_of_converged << "/" << number_of_functions << std::endl;

    std::cout << "Total time random network timer version: " << random_stats.t_allgatherv << "/"  << random_stats.total_time << " s\n";
    std::cout << "Convergence rate random network: " << random_stats.number_of_converged << "/" << number_of_functions << std::endl << std::endl;
    */

    //uniform output format for benchmarking
    //std::cout << "\n";
    
    // Intestazione:
    // std::cout << "RESULT,method,dim,n_points,max_iter,delta_x,seed,time_total,time_comm,converged,total\n";

    std::cout << "RESULT,classic," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
              << classic_stats.total_time << "," << classic_stats.t_allgatherv << "," 
              << classic_stats.number_of_converged << "," << number_of_functions << "\n";

    std::cout << "RESULT,scale_free," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
              << scale_stats.total_time << "," << scale_stats.t_allgatherv << "," 
              << scale_stats.number_of_converged << "," << number_of_functions << "\n";

    std::cout << "RESULT,small_world," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
              << small_stats.total_time << "," << small_stats.t_allgatherv << "," 
              << small_stats.number_of_converged << "," << number_of_functions << "\n";

    std::cout << "RESULT,random," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
              << random_stats.total_time << "," << random_stats.t_allgatherv << "," 
              << random_stats.number_of_converged << "," << number_of_functions << "\n";

    std::cout << "RESULT,not_converged," << n << "," << number_of_functions << "\n";

    std::cout << "\n";
}
  MPI_Finalize();
  return 0;
}
