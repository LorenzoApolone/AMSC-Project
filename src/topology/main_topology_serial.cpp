/**
 * @file main_topology_serial.cpp
 * @brief Main file for serial topology-based PSO experiments.
 *
 * This file executes and compares multiple serial PSO variants 
 * (Small-World, Scale-Free, Random) 
 * over a test suite of benchmark functions available in the "functions.hpp" file.
 * All the PSO variants are executed sequentially (in serial). 
 *
 * All setups, including graph generation, and statistics collections, 
 * are taken from the "benchmark_utils.hpp" module.
 * 
 * The program:
 * 1. initializes MPI (for internal benchmark utils timing consistency),
 * 2. parses command-line parameters,
 * 3. imports the function factory from benchmark utilities,
 * 4. runs the selected serial PSO variants,
 * 5. collects convergence and timing metrics,
 * 6. prints machine-readable summaries.
 */

#include "benchmark_utils.hpp"

/**
 * @brief Main of the serial benchmark program.
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

  if (rank != 0) {
      MPI_Finalize();
      return 0;
  }

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

  if (dim == 0 || n_points == 0 || max_iter == 0) {
    std::cerr << "Error: <dim>, <n_points> and <max_iter> must be strictly greater than 0!\n";
    return 1;
  }
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

  ExperimentStats serial_small_stats = run_serial_topology_experiment(
      TopologyMode::SMALL_WORLD,
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

  ExperimentStats serial_scale_stats = run_serial_topology_experiment(
      TopologyMode::SCALE_FREE,
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

  ExperimentStats serial_random_stats = run_serial_topology_experiment(
      TopologyMode::RANDOM,
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

  ExperimentStats serial_classic_stats = run_serial_classic_experiment(
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

  std::array<std::vector<std::string>, 5> all = {
      serial_small_stats.functions_converged,
      serial_scale_stats.functions_converged,
      serial_random_stats.functions_converged,
      serial_classic_stats.functions_converged,
      function_names
  };
  
  int n = not_converged(all, false);

  std::cout << "RESULT,classic," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
            << serial_classic_stats.total_time << ",0," 
            << serial_classic_stats.number_of_converged << "," << number_of_functions << "\n";

  std::cout << "RESULT,scale_free," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
            << serial_scale_stats.total_time << ",0," 
            << serial_scale_stats.number_of_converged << "," << number_of_functions << "\n";

  std::cout << "RESULT,small_world," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
            << serial_small_stats.total_time << ",0," 
            << serial_small_stats.number_of_converged << "," << number_of_functions << "\n";

  std::cout << "RESULT,random," << dim << "," << n_points << "," << max_iter << "," << delta_x << "," << seed << "," 
            << serial_random_stats.total_time << ",0," 
            << serial_random_stats.number_of_converged << "," << number_of_functions << "\n";

  std::cout << "RESULT,not_converged," << n << "," << number_of_functions << "\n";

  std::cout << "\n";

  MPI_Finalize();
  return 0;
}
