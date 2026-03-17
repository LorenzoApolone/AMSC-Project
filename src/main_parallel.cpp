#include "functions.hpp"
#include "interfaces/StoppingCriteriaManager.hpp"
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
  int stopped_by_maxiter_and_incorrect = 0;   
  int stopped_by_maxiter_and_correct = 0;   
  int incorrect_when_early_stop = 0;  
  int correct_when_early_stop = 0;
  int correct_total = 0;
  int number_of_functions = 0;
  int iterations_stagnation = 800;
  double stagnation_tol = 1e-9; // delta x for
  double diversity_tol = 1e-7; // diversity tolerance for stopping criteria
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
  factory["Levy"] = [](unsigned int dim){ return std::make_unique<Levy>(dim); };
  factory["Michalewicz"] = [](unsigned int dim){ return std::make_unique<Michalewicz>(dim); };
  factory["Bohachevsky"] = [](unsigned int dim){ return std::make_unique<Bohachevsky>(dim); };
  factory["Powell"] = [](unsigned int dim){ return std::make_unique<Powell>(dim); };
  factory["DixonPrice"] = [](unsigned int dim){ return std::make_unique<DixonPrice>(dim); };
  factory["StyblinskiTang"] = [](unsigned int dim){ return std::make_unique<StyblinskiTang>(dim); };
  factory["QuinticFunction"] = [](unsigned int dim){ return std::make_unique<QuinticFunction>(dim); };

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
                                             "ModifiedXinSheYang5",
                                             "Levy",
                                             "Michalewicz",
                                             "Bohachevsky",
                                             "Powell",
                                             "DixonPrice",
                                             "StyblinskiTang", 
                                             "QuinticFunction"};


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //--------------------------------------timer MPI evaluate to delete it on the final version----------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  double t_start = MPI_Wtime();
  //-------------------------------------------end of the timer-------------------------------------------------------
  // Run the solver
  for (const auto &name : function_names)
  {
    number_of_functions++;
    
    auto f_ptr = factory[name](dim);
    StoppingCriteriaManager stop(max_iter, iterations_stagnation, stagnation_tol, diversity_tol);

    OutputObject result = pso_mpi(*f_ptr, dim, stop, n_points);
    if (rank == 0)
    {
      const double final_fitness = result.get_best_fitness();
      const double f_star = f_ptr->value(f_ptr->get_true_solution());

      const bool is_correct = (std::abs(final_fitness - f_star) <= delta_x);
      if (is_correct) {
          correct_total++;
      }

      const bool bool_stopped_by_maxiter = (result.iterations >= (int)max_iter);

      if (bool_stopped_by_maxiter && !is_correct) {
          stopped_by_maxiter_and_incorrect++;   
      }
      if (bool_stopped_by_maxiter && is_correct) {
          stopped_by_maxiter_and_correct++;   
      }
      if (!is_correct && !bool_stopped_by_maxiter) {
          incorrect_when_early_stop++;
      }
      if (is_correct && !bool_stopped_by_maxiter) {
          correct_when_early_stop++;
      }
  //    result.terminal_info();
  //    result.output_to_file();
      
    }
  }

  //------------------------------------------- new part added to confront the time-------------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  double t_end = MPI_Wtime();

  if (rank == 0) {
    std::cout << "Number of functions: " << number_of_functions << std::endl;
    std::cout << "Total time: " << (t_end - t_start) << " s\n";
    std::cout << "Stopped by max iter and incorrect: " << stopped_by_maxiter_and_incorrect << std::endl;
    std::cout << "Stopped by max iter and correct: " << stopped_by_maxiter_and_correct << std::endl;
    std::cout << "Incorrect when early stop: " << incorrect_when_early_stop << std::endl;
    std::cout << "Correct when early stop: " << correct_when_early_stop << std::endl;
    std::cout << "Correct total: " << correct_total << std::endl;
  }
  //------------------------------------------- end of the new part-------------------------------------------------------
  
  MPI_Finalize();
  return 0;
}
