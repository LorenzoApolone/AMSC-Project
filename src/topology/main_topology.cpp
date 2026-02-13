#include "create_network.hpp"
#include "confront.hpp"
#include "pso_topology.hpp"
#include "../methods.hpp"
#include "../functions.cpp" 
    
#include <mpi.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// --- helper: broadcast di vector<vector<int>> ---
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

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Argouments: <dim> <n_points> <p> <max_iter> <delta_x>
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
  double p             = 0.05; // rewiring probability for small-world network
  int m = 3; // for scale-free network, number of edges of each new node
  int number_of_converged_small = 0;
  int number_of_converged_scale = 0;
  int number_of_converged_random = 0;
  int number_of_converged_classic = 0;
  int number_of_converged_complete = 0;
  int number_of_functions = 0;

  int number_of_converged_small_v1 = 0;
  int number_of_converged_scale_v1 = 0;
  int number_of_converged_random_v1 = 0;
  double p_rewiring = 0.05; // rewiring probability for small-world network
  double p_random = 0.03; // edge probability for random network
  std::vector<std::string> functions_converged_small;
  std::vector<std::string> functions_converged_scale;
  std::vector<std::string> functions_converged_random;
  std::vector<std::string> functions_converged_classic;
  std::vector<std::string> functions_converged_complete;
  StopCriterion stop(max_iter, delta_x);

  // Factory Definition 
  std::unordered_map<std::string,
                     std::function<std::unique_ptr<TestFunction>(unsigned int)>> factory;

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
  /* funzioni che convergono difficilmente 
  std::vector<std::string> function_names = {
   "DropWave","Alpine1",
    "HGBat","ExpSchafferF6","Schwefel",
    "XinSheYang1", "Salomon",
    "ModifiedXinSheYang3","ModifiedXinSheYang5"
  };
*/
  
  std::vector<std::string> function_names = {
    "Sphere","Ellipsoid","SumOfDiffPowers","DropWave","Weierstrass","Alpine1","Ackley",
    "Griewank","Rastrigin","HappyCat","HGBat","Rosenbrock","HighCondElliptic","Discus",
    "BentCigar","PermdbFunc","Schafferf7Func","ExpSchafferF6","RotatedHyper","Schwefel",
    "SumOfDifferentPowers2","XinSheYang1","Schwefel221","Schwefel222","Salomon",
    "ModifiedRidge","Zakharov","ModifiedXinSheYang3","ModifiedXinSheYang5", 
    "Levy","Michalewicz","Bohachevsky","Powell","DixonPrice","StyblinskiTang"
  };

  
  //+++++++++++++++++++++++++++++++++++++++++++++++++Timer version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  timer MPI evaluate to delete it on the final version
  MPI_Barrier(MPI_COMM_WORLD);
  double t_start_small = MPI_Wtime();
  double t_allgatherv_small = 0.0;
  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr = factory[name](dim);
    std::vector<std::vector<int>> adjacency_list;
  
    if (rank == 0) {
     // create_small_world_network(static_cast<int>(n_points), m, adjacency_list);
      create_network(static_cast<int>(n_points), p_rewiring, adjacency_list);
      number_of_functions++;
    }
    bcast_adjacency_list(adjacency_list, static_cast<int>(n_points), rank);
    OutputObject result = pso_small_timerv(*f_ptr, dim, stop, n_points, adjacency_list, converged, t_allgatherv_small);
    if (rank == 0) {
  //    result.terminal_info();
  //  result.output_to_file();     
      if(converged == true){
        functions_converged_small.push_back(name);
        number_of_converged_small++;
      }

    }
      
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_small = MPI_Wtime();

//+++++++++++++++++++++++++++++++++++++++++++++++Timer version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*
  //+++++++++++++++++++++++++++++++++++++++++++++++++Normal version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  timer MPI evaluate to delete it on the final version
  MPI_Barrier(MPI_COMM_WORLD);
  double t_start_small_v1 = MPI_Wtime();
  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr = factory[name](dim);
    std::vector<std::vector<int>> adjacency_list_v1;
  
    if (rank == 0) {
     // create_small_world_network(static_cast<int>(n_points), m, adjacency_list);
      create_network(static_cast<int>(n_points), p_rewiring, adjacency_list_v1);
      
    }
    bcast_adjacency_list(adjacency_list_v1, static_cast<int>(n_points), rank);
    OutputObject result = pso_normal(*f_ptr, dim, stop, n_points, adjacency_list_v1, converged);
    if (rank == 0) {
  //    result.terminal_info();
  //  result.output_to_file();     
      if(converged == true){
        number_of_converged_small_v1++;
      }

    }
      
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_small_v1 = MPI_Wtime();

//+++++++++++++++++++++++++++++++++++++++++++++++Normal version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/


//+++++++++++++++++++++++++++++++++++++++++++++++++Timer version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MPI_Barrier(MPI_COMM_WORLD);
  double t_start_scale = MPI_Wtime();
  double t_allgatherv_scale = 0.0;
  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr1 = factory[name](dim);
    std::vector<std::vector<int>> adjacency_list2;
  
    if (rank == 0) {
      create_scale_free_network(static_cast<int>(n_points), m, adjacency_list2);
    }
    bcast_adjacency_list(adjacency_list2, static_cast<int>(n_points), rank);
    OutputObject result = pso_small_timerv(*f_ptr1, dim, stop, n_points, adjacency_list2, converged, t_allgatherv_scale);
    if (rank == 0) {
 //   result.terminal_info();
 //   result.output_to_file();
      if(converged == true){
        number_of_converged_scale++;
        functions_converged_scale.push_back(name);
      }
    }    
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_scale = MPI_Wtime();

//++++++++++++++++++++++++++++++++++++++++++++++Timer Version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*
//+++++++++++++++++++++++++++++++++++++++++++++++++Normal version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MPI_Barrier(MPI_COMM_WORLD);
  double t_start_scale_v1 = MPI_Wtime();
  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr1 = factory[name](dim);
    std::vector<std::vector<int>> adjacency_list_v2;
  
    if (rank == 0) {
      create_scale_free_network(static_cast<int>(n_points), m, adjacency_list_v2);
    }
    bcast_adjacency_list(adjacency_list_v2, static_cast<int>(n_points), rank);
    OutputObject result = pso_normal(*f_ptr1, dim, stop, n_points, adjacency_list_v2, converged);
    if (rank == 0) {
 //   result.terminal_info();
 //   result.output_to_file();
      if(converged == true){
        number_of_converged_scale_v1++;
      }
    }    
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_scale_v1 = MPI_Wtime();

//++++++++++++++++++++++++++++++++++++++++++++++Normal Version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*/
//+++++++++++++++++++++++++++++++++++++++++++++++Timer Version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
MPI_Barrier(MPI_COMM_WORLD);
double t_start_random = MPI_Wtime();
double t_allgatherv_random = 0.0;

for (const auto& name : function_names) {
  bool converged = false;
  auto f_ptr2 = factory[name](dim);
  std::vector<std::vector<int>> adjacency_list2;

  if (rank == 0) {
    create_random_network(static_cast<int>(n_points), p_random, adjacency_list2);
  }

  bcast_adjacency_list(adjacency_list2, static_cast<int>(n_points), rank);
  OutputObject result = pso_small_timerv(*f_ptr2, dim, stop, n_points, adjacency_list2, converged, t_allgatherv_random);

  if (rank == 0 && converged) {
    number_of_converged_random++;
    functions_converged_random.push_back(name);
  }
}

MPI_Barrier(MPI_COMM_WORLD);
double t_end_random = MPI_Wtime();


//+++++++++++++++++++++++++++++++++++++++++++++++Timer Version+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/*
//+++++++++++++++++++++++++++++++++++++++++++++++++Normal version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 
MPI_Barrier(MPI_COMM_WORLD);
double t_start_random_v1 = MPI_Wtime();

for (const auto& name : function_names) {
  bool converged = false;
  auto f_ptr2 = factory[name](dim);
  std::vector<std::vector<int>> adjacency_list_v2;

  if (rank == 0) {
    create_random_network(static_cast<int>(n_points), p_random, adjacency_list_v2);
  }

  bcast_adjacency_list(adjacency_list_v2, static_cast<int>(n_points), rank);
  OutputObject result = pso_normal(*f_ptr2, dim, stop, n_points, adjacency_list_v2, converged);

  if (rank == 0 && converged) {
    number_of_converged_random_v1++;
  }
}

MPI_Barrier(MPI_COMM_WORLD);
double t_end_random_v1 = MPI_Wtime();

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++Normal version++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
 MPI_Barrier(MPI_COMM_WORLD);


MPI_Barrier(MPI_COMM_WORLD);
  double t_start_classic = MPI_Wtime();

  for (const auto &name : function_names)
  {
    bool converged = false;
    auto f_ptr = factory[name](dim);
    OutputObject result = pso_mpi(*f_ptr, dim, stop, n_points, converged);
    if (rank == 0)
    {
 //     result.terminal_info();
 //     result.output_to_file();
      if (converged){
        number_of_converged_classic++;
        functions_converged_classic.push_back(name);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_classic = MPI_Wtime();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*

  MPI_Barrier(MPI_COMM_WORLD);
  double t_start_complete = MPI_Wtime();

  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr3 = factory[name](dim);
    std::vector<std::vector<int>> adjacency_list3;
  
    if (rank == 0) {
      create_fully_connected_network(static_cast<int>(n_points), adjacency_list3);
    }
    bcast_adjacency_list(adjacency_list3, static_cast<int>(n_points), rank);
    OutputObject result = pso_normal(*f_ptr3, dim, stop, n_points, adjacency_list3, converged);
    if (rank == 0) {
 //     result.terminal_info();
  //  result.output_to_file();
      if(converged == true){
        number_of_converged_complete++;
        functions_converged_complete.push_back(name);
      }
    }
      
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_complete = MPI_Wtime();
*/
//--------------------------------------------Codice di prova---------------------------------------------------------------------------

/*
 MPI_Barrier(MPI_COMM_WORLD);


  double t_start_complete = MPI_Wtime();
  std::vector<std::vector<int>> adjacency_list4;
  if (rank == 0) {
      create_fully_connected_network(static_cast<int>(n_points), adjacency_list4);
    }
    bcast_adjacency_list(adjacency_list4, static_cast<int>(n_points), rank);
  for (const auto& name : function_names) {
    bool converged = false;
    auto f_ptr4 = factory[name](dim);
   
  
    
    OutputObject result = pso_small_debugger(*f_ptr4, dim, stop, n_points, adjacency_list4, converged);
    if (rank == 0) {
 //     result.terminal_info();
  //  result.output_to_file();
      if(converged == true){
        number_of_converged_complete++;
        functions_converged_complete.push_back(name);
      }
    }
      
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double t_end_complete = MPI_Wtime();

  */
//--------------------------------------------Fine codice di prova---------------------------------------------------------------------------
if (rank == 0) {

    std::cout << "Total time classic PSO: " << (t_end_classic - t_start_classic) << " s\n";
    std::cout << "Convergence rate classic PSO: " << number_of_converged_classic << "/" << number_of_functions << std::endl << std::endl;

    std::cout << "Total time scale-free network timer version: " << t_allgatherv_scale << "/" << (t_end_scale - t_start_scale) << " s\n";
    std::cout << "Convergence rate scale-free network: " << number_of_converged_scale << "/" << number_of_functions << std::endl;
 //   std::cout << "Total time scale-free network normal version: " <<  (t_end_scale_v1 - t_start_scale_v1) << " s\n \n";


    std::cout << "Total time small-world network timer version: "  << t_allgatherv_small << "/" << (t_end_small - t_start_small) << " s\n";
    std::cout << "Convergence rate small-world network: " << number_of_converged_small << "/" << number_of_functions << std::endl;
 //   std::cout << "Total time small-world network normal version: "  << (t_end_small_v1 - t_start_small_v1) << " s\n \n";

    std::cout << "Total time random network timer version: " << t_allgatherv_random << "/" << (t_end_random - t_start_random) << " s\n\n";
    std::cout << "Convergence rate random network: " << number_of_converged_random << "/" << number_of_functions << std::endl << std::endl;
 //   std::cout << "Total time random network normal version: " << (t_end_random_v1 - t_start_random_v1) << " s\n\n";
    /*
    std::cout << "Total time complete network: " << (t_end_complete - t_start_complete) << " s\n";
    std::cout << "Convergence rate complete network: " << number_of_converged_complete << "/" << number_of_functions << std::endl << std::endl;
  
  */
    std::array<std::vector<std::string>, 5> all = {
        functions_converged_small,
        functions_converged_scale,
        functions_converged_random,
        functions_converged_classic,
        function_names
    };
    
    not_converged(all, 1);

  
 
    
//    uniqueness(all);
  }
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  MPI_Finalize();
  return 0;
}
