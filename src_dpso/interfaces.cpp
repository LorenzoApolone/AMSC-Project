/**
 * @file output_object.cpp
 * @brief Implementation of OutputObject methods for PSO benchmarking framework.
 *
 * This file defines functions to print results to the terminal and
 * to export them into structured output files under the "tests" directory.
 *
 */

#include "interfaces.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>

using namespace std;
/**
 * @brief Print summary information about a completed optimization run to the terminal.
 *
 * Displays key metrics such as:
 * - Test function name
 * - Problem dimension
 * - Number of particles
 * - Number of cores used
 * - Final convergence value (last Δx)
 * - Execution time
 * - Number of iterations performed
 *
 * Call it after a PSO run to get a quick overview of results. not suggested for massive testing, in this case use output_to_file() 
 */
void OutputObject::terminal_info(){
    std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Test function: " << function_name << std::endl;
    std::cout << "Dimension: " << x_best.size() << std::endl;
    std::cout << "Number of points: " << n_points << std::endl;
    std::cout << "Number of cores: " << n_cores << std::endl;
    std::cout << "Stopping criterion: " << std::endl;
    std::cout << "Final delta x: " << conv_history[conv_history.size() - 1] << std::endl;
    std::cout << "Execution time: " << execution_time << std::endl;
    std::cout << "Iterations: " << conv_history.size() << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << std::endl;
}

int get_max_test_number(const std::filesystem::path& dir) {
    int max_num = -1;
    const std::string prefix = "test_";
    const std::string suffix = ".txt";

    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;

        auto name = entry.path().filename().string();

        // Must start with "test_" and end with ".txt"
        if (name.size() <= prefix.size() + suffix.size()) continue;
        if (name.rfind(prefix, 0) != 0) continue;                               // does not start with prefix
        if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;                                                          // does not end with suffix

        // Extract numeric part
        std::string num_str = name.substr(prefix.size(),
                                          name.size() - prefix.size() - suffix.size());
        try {
            int n = std::stoi(num_str);
            if (n > max_num) max_num = n;
        } catch (...) {
            // Ignore malformed numbers
        }
    }

    return max_num;
}
/**
 * @brief Write PSO run results to a structured output file.
 *
 * Creates a nested folder structure under the `tests/` directory according to:
 * ```
 * tests/<function_name>/<dimension>/<n_points>/<n_cores>/
 * ```
 * Each file is named sequentially as `test_<N>.txt`, where `<N>` increments automatically.
 *
 * The output file contains one line per iteration, with the following columns:
 * ```
 * max_iter  tol  it_n  delta_x  final_t
 * ```
 *
 */
void OutputObject::output_to_file(){
    // Creating and open a text file (and folders, if needed)
    std::filesystem::path save_dir =
        std::filesystem::path("tests") /
        function_name /
        std::to_string(x_best.size()) /
        std::to_string(n_points) /
        std::to_string(n_cores);
        
    std::filesystem::create_directories(save_dir);

    int file_n = get_max_test_number(save_dir) + 1;     // the files are named test_0, test_1, test_2 and so ons
    std::string filename = "test_" + std::to_string(file_n) + ".txt";
    ofstream MyFile(save_dir / filename);

    // Write to the file
    MyFile << "max_iter tol it_n delta_x final_t\n";

    for (size_t i = 0; i < conv_history.size(); i++){
        MyFile << stopcriterion.get_max_iter() << " "
        << stopcriterion.get_tolerance() << " "
        << i << " "
        << conv_history[i] << " "
        << execution_time << "\n";
    }

    // Close the file
    MyFile.close();
};