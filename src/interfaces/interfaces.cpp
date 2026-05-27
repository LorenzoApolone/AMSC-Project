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
 * Call it after a PSO run to get a quick 
 * overview of results. not suggested for massive testing, in this case use output_to_file() 
 */

static std::string sanitize_filename(std::string s)
{
    for (char& c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c=='-' || c=='_'))
            c = '_';
    }
    return s;
}
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

