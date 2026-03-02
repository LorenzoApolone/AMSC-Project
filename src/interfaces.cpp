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
        MyFile << stopcriterion.get_max_iters() << " "
        
        << i << " "
        << conv_history[i] << " "
        << execution_time << "\n";
     }

    // Close the file
    MyFile.close();
};


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rifai, hai tolto .get_tollerance

NON USARE ASSOLUTAMENTE PRIMA DI AVER CHIESTO A LORENZO 

questa funzione salva su file .csv i risultati della soluzione di pso per ogni funzione testata. 
I risultati sono salvati in un file .csv con nome che identifica il metodo usato, la dimensione, il numero di particelle, il numero di core e il numero di ripetizione del test.
Se tutti questi parametri rimangono uguali in due esecuzioni diverse, i risultati vengono appesi in fondo nello stesso file. NON VA BENE. 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void OutputObject::append_summary_csv_by_method(const std::string& method_name, int rep) const
{
    namespace fs = std::filesystem;

    fs::path out_dir = fs::path("tests") / "summary";
    fs::create_directories(out_dir);

    const std::string m = sanitize_filename(method_name);

    fs::path csv_path =
        out_dir / (m +
                   "_d"  + std::to_string(d) +
                   "_np" + std::to_string(n_points) +
                   "_c"  + std::to_string(n_cores) +
                   "_rep" + std::to_string(rep) +
                   ".csv");

    const bool file_exists = fs::exists(csv_path);

    std::ofstream out(csv_path, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Errore apertura CSV: " << csv_path << "\n";
        return;
    }

    if (!file_exists) {
        out << "# method=" << method_name << "\n";
        out << "# dim=" << d << "\n";
        out << "# n_points=" << n_points << "\n";
        out << "# n_cores=" << n_cores << "\n";
        out << "# max_iter=" << stopcriterion.get_max_iters() << "\n";
        out << "# rep=" << rep << "\n\n";

        out << "function,converged,iters,final_err,final_fitness,time_s\n";
    }

    const size_t iters = conv_history.size();
    const double final_err = (iters > 0 ? conv_history.back()
                                        : std::numeric_limits<double>::quiet_NaN());

  
    

    out << function_name << ","
        
        << iters << ","
        << final_err << ","
        << f_val << ","
        << execution_time
        << "\n";
}