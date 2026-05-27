// confront.cpp
#include <vector>
#include <string>
#include <iostream>
#include "confront.hpp"
#include <array>
#include <algorithm>
#include <unordered_set>


/**
 * @file confront.cpp
 * @brief Implementation of utilities to compare convergence across PSO topologies.
 */

void uniqueness(const std::array<std::vector<std::string>, 5>& vectors){

    // enum declared in hpp file  carefull to use the correct index to access the vectors
    const auto& small  = vectors[static_cast<std::size_t>(Topology::Small)];
    const auto& scale  = vectors[static_cast<std::size_t>(Topology::Scale)];
    const auto& random = vectors[static_cast<std::size_t>(Topology::Random)];
    const auto& classic = vectors[static_cast<std::size_t>(Topology::Classic)];
    const auto& functions_names = vectors[static_cast<std::size_t>(Topology::FunctionsNames)];

    std::unordered_set<std::string> small_set(small.begin(), small.end());
    std::unordered_set<std::string> scale_set(scale.begin(), scale.end());
    std::unordered_set<std::string> random_set(random.begin(), random.end());
    std::unordered_set<std::string> classic_set(classic.begin(), classic.end());

    bool is_found_small = false;
    bool is_found_scale = false;
    bool is_found_random = false;
    bool is_found_classic = false;

    for (const auto& name : small) {
        if (scale_set.find(name) == scale_set.end() &&
            random_set.find(name) == random_set.end() &&
            classic_set.find(name) == classic_set.end()) {
            std::cout << "Function " << name << " converged only in small-world topology.\n";
            is_found_small = true;
        }
    }
    if (is_found_small == false) {
        std::cout << "No function converged only in small-world topology.\n";
    }   

    for (const auto& name : scale) {
        if (small_set.find(name) == small_set.end() &&
            random_set.find(name) == random_set.end() &&
            classic_set.find(name) == classic_set.end()) {
            std::cout << "Function " << name << " converged only in scale-free topology.\n";
            is_found_scale = true;
        }
    }
    if (is_found_scale == false) {
        std::cout << "No function converged only in scale-free topology.\n";
    } 
    for (const auto& name : random) {
        if (small_set.find(name) == small_set.end() &&
            scale_set.find(name) == scale_set.end() &&
            classic_set.find(name) == classic_set.end()) {
            std::cout << "Function " << name << " converged only in random topology.\n";
            is_found_random = true;
        }
    }
    if (is_found_random == false) {
        std::cout << "No function converged only in random topology.\n";
    }

    for (const auto& name : classic) {
        if (small_set.find(name) == small_set.end() &&
            scale_set.find(name) == scale_set.end() &&
            random_set.find(name) == random_set.end()) {
            std::cout << "Function " << name << " converged only in classic topology.\n";
            is_found_classic = true;
        }
    }
    if (is_found_classic == false) {
        std::cout << "No function converged only in classic topology.\n";
    }

}

int not_converged(const std::array<std::vector<std::string>, 5>& vectors, bool print){

    const auto& small  = vectors[static_cast<std::size_t>(Topology::Small)];
    const auto& scale  = vectors[static_cast<std::size_t>(Topology::Scale)];
    const auto& random = vectors[static_cast<std::size_t>(Topology::Random)];
    const auto& classic = vectors[static_cast<std::size_t>(Topology::Classic)];
    const auto& functions_names = vectors[static_cast<std::size_t>(Topology::FunctionsNames)];

    std::unordered_set<std::string> small_set(small.begin(), small.end());
    std::unordered_set<std::string> scale_set(scale.begin(), scale.end());
    std::unordered_set<std::string> random_set(random.begin(), random.end());
    std::unordered_set<std::string> classic_set(classic.begin(), classic.end());

    int count = 0;
    for (const auto& name : functions_names) {
        if (small_set.find(name) == small_set.end() &&
            scale_set.find(name) == scale_set.end() &&
            random_set.find(name) == random_set.end() &&
            classic_set.find(name) == classic_set.end()) {
            if (print) {
                std::cout << "Function " << name << " did not converge in any topology.\n";
            }
            count++;
        }
    }
    if (print) {
        std::cout << std::endl;
    }
    return count;
}
