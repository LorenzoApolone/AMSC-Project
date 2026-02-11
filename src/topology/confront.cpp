#include <vector>
#include <string>
#include <iostream>
#include "confront.hpp"
#include <array>
#include <algorithm>

// this function checks if there is functions that converges only in one topology. 
void uniqueness(const std::array<std::vector<std::string>, 5>& vectors){

    // enum declared in hpp file  carefull to use the correct index to access the vectors
    const auto& small  = vectors[static_cast<std::size_t>(Topology::Small)];
    const auto& scale  = vectors[static_cast<std::size_t>(Topology::Scale)];
    const auto& random = vectors[static_cast<std::size_t>(Topology::Random)];
    const auto& classic = vectors[static_cast<std::size_t>(Topology::Classic)];
    const auto& functions_names = vectors[static_cast<std::size_t>(Topology::FunctionsNames)];
    

    bool is_found_small = false;
    bool is_found_scale = false;
    bool is_found_random = false;
    bool is_found_classic = false;
    // 'find' search an element in a vector, if the element is found than != to vector.end() otherwise it is not found 
    for (const auto& name : small) {
        if (std::find(scale.begin(), scale.end(), name) == scale.end() &&
            std::find(random.begin(), random.end(), name) == random.end() &&
            std::find(classic.begin(), classic.end(), name) == classic.end()) {
            std::cout << "Function " << name << " converged only in small-world topology.\n";
            is_found_small = true;
        }
    }
    if (is_found_small == false) {
        std::cout << "No function converged only in small-world topology.\n";
    }   

    for (const auto& name : scale) {
        if (std::find(small.begin(), small.end(), name) == small.end() &&
            std::find(random.begin(), random.end(), name) == random.end() &&
            std::find(classic.begin(), classic.end(), name) == classic.end()) {
            std::cout << "Function " << name << " converged only in scale-free topology.\n";
            is_found_scale = true;
        }
    }
    if (is_found_scale == false) {
        std::cout << "No function converged only in scale-free topology.\n";
    } 
    for (const auto& name : random) {
        if (std::find(small.begin(), small.end(), name) == small.end() &&
            std::find(scale.begin(), scale.end(), name) == scale.end() &&
            std::find(classic.begin(), classic.end(), name) == classic.end()) {
            std::cout << "Function " << name << " converged only in random topology.\n";
            is_found_random = true;
        }
    }
    if (is_found_random == false) {
        std::cout << "No function converged only in random topology.\n";
    }

    for (const auto& name : classic) {
        if (std::find(small.begin(), small.end(), name) == small.end() &&
            std::find(scale.begin(), scale.end(), name) == scale.end() &&
            std::find(random.begin(), random.end(), name) == random.end()) {
            std::cout << "Function " << name << " converged only in classic topology.\n";
            is_found_classic = true;
        }
    }
    if (is_found_classic == false) {
        std::cout << "No function converged only in classic topology.\n";
    }

    for (const auto& name : functions_names) {
        if (std::find(small.begin(), small.end(), name) == small.end() &&
            std::find(scale.begin(), scale.end(), name) == scale.end() &&
            std::find(random.begin(), random.end(), name) == random.end() &&
            std::find(classic.begin(), classic.end(), name) == classic.end()) {
            std::cout << "Function " << name << " did not converge in any topology.\n";
        }
    }
}


void not_converged(const std::array<std::vector<std::string>, 5>& vectors, int v){

    const auto& small  = vectors[static_cast<std::size_t>(Topology::Small)];
    const auto& scale  = vectors[static_cast<std::size_t>(Topology::Scale)];
    const auto& random = vectors[static_cast<std::size_t>(Topology::Random)];
    const auto& classic = vectors[static_cast<std::size_t>(Topology::Classic)];
    const auto& functions_names = vectors[static_cast<std::size_t>(Topology::FunctionsNames)];
    std::cout << "Version " << v << std::endl;
    for (const auto& name : functions_names) {
        if (std::find(small.begin(), small.end(), name) == small.end() &&
            std::find(scale.begin(), scale.end(), name) == scale.end() &&
            std::find(random.begin(), random.end(), name) == random.end() &&
            std::find(classic.begin(), classic.end(), name) == classic.end()) {
            std::cout << "Function " << name << " did not converge in any topology.\n";
        }
    }
    std:: cout << std::endl;
}