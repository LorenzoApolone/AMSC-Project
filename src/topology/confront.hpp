// confront.hpp
#pragma once
#include <cstddef>
#include <vector>
#include <string>
#include <array>
/**
 * @file confront.hpp
 * @brief Utilities to compare convergence results across PSO topologies.
 *
 * This module provides helper functions that, given the list of benchmark
 * functions that converged under each topology, prints:
 * - functions that converged only in one specific topology
 * - functions that did not converge in any topology
 *
 * The input is an array of 5 vectors:
 * - index Topology::Small:        function names that converged in small-world PSO
 * - index Topology::Scale:        function names that converged in scale-free PSO
 * - index Topology::Random:       function names that converged in random-network PSO
 * - index Topology::Classic:      function names that converged in classic (global/neighborhood) PSO
 * - index Topology::FunctionsNames: full list of benchmark function names
 */

/**
 * @enum Topology
 * @brief Indices used to access topology-specific vectors in the input array.
 *
 * The enum values are intended to be used as indices 
 * 
 */


enum class Topology : std::size_t {
    Small  = 0,
    Scale  = 1,
    Random = 2,
    Classic = 3, 
    FunctionsNames = 4
};
/**
 * @brief Print functions that converge only in a single topology, and functions that never converge.
 *
 * Given a set of vectors containing function names that converged for each topology,
 * this function prints:
 * - functions that appear only in one topology vector (exclusive convergence),
 * - functions that appear in none of the topology vectors (non convergence).
 *
 * @param vectors Array containing:
 *   - vectors[Topology::Small]        converged in small-world
 *   - vectors[Topology::Scale]        converged in scale-free
 *   - vectors[Topology::Random]       converged in random
 *   - vectors[Topology::Classic]      converged in classic
 *   - vectors[Topology::FunctionsNames] all benchmark names
 *
 * @note The function uses linear searches do not use it for large datasets
 *       
 *
 * 
 */
void uniqueness(const std::array<std::vector<std::string>, 5>& vectors);

/**
 * @brief Print the name of functions that did not converge in any topology.
 *
 * Iterates over `vectors[Topology::FunctionsNames]` and checks whether each name is absent
 * from all topology-specific vectors. Every such function is printed and counted.
 *  
 *
 * @param vectors Same structure as in uniqueness().
 * @return Number of benchmark functions that did not converge in any topology.
 *
 * @warning The function performs output on `std::cout` (side effects).
 */
int not_converged(const std::array<std::vector<std::string>, 5>& vectors);