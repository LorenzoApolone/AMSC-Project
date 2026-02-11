#pragma once
#include <cstddef>
#include <vector>
#include <string>
#include <array>
enum class Topology : std::size_t {
    Small  = 0,
    Scale  = 1,
    Random = 2,
    Classic = 3, 
    FunctionsNames = 4
};

void uniqueness(const std::array<std::vector<std::string>, 5>& vectors);
void not_converged(const std::array<std::vector<std::string>, 5>& vectors, int v);