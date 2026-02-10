#pragma once
#include <vector>

void create_network(int N, double p, std::vector<std::vector<int>>& adjacency_list);
void create_scale_free_network(int N, int m, std::vector<std::vector<int>>& adjacency_list);
void create_random_network(int N, double p, std::vector<std::vector<int>>& adjacency_list);
void create_fully_connected_network(int N, std::vector<std::vector<int>>& adjacency_list);