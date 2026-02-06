#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include "create_network.hpp"
using namespace std;
// Scale-free network using Barabási–Albert model, N is the number of nodes, m is the number of edges of each new node
void create_scale_free_network(int N, int m, std::vector<std::vector<int>>& adjacency_list)
{
    adjacency_list.assign(N, {});
    if (N <= 0) return;

    if (m < 1) m = 1;
    if (m >= N) m = N - 1;

    std::random_device rd;
    std::mt19937 gen(rd());

    int m0 = m + 1;
    if (m0 > N) m0 = N;

    std::vector<int> degree(N, 0);
    std::vector<int> pool;

    for (int i = 0; i < m0; ++i) {
        for (int j = i + 1; j < m0; ++j) {
            adjacency_list[i].push_back(j);
            adjacency_list[j].push_back(i);
            degree[i]++;
            degree[j]++;
        }
    }

    for (int i = 0; i < m0; ++i)
        for (int k = 0; k < degree[i]; ++k)
            pool.push_back(i);

    for (int new_node = m0; new_node < N; ++new_node) {
        std::vector<int> targets;
        targets.reserve(m);

        if (pool.empty()) {
            for (int t = 0; t < new_node && (int)targets.size() < m; ++t)
                targets.push_back(t);
        } else {
            std::uniform_int_distribution<> pick(0, (int)pool.size() - 1);
            while ((int)targets.size() < m) {
                int cand = pool[pick(gen)];
                if (cand == new_node) continue;
                if (std::find(targets.begin(), targets.end(), cand) != targets.end()) continue;
                targets.push_back(cand);
            }
        }

        for (int t : targets) {
            adjacency_list[new_node].push_back(t);
            adjacency_list[t].push_back(new_node);
            degree[new_node]++;
            degree[t]++;
            pool.push_back(t);
        }

        for (int k = 0; k < degree[new_node]; ++k)
            pool.push_back(new_node);
    }
}