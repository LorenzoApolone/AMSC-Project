#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include "create_network.hpp"
using namespace std;

void create_network(int N, double p, vector<vector<int>> &adjacency_list)
{
    adjacency_list.clear();
    adjacency_list.resize(N);

    // Ring lattice (K = 2)
    for (int i = 0; i < N; ++i)
    {
        adjacency_list[i].push_back((i + 1) % N);
        adjacency_list[i].push_back((i - 1 + N) % N);
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1.0);
    uniform_int_distribution<> node_dis(0, N - 1);

    // Rewiring
    for (int i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < adjacency_list[i].size(); ++j)
        {
            if (dis(gen) < p)
            {
                int new_node;
                do
                {
                    new_node = node_dis(gen);
                }
                while (
                    new_node == i ||
                    find(adjacency_list[i].begin(),
                         adjacency_list[i].end(),
                         new_node) != adjacency_list[i].end()
                );

                adjacency_list[i][j] = new_node;
            }
        }
    }
 
}

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

void create_random_network(int N, double p, std::vector<std::vector<int>>& adjacency_list)
{
    adjacency_list.assign(N, {});
    if (N <= 0 || p < 0.0 || p > 1.0) return;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (dis(gen) < p) {
                adjacency_list[i].push_back(j);
                adjacency_list[j].push_back(i);
            }
        }
    }
}

void create_fully_connected_network(int N, std::vector<std::vector<int>>& adjacency_list)
{
    adjacency_list.assign(N, {});
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                adjacency_list[i].push_back(j);
            }
        }
    }
}