// create_network.cpp
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <set>
#include "create_network.hpp"

void create_network(int N, double p, std::vector<std::vector<int>> &adjacency_list) {
    adjacency_list.assign(N, std::vector<int>());
    std::set<std::pair<int, int>> edges;

    for (int i = 0; i < N; ++i) {
        int r = (i + 1) % N;
        edges.insert({std::min(i, r), std::max(i, r)});
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> node_dis(0, N - 1);

    for (auto it = edges.begin(); it != edges.end(); ) {
        if (dis(gen) < p) {
            int u = it->first;
            int v = it->second;
            int new_v;

            do {
                new_v = node_dis(gen);
            } while (new_v == u || new_v == v ||
                     edges.count({std::min(u, new_v), std::max(u, new_v)}));

            it = edges.erase(it);
            edges.insert({std::min(u, new_v), std::max(u, new_v)});
        } else {
            ++it;
        }
    }

    for (const auto &e : edges) {
        adjacency_list[e.first].push_back(e.second);
        adjacency_list[e.second].push_back(e.first);
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

    for (int i = 0; i < m0; ++i) {
        for (int k = 0; k < degree[i]; ++k) {
            pool.push_back(i);
        }
    }

    for (int new_node = m0; new_node < N; ++new_node) {
        std::vector<int> targets;
        targets.reserve(m);

        if (pool.empty()) {
            for (int t = 0; t < new_node && static_cast<int>(targets.size()) < m; ++t) {
                targets.push_back(t);
            }
        } else {
            std::uniform_int_distribution<> pick(0, static_cast<int>(pool.size()) - 1);
            while (static_cast<int>(targets.size()) < m) {
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

        for (int k = 0; k < degree[new_node]; ++k) {
            pool.push_back(new_node);
        }
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