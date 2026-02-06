#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>

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
