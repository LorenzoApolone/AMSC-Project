/**
 * @file create_network.hpp
 * @brief Network topology generators for PSO communication graphs.
 *
 * This module provides functions that build an adjacency list representing
 * an undirected graph of N nodes, where each node stores the indices of its
 * neighbors.
 *
 * @note The adjacency list uses 0-based node indices in range [0, N-1].
 * @note `NetworkType` is declared in `NetworkType.hpp` and re-used here by
 *       higher-level topology factories.
 */
#pragma once

#include "NetworkType.hpp"
#include <random>
#include <vector>

/**
 * @brief Generate a small-world network using an internal random engine.
 *
 * Builds a ring lattice with degree 2 and rewires each existing edge with
 * probability @p p to a valid non-neighbor target.
 *
 * @param[in] N Number of nodes.
 * @param[in] p Rewiring probability in the range [0, 1].
 * @param[out] adjacency_list Output adjacency list of size N.
 */
void create_network(int N, double p,
                    std::vector<std::vector<int>> &adjacency_list);

/**
 * @brief Generate a small-world network using a caller-supplied random engine.
 *
 * This overload exposes the RNG so callers can obtain deterministic topologies
 * from a controlled seed.
 *
 * @param[in] N Number of nodes.
 * @param[in] p Rewiring probability in the range [0, 1].
 * @param[out] adjacency_list Output adjacency list of size N.
 * @param[in,out] gen Random engine used for rewiring decisions.
 */
void create_network(int N, double p,
                    std::vector<std::vector<int>> &adjacency_list,
                    std::mt19937 &gen);

/**
 * @brief Generate a scale-free network via preferential attachment using an
 * internal random engine.
 *
 * The construction starts from a fully connected seed of size m + 1 and then
 * adds one node at a time, attaching it to m distinct existing nodes sampled
 * proportionally to degree.
 *
 * @param[in] N Number of nodes.
 * @param[in] m Number of links added by each new node.
 * @param[out] adjacency_list Output adjacency list of size N.
 */
void create_scale_free_network(int N, int m,
                               std::vector<std::vector<int>> &adjacency_list);

/**
 * @brief Generate a scale-free network via preferential attachment using a
 * caller-supplied random engine.
 *
 * @param[in] N Number of nodes.
 * @param[in] m Number of links added by each new node.
 * @param[out] adjacency_list Output adjacency list of size N.
 * @param[in,out] gen Random engine used for preferential sampling.
 */
void create_scale_free_network(int N, int m,
                               std::vector<std::vector<int>> &adjacency_list,
                               std::mt19937 &gen);

/**
 * @brief Generate an Erdos-Renyi random graph G(N, p) using an internal
 * random engine.
 *
 * For each unordered pair of distinct nodes (i, j), an undirected edge is
 * inserted with probability @p p.
 *
 * @param[in] N Number of nodes.
 * @param[in] p Edge probability in the range [0, 1].
 * @param[out] adjacency_list Output adjacency list of size N.
 */
void create_random_network(int N, double p,
                           std::vector<std::vector<int>> &adjacency_list);

/**
 * @brief Generate an Erdos-Renyi random graph G(N, p) using a caller-supplied
 * random engine.
 *
 * @param[in] N Number of nodes.
 * @param[in] p Edge probability in the range [0, 1].
 * @param[out] adjacency_list Output adjacency list of size N.
 * @param[in,out] gen Random engine used for edge sampling.
 */
void create_random_network(int N, double p,
                           std::vector<std::vector<int>> &adjacency_list,
                           std::mt19937 &gen);