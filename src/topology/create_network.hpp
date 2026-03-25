/**
 * @file create_network.hpp
 * @brief Network topology generators for PSO communication graphs.
 *
 * This module provides functions that build an adjacency list representing
 * an undirected graph of N nodes, where each node stores the indices of its
 * neighbors.
 *
 * @note The adjacency list uses 0-based node indices in range [0, N-1].
 */
#pragma once
#include <vector>

/**
 * @enum NetworkType
 * @brief Supported network topology types.
 */
enum class NetworkType {
  SMALL_WORLD,      /**< Watts–Strogatz-like small-world (ring lattice + rewiring). */
  SCALE_FREE,       /**< Barabási–Albert-like preferential attachment. */
  RANDOM,           /**< Erdős–Rényi G(N,p) random graph. */
};

/**
 * @brief Generate a small-world network (ring lattice with rewiring).
 *
 * Builds a ring lattice with degree K=2 (each node connected to its two
 * immediate neighbors) and then rewires each existing edge endpoint with
 * probability @p p to a uniformly chosen new node that is not equal to the
 * source and is not already a neighbor of the source.
 * @param[in]  N              Number of nodes (must be > 0).
 * @param[in]  p              Rewiring probability in [0, 1].
 * 
 * @param[out] adjacency_list Adjacency list of size N containing the resulting graph.
 *
 * @pre N > 0
 * @pre 0.0 <= p <= 1.0
 *
 */
void create_network(int N, double p, std::vector<std::vector<int>> &adjacency_list, unsigned int seed);

/**
 * @brief Generate a scale-free network via preferential attachment.
 *
 * Uses a Barabási–Albert-style construction:
 * - Start from a fully connected seed of size m0 = m + 1.
 * - Iteratively add new nodes; each new node attaches to m distinct existing
 *   nodes, sampled with probability proportional to current degree
 *
 * @param[in]  N              Number of nodes.
 * @param[in]  m              Number of edges to attach from each new node
 *                            
 * @param[out] adjacency_list Adjacency list of size N containing the resulting graph.
 *
 * @pre N >= 0
 * @pre m >= 1
 * @pre m < N 
 */
void create_scale_free_network(int N, int m, std::vector<std::vector<int>> &adjacency_list, unsigned int seed);

/**
 * @brief Generate an Erdős–Rényi random graph G(N, p).
 *
 * For each unordered pair (i, j) with i < j, an edge is added with probability
 * @p p. The adjacency list is symmetric.
 *
 * @param[in]  N              Number of nodes.
 * @param[in]  p              Edge probability in [0, 1]. 
 * 
 * @param[out] adjacency_list Adjacency list of size N.
 *
 * @pre N >= 0
 * @pre 0.0 <= p <= 1.0
 *
 */
void create_random_network(int N, double p, std::vector<std::vector<int>> &adjacency_list, unsigned int seed);
