/**
 * @file SubSwarmOwnershipUtils.hpp
 * @brief Utility functions that map CPSO sub-swarms to MPI ranks.
 *
 * The parallel CPSO assigns the global list of sub-swarms to the available
 * ranks as contiguous blocks. These helpers compute those ownership blocks and
 * translate a local batch position back to the corresponding global sub-swarm
 * index.
 */
#pragma once

#include <algorithm>
#include <utility>
#include <vector>

/**
 * @brief Splits the global sub-swarm list into contiguous ranges.
 * The distribution is as balanced as possible: the first remainder ranks get
 * one extra sub-swarm, while ranks beyond `total_swarms` receive an empty
 * interval.
 *
 * @param total_swarms Total number of CPSO sub-swarms.
 * @param mpi_size Total number of MPI ranks.
 * @return One half-open ownership interval per MPI rank.
 */
inline std::vector<std::pair<int, int>>
compute_all_swarm_ranges(int total_swarms, int mpi_size) {
  if (total_swarms <= 0 || mpi_size <= 0) {
    return {};
  }

  // Only the first min(mpi_size, total_swarms) ranks can own real work.
  const int effective_mpi_size = std::min(mpi_size, total_swarms);
  std::vector<int> swarms_per_proc(mpi_size, 0);
  const int base_count = total_swarms / effective_mpi_size;
  const int remainder_swarms = total_swarms % effective_mpi_size;

  // Distribute the leftover sub-swarms one by one to the first ranks.
  for (int i = 0; i < effective_mpi_size; ++i) {
    swarms_per_proc[i] = base_count + (i < remainder_swarms ? 1 : 0);
  }

  std::vector<std::pair<int, int>> ranges(mpi_size, {0, 0});
  int current_start = 0;
  for (int rank = 0; rank < mpi_size; ++rank) {
    ranges[rank] = {current_start, current_start + swarms_per_proc[rank]};
    current_start += swarms_per_proc[rank];
  }

  return ranges;
}

/**
 * @brief Returns the contiguous ownership interval of one MPI rank.
 *
 * @param total_swarms Total number of CPSO sub-swarms.
 * @param mpi_rank Rank whose ownership interval is requested.
 * @param mpi_size Total number of MPI ranks.
 * @return The half-open range `[start, end)` owned by that rank, or `{0, 0}`
 *         when the rank is invalid.
 */
inline std::pair<int, int> compute_local_swarm_range(int total_swarms,
                                                      int mpi_rank,
                                                      int mpi_size) {
  std::vector<std::pair<int, int>> ranges =
      compute_all_swarm_ranges(total_swarms, mpi_size);

  if (mpi_rank < 0 || mpi_rank >= static_cast<int>(ranges.size())) {
    return {0, 0};
  }

  return ranges[mpi_rank];
}

/**
 * @brief Maps a local batch position of one rank back to the global sub-swarm index.
 *
 * During the parallel CPSO loop, each rank iterates over its owned sub-swarms by
 * local batch position. This helper converts that local position
 * into the global sub-swarm index used by the shared arrays.
 *
 * @param ranges Ownership intervals previously returned by `compute_all_swarm_ranges`.
 * @param rank Rank whose local batch position is being converted.
 * @param batch_index Zero-based position inside the owned local interval.
 * @return The corresponding global sub-swarm index, or `-1` when that batch does
 *         not exist for the selected rank.
 */
inline int compute_batch_swarm_index(
    const std::vector<std::pair<int, int>> &ranges, int rank, int batch_index) {
  if (rank < 0 || rank >= static_cast<int>(ranges.size()) || batch_index < 0) {
    return -1;
  }

  // Convert the local batch position to the matching global swarm index.
  const int swarm_idx = ranges[rank].first + batch_index;
  return swarm_idx < ranges[rank].second ? swarm_idx : -1;
}
