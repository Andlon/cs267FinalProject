#pragma once

#include <mpi.h>
#include <algorithm>

/**
 * @brief rank_in_group Returns the rank of the calling processor in the given group,
 * or -1 if it does not belong to the group.
 * @param group
 * @return
 */
inline int rank_in_group(MPI_Group group)
{
    int rank;
    MPI_Group_rank(group, &rank);

    return rank != MPI_UNDEFINED
            ? rank
            : -1;
}

inline std::vector<int> all_ranks_in_group(MPI_Group group)
{
    int size;
    MPI_Group_size(group, &size);
    std::vector<int> ranks(size);
    std::iota(ranks.begin(), ranks.end(), 0);
    return ranks;
}

inline std::vector<int> all_global_ranks()
{
    MPI_Group global_group;
    MPI_Comm_group(MPI_COMM_WORLD, &global_group);
    return all_ranks_in_group(global_group);
}
