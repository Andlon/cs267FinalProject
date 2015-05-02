#pragma once

#include <mpi.h>
#include <algorithm>
#include <stdexcept>

inline MPI_Group group_from_comm(MPI_Comm comm)
{
    MPI_Group group;
    MPI_Comm_group(comm, &group);
    return group;
}

inline MPI_Group consecutive_subset_group(MPI_Group group, int a, int b)
{
    if (a < 0)
        throw std::logic_error("given interval [a, b) must only contain non-negative members");
    if (b < a)
        throw std::logic_error("given interval [a, b) must have b >= a");

    std::vector<int> ranks(b - a);
    std::iota(ranks.begin(), ranks.end(), a);

    MPI_Group new_group;
    MPI_Group_incl(group, ranks.size(), ranks.data(), &new_group);
    return new_group;
}

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
