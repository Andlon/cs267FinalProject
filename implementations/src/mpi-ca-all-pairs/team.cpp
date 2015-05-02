#include "team.h"
#include <mpi_util.h>
#include <algorithm>


team::team()
    : _comm(MPI_COMM_NULL)
{
}

team::team(MPI_Comm active_comm, parallel_options options, int team_index)
{
    int c = options.replication_factor();

    // Number of team members is equal to replication factor, c
    std::vector<int> team_ranks(c);
    std::iota(team_ranks.begin(), team_ranks.end(), 0);
    std::transform(team_ranks.begin(), team_ranks.end(), team_ranks.begin(),
                   [c, team_index] (auto i) { return team_index + i * c; });

    MPI_Group active_group;
    MPI_Comm_group(active_comm, &active_group);

    MPI_Group team_group;
    MPI_Group_incl(active_group, c, team_ranks.data(), &team_group);
    MPI_Comm_create(active_comm, team_group, &_comm);
}

void team::broadcast(parallel_read_result &result, MPI_Datatype data_point_type)
{
    MPI_Bcast(&result.global_point_count, 1, MPI_UINT64_T, 0, _comm);
    MPI_Bcast(&result.max_distance, 1, MPI_DOUBLE, 0, _comm);

    u_int64_t data_size = result.data.size();

    // TODO: Remove this broadcast. Can compute the data size using uniform_distribution
    MPI_Bcast(&data_size, 1, MPI_UINT64_T, 0, _comm);
    result.data.resize(data_size);
    MPI_Bcast(result.data.data(), data_size, data_point_type, 0, _comm);
}

int team::my_rank() const
{
    MPI_Group team_group;
    MPI_Comm_group(_comm, &team_group);
    return rank_in_group(team_group);
}

bool team::my_node_is_leader() const
{
    // Leader has relative rank 0
    return my_rank() == 0;
}

bool team::my_node_belongs() const
{
    return _comm != MPI_COMM_NULL;
}

std::vector<int> team::ranks() const
{
    MPI_Group global_group;
    MPI_Comm_group(MPI_COMM_WORLD, &global_group);

    MPI_Group team_group;
    MPI_Comm_group(_comm, &team_group);

    std::vector<int> team_ranks = all_ranks_in_group(team_group);
    std::vector<int> global_ranks(team_ranks.size());
    MPI_Group_translate_ranks(team_group, team_ranks.size(), team_ranks.data(),
                              global_group, global_ranks.data());

    return global_ranks;
}
