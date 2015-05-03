#include "team.h"
#include <mpi_util.h>
#include <algorithm>


column_team::column_team()
    : _comm(MPI_COMM_NULL)
{
}

column_team::column_team(MPI_Comm active_comm, MPI_Datatype datatype, int replication_factor, int team_index)
    :   _datatype(datatype)
{
    int c = replication_factor;

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

column_team::column_team(MPI_Comm active_comm, MPI_Datatype datatype, parallel_options options, int team_index)
    : column_team(active_comm, datatype, options.replication_factor(), team_index)
{
}

void column_team::broadcast(parallel_read_result &result)
{
    u_int64_t data_size = result.data.size();
    MPI_Bcast(&result.global_point_count, 1, MPI_UINT64_T, 0, _comm);
    MPI_Bcast(&result.max_distance, 1, MPI_DOUBLE, 0, _comm);
    MPI_Bcast(&data_size, 1, MPI_UINT64_T, 0, _comm);
    result.data.resize(data_size);
    MPI_Bcast(result.data.data(), data_size, _datatype, 0, _comm);
}

int column_team::my_rank() const
{
    MPI_Group team_group;
    MPI_Comm_group(_comm, &team_group);
    return rank_in_group(team_group);
}

bool column_team::my_node_is_leader() const
{
    // Leader has relative rank 0
    return my_rank() == 0;
}

bool column_team::my_node_belongs() const
{
    return _comm != MPI_COMM_NULL;
}

std::vector<int> column_team::ranks() const
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
