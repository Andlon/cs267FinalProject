#pragma once

#include <mpi.h>

#include "data.h"
#include "parallel_options.h"

class column_team
{
public:
    column_team();
    column_team(MPI_Comm active_comm, MPI_Datatype datatype, int replication_factor, int team_index);
    column_team(MPI_Comm active_comm, MPI_Datatype datatype, parallel_options options, int team_index);

    void broadcast(parallel_read_result & result);

    int my_rank() const;
    bool my_node_is_leader() const;
    bool my_node_belongs() const;

    std::vector<int> ranks() const;

private:
    MPI_Comm _comm;
    MPI_Datatype _datatype;
};
