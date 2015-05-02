#pragma once

#include <mpi.h>

#include "data.h"
#include "parallel_options.h"

class team
{
public:
    team();
    team(MPI_Comm active_comm, parallel_options options, int team_index);

    void broadcast(parallel_read_result & result, MPI_Datatype data_point_type);

    int my_rank() const;
    bool my_node_is_leader() const;
    bool my_node_belongs() const;

    std::vector<int> ranks() const;

private:
    MPI_Comm _comm;
};
