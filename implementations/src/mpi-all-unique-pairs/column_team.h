#pragma once

#include <mpi.h>

#include "data.h"
#include "parallel_options.h"

namespace pev {

/**
 * @brief Represents a column computation team, in the sense that the members of the
 * team are located in the same column of the computational grid. The topmost
 * processor is designated as the team leader, whose responsibility it is to
 * distribute the data to his team. All members of the same team share
 * the same data points for which they together interact with all other teams to
 * ensure that every point in the local team buffer interacts with every other point
 * in the other teams' buffers.
 */
class column_team
{
public:
    column_team();
    column_team(MPI_Comm active_comm, MPI_Datatype datatype, int replication_factor, int team_index);
    column_team(MPI_Comm active_comm, MPI_Datatype datatype, parallel_options options, int team_index);

    /**
     * @brief broadcast Broadcasts the results of a chunked read from the leader to the rest of the team.
     * @param result The result, which is modified to hold the results of the team leader's read operation.
     */
    void broadcast(chunked_read_result & result);

    /**
     * @brief my_rank Gives the local rank of the current processor within the team.
     * @return The local rank of the processor within the team, or -1 if not part of the team.
     */
    int my_rank() const;

    /**
     * @brief my_node_is_leader Determines whether the current processor is the leader of the team.
     * @return True if leader, false otherwise.
     */
    bool my_node_is_leader() const;

    std::vector<int> ranks() const;

private:
    MPI_Comm _comm;
    MPI_Datatype _datatype;
};

}
