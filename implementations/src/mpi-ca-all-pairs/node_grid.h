#pragma once

#include <rect.h>
#include <mpi.h>

#include "parallel_options.h"
#include "column_team.h"
#include "data.h"
#include "variogram.h"

namespace pev {

/**
 * @brief An abstraction of the (P / c) x c computation grid used in the
 * Communication-Avoiding N-body algorithm. Handles all communication with MPI,
 * so that the high-level algorithm can succinctly be described in terms
 * of function calls into node_grid.
 *
 * The grid is represented from the perspective of the processor the code is
 * currently executing on. Hence function calls give information appropriate
 * to the current node.
 *
 * Terminology:
 * The grid is divided into (P / c) teams of c processors each,
 * with the top-most processor being the leader of that team.
 * Here P is referring to the number of active processors, and
 * c is the replication factor.
 */
class node_grid {
public:
    explicit node_grid(parallel_options options);

    /**
     * @brief node_is_active Determines whether or not the current node is an active
     * participant in the computation.
     * @return True if active, otherwise false.
     */
    bool node_is_active() const;

    /**
     * @brief node_is_leader Determines whether or not the current node is one of the
     * team leaders.
     * @return True if team leader, otherwise not.
     */
    bool node_is_leader() const;

    /**
     * @brief node_rank Gives the rank of the current processor in the active communicator.
     * @return The rank of the current processor, or -1 if the current processor is not active.
     */
    int node_rank() const;

    /**
     * @brief node_column Returns the column in the grid corresponding to the current processor.
     * @return The column of the current processor, or -1 if not active.
     */
    int node_column() const;

    /**
     * @brief team_count Gives the number of teams (columns) in the grid.
     * @return Number of teams/columns.
     */
    int team_count() const;

    /**
     * @brief create_team Creates a team for the current processor.
     * @return
     */
    column_team create_team();

    /**
     * @brief reduce_bounds Performs an all-to-all reduction of each
     * leader's bounds, returning global bounds for the complete
     * data set.
     * @param local_bounds Bounds for the local data set.
     * @return Bounds for the global data set.
     */
    pev::rect<double> reduce_bounds(pev::rect<double> local_bounds);

    /**
     * @brief shift_along_row Exchanges particles with other processors in the grid
     * by sending the contents of the current exchange_buffer to the processor
     * determined by shifting shift_amount to the right within the same row,
     * and receiving the contents of the processor determined by shifting
     * shift_amount to the left within the same row. The returned buffer contains
     * the received data points.
     *
     * Note that this overload takes an rvalue for the exchange_buffer,
     * reusing the allocated memory in the given buffer if possible.
     *
     * @param exchange_buffer The buffer containing data to send.
     * @param shift_amount The number of columns to shift along the row.
     * @return A buffer containing the data received.
     */
    std::vector<data_point> shift_along_row(std::vector<data_point> && exchange_buffer,
                                            int shift_amount);

    /**
     * @brief shift_along_row Exchanges particles with other processors in the grid
     * by sending the contents of the current exchange_buffer to the processor
     * determined by shifting shift_amount to the right within the same row,
     * and receiving the contents of the processor determined by shifting
     * shift_amount to the left within the same row. The returned buffer contains
     * the received data points.
     *
     * Note that this overload internally copies the contents of the exchange_buffer.
     * Use the rvalue overload if you want to avoid redundant allocations.
     *
     * @param exchange_buffer The buffer containing data to send.
     * @param shift_amount The number of columns to shift along the row.
     * @return A buffer containing the data received.
     */
    std::vector<data_point> shift_along_row(const std::vector<data_point> & exchange_buffer,
                                            int shift_amount);

    /**
     * @brief reduce_variogram Performs a reduction of partial variograms across
     * all active processors to produce a single (unfinalized) variogram on the
     * root node (rank 0). The output on all other processors is undefined
     * and should be discarded.
     *
     * Note that the variogram still needs to be finalized on the root node
     * (i.e. values averaged etc.) before the output is useful.
     *
     * @param local_variogram
     * @return A global (partial) variogram if the local node is the root node, otherwise the
     * contents of the variogram is undefined.
     */
    variogram_data reduce_variogram(variogram_data local_variogram);


private:
    MPI_Datatype _data_point_type;
    MPI_Comm _active_comm;
    MPI_Comm _leader_comm;

    int P;
    int c;

    std::vector<data_point> _buffer;
};

}
