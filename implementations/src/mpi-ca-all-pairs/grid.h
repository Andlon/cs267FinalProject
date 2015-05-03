#pragma once

#include <rect.h>
#include <mpi.h>

#include "parallel_options.h"
#include "team.h"
#include "data.h"
#include "variogram.h"

class node_grid {
public:
    explicit node_grid(parallel_options options);

    bool local_is_active() const;
    bool local_is_leader() const;
    int local_rank() const;
    int local_column() const;

    int team_count() const;

    column_team create_team();

    /**
     * @brief reduce_bounds Performs an all-to-all reduction of each
     * leader's bounds, returning global bounds for the complete
     * data set.
     * @param local_bounds Bounds for the local data set.
     * @return Bounds for the global data set.
     */
    custom::rect<double> reduce_bounds(custom::rect<double> local_bounds);

    std::vector<data_point> shift_along_row(std::vector<data_point> && exchange_buffer,
                                            int shift_amount);

    variogram_data reduce_variogram(variogram_data local_variogram);


private:
    MPI_Datatype _data_point_type;
    MPI_Comm _active_comm;
    MPI_Comm _leader_comm;

    int P;
    int c;

    std::vector<data_point> _buffer;
};
