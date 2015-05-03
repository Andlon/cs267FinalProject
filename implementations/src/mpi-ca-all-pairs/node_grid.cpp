#include "node_grid.h"
#include <mpi_util.h>
#include <vector>
#include <cassert>

namespace pev {

namespace {

int mod(int a, int b)
{
    if (b == 0)
        throw std::logic_error("Can not take (a mod b) with b zero.");
    int r = a % b;
    return r >= 0 ? r : r + std::abs(b);
}

MPI_Comm create_active_comm(size_t active_node_count)
{
    std::vector<int> active_ranks(active_node_count);
    std::iota(active_ranks.begin(), active_ranks.end(), 0);

    // Create a group that is the active subset of the complete set of available processors,
    // and associate a communicator with this group
    MPI_Group global_group = group_from_comm(MPI_COMM_WORLD);
    MPI_Group active_group = consecutive_subset_group(global_group, 0, active_node_count);

    MPI_Comm active_comm;
    MPI_Comm_create(MPI_COMM_WORLD, active_group, &active_comm);

    return active_comm;
}

MPI_Comm create_leader_comm(MPI_Comm active_comm, size_t leader_count)
{
    // Leaders are the first leader_count ranks in the set of active ranks
    MPI_Group active_group = group_from_comm(active_comm);
    MPI_Group leader_group = consecutive_subset_group(active_group, 0, leader_count);

    MPI_Comm leader_comm;
    MPI_Comm_create(active_comm, leader_group, &leader_comm);
    return leader_comm;
}

}

node_grid::node_grid(parallel_options options)
    : P(options.active_processor_count()),
      c(options.replication_factor()),
      _active_comm(MPI_COMM_NULL),
      _leader_comm(MPI_COMM_NULL),
      _data_point_type(create_mpi_data_point_type())
{
    _active_comm = create_active_comm(P);

    if (_active_comm != MPI_COMM_NULL)
        _leader_comm = create_leader_comm(_active_comm, P / c);
}

int node_grid::node_rank() const
{
    if (node_is_active())
    {
        int rank;
        MPI_Comm_rank(_active_comm, &rank);
        return rank;
    }
    else
        return -1;
}

int node_grid::node_column() const
{
    return node_is_active()
            ? node_rank() % team_count()
            : -1;
}

int node_grid::team_count() const
{
    return P / c;
}

bool node_grid::node_is_active() const
{
    return _active_comm != MPI_COMM_NULL;
}

bool node_grid::node_is_leader() const
{
    return _leader_comm != MPI_COMM_NULL;
}

column_team node_grid::create_team()
{
    int team_index = node_rank() % team_count();
    return column_team(_active_comm, _data_point_type, c, team_index);
}

pev::rect<double> node_grid::reduce_bounds(pev::rect<double> local_bounds)
{
    assert(_leader_comm != MPI_COMM_NULL);

    // Form a rectangle that covers the whole domain by reducing across leaders
    double x_min = local_bounds.x();
    double x_max = local_bounds.x() + local_bounds.width();
    double y_min = local_bounds.y();
    double y_max = local_bounds.y() + local_bounds.height();

    MPI_Allreduce(MPI_IN_PLACE, &x_min, 1, MPI_DOUBLE, MPI_MIN, _leader_comm);
    MPI_Allreduce(MPI_IN_PLACE, &x_max, 1, MPI_DOUBLE, MPI_MAX, _leader_comm);
    MPI_Allreduce(MPI_IN_PLACE, &y_min, 1, MPI_DOUBLE, MPI_MIN, _leader_comm);
    MPI_Allreduce(MPI_IN_PLACE, &y_max, 1, MPI_DOUBLE, MPI_MAX, _leader_comm);

    return pev::rect<double>::from_corners(x_min, x_max, y_min, y_max);
}

std::vector<data_point> node_grid::shift_along_row(std::vector<data_point> &&exchange_buffer, int shift_amount)
{
    int rank = node_rank();

    if (shift_amount == 0)
        return std::move(exchange_buffer);

    int col = rank % team_count();
    int row = rank / team_count();

    // Note: use mod function instead of % operator here,
    // as % is non-Euclidean
    int send_col = mod(col + shift_amount, team_count());
    int send_rank = row * team_count() + send_col;
    int recv_col = mod(col - shift_amount, team_count());
    int recv_rank = row * team_count() + recv_col;

    MPI_Request send_request = MPI_REQUEST_NULL;
    MPI_Isend(exchange_buffer.data(), exchange_buffer.size(),
              _data_point_type, send_rank, 0, _active_comm, &send_request);

    // Probe the received message to determine
    // how big to make our exchange buffer
    MPI_Status recv_status;
    MPI_Probe(recv_rank, 0, _active_comm, &recv_status);

    int recv_size;
    MPI_Get_count(&recv_status, _data_point_type, &recv_size);
    if (recv_size == MPI_UNDEFINED)
        throw std::runtime_error("unable to recover from undefined receive size");

    // Re-use member buffer
    //std::vector<data_point> & recv_buffer = _buffer;
    std::vector<data_point> recv_buffer;
    recv_buffer.resize(recv_size);

    MPI_Recv(recv_buffer.data(), recv_buffer.size(),
             _data_point_type, recv_rank, 0, _active_comm, MPI_STATUS_IGNORE);
    MPI_Wait(&send_request, MPI_STATUS_IGNORE);

    // recv_buffer contains the content we want to return, so swap
    //exchange_buffer.swap(recv_buffer);
    return recv_buffer;
}

std::vector<data_point> node_grid::shift_along_row(const std::vector<data_point> &exchange_buffer, int shift_amount)
{
    return shift_along_row(std::vector<data_point>(exchange_buffer), shift_amount);
}

variogram_data node_grid::reduce_variogram(variogram_data local_variogram)
{
    int num_bins = local_variogram.bin_count;
    variogram_data global_variogram = local_variogram;

    MPI_Reduce(local_variogram.distance_averages.data(), global_variogram.distance_averages.data(),
               num_bins, MPI_DOUBLE, MPI_SUM, 0, _active_comm);
    MPI_Reduce(local_variogram.gamma.data(), global_variogram.gamma.data(),
               num_bins, MPI_DOUBLE, MPI_SUM, 0, _active_comm);
    MPI_Reduce(local_variogram.num_pairs.data(), global_variogram.num_pairs.data(),
               num_bins, MPI_UINT64_T, MPI_SUM, 0, _active_comm);

    return global_variogram;
}

} // End namespace pev
