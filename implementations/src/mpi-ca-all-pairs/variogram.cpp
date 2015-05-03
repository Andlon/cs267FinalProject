#include "variogram.h"
#include <pair_index_set.h>
#include <mpi_util.h>

#include <limits>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <utility>
#include <mpi.h>
#include <functional>
#include <rect.h>
#include "team.h"

// Anonymous namespace to hold functions that should not be exported,
// i.e. functions that are local to this file only.
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
    MPI_Group global_group, active_group;
    MPI_Comm_group(MPI_COMM_WORLD, &global_group);
    MPI_Group_incl(global_group, active_ranks.size(), active_ranks.data(), &active_group);

    MPI_Comm active_comm;
    MPI_Comm_create(MPI_COMM_WORLD, active_group, &active_comm);

    return active_comm;
}

MPI_Comm create_leader_comm(MPI_Comm active_comm, size_t leader_count)
{
    MPI_Group active_group = group_from_comm(active_comm);

    // Leaders are the first leader_count ranks in the set of active ranks
    MPI_Group leader_group = consecutive_subset_group(active_group, 0, leader_count);

    MPI_Comm leader_comm;
    MPI_Comm_create(active_comm, leader_group, &leader_comm);
    return leader_comm;
}

int shift_right(std::vector<data_point> & exchange_buffer,
                MPI_Datatype data_point_type,
                MPI_Comm comm,
                int rank,
                int team_count,
                int shift_amount)
{
    if (shift_amount == 0)
        return rank;

    int col = rank % team_count;
    int row = rank / team_count;

    // Note: use mod function instead of % operator here,
    // as % is non-Euclidean
    int send_col = mod(col + shift_amount, team_count);
    int send_rank = row * team_count + send_col;
    int recv_col = mod(col - shift_amount, team_count);
    int recv_rank = row * team_count + recv_col;

    // TODO: Consider taking a send buffer in so we don't
    // need to allocate one for every shift. Probably
    // miniscule on performance though, need to profile.
    std::vector<data_point> send_buffer = exchange_buffer;

    MPI_Request send_request = MPI_REQUEST_NULL;
    MPI_Isend(send_buffer.data(), send_buffer.size(),
              data_point_type, send_rank, 0, comm, &send_request);

    // Probe the received message to determine
    // how big to make our exchange buffer
    MPI_Status recv_status;
    MPI_Probe(recv_rank, 0, comm, &recv_status);

    int recv_size;
    MPI_Get_count(&recv_status, data_point_type, &recv_size);
    if (recv_size == MPI_UNDEFINED)
        throw std::runtime_error("unable to recover from undefined receive size");

    exchange_buffer.resize(recv_size);

    MPI_Recv(exchange_buffer.data(), exchange_buffer.size(),
             data_point_type, recv_rank, 0, comm, MPI_STATUS_IGNORE);
    MPI_Wait(&send_request, MPI_STATUS_IGNORE);

    return recv_rank;
}

double compute_gamma_contribution(const data_point &p1, const data_point &p2)
{
    return pow(p2.value - p1.value, 2);
}

void compute_contribution(variogram_data & variogram,
                          const std::vector<data_point> & local_buffer,
                          const std::vector<data_point> & exchange_buffer)
{
    // Note that we perturb the size of the interval slightly to make sure that
    // the maximum distance falls within a valid interval. In this case we use some
    // arbitrary factor of machine epsilon.
    const auto eps = 10.0 * std::numeric_limits<double>::epsilon();
    const auto interval = variogram.max_distance / variogram.num_bins + eps;
    auto compute_bin = [interval] (double distance) -> size_t {
        return floor(distance / interval);
    };

    for (size_t i = 0; i < local_buffer.size(); ++i)
    {
        const auto local_point = local_buffer[i];
        for (size_t j = 0; j < exchange_buffer.size(); ++j)
        {
            if (i != j)
            {
                const auto exchange_point = exchange_buffer[j];
                double dist = distance(local_point, exchange_point);
                double gamma = compute_gamma_contribution(local_point, exchange_point);

                auto bin = compute_bin(dist);
                assert(bin < variogram.num_bins);
                variogram.distance_averages[bin] += dist;
                variogram.num_pairs[bin] += 1;
                variogram.gamma[bin] += gamma;
            }
        }
    }
}

void finalize_data(variogram_data & data)
{
    for (size_t i = 0; i < data.num_bins; ++i) {
        // Correct the fact that we've counted pairs twice
        data.num_pairs[i] /= 2;

        // Note the extra factor of 2, which appears because
        // we're doing redundant work and the value of gamma is
        // twice as large as it should be
        data.gamma[i] /= (2 * 2 * data.num_pairs[i]);
        data.distance_averages[i] /= data.num_pairs[i];
    }
}

} // End anonymous namespace

variogram_data::variogram_data()
    : num_bins(0), max_distance(0.0)
{

}

variogram_data::variogram_data(size_t num_bins)
    : num_bins(num_bins), max_distance(0.0)
{
    this->gamma = std::vector<double>(num_bins, 0);
    this->distance_averages = std::vector<double>(num_bins, 0);
    this->num_pairs = std::vector<size_t>(num_bins, 0u);
    this->num_bins = num_bins;
}

variogram_data empirical_variogram_parallel(const std::string &input_file, parallel_options options, size_t num_bins)
{
    MPI_Datatype data_point_type = create_mpi_data_point_type();

    // The variable names here are consistent with the terminology used in
    // Driscoll et al. "A Communication-Optimal N-Body Algorithm for Direct Interactions"
    int P = options.active_processor_count();
    int c = options.replication_factor();
    assert(P % c == 0);

    MPI_Comm active_comm = create_active_comm(P);

    // Kick any inactive processors out
    if (active_comm == MPI_COMM_NULL)
        return variogram_data(num_bins);

    MPI_Comm leader_comm = create_leader_comm(active_comm, P / c);

    int active_rank;
    MPI_Comm_rank(active_comm, &active_rank);

    // Create my team
    int team_count = P / c;
    int my_team_index = active_rank % (team_count);
    team my_team(active_comm, options, my_team_index);

    parallel_read_result read_result;
    if (my_team.my_node_is_leader()) {
        read_result = std::move(read_file_chunk_parallel(input_file, team_count, my_team_index));

        // Find local max distance and max-reduce among all leaders
        auto bounds = bounding_rectangle(read_result.data);

        // Form a rectangle that covers the whole domain by reducing across leaders
        double x_min = bounds.x();
        double x_max = bounds.x() + bounds.width();
        double y_min = bounds.y();
        double y_max = bounds.y() + bounds.height();

        MPI_Allreduce(MPI_IN_PLACE, &x_min, 1, MPI_DOUBLE, MPI_MIN, leader_comm);
        MPI_Allreduce(MPI_IN_PLACE, &x_max, 1, MPI_DOUBLE, MPI_MAX, leader_comm);
        MPI_Allreduce(MPI_IN_PLACE, &y_min, 1, MPI_DOUBLE, MPI_MIN, leader_comm);
        MPI_Allreduce(MPI_IN_PLACE, &y_max, 1, MPI_DOUBLE, MPI_MAX, leader_comm);

        auto global_bounds = custom::rect<double>::from_corners(x_min, x_max, y_min, y_max);
        read_result.max_distance = global_bounds.diagonal();
    }

    my_team.broadcast(read_result, data_point_type);

    const std::vector<data_point> & local_buffer = read_result.data;
    std::vector<data_point> exchange_buffer = local_buffer;
    variogram_data local_variogram(num_bins);
    local_variogram.max_distance = read_result.max_distance;

    // Shift right by k where k is equal to the row of the current processor
    shift_right(exchange_buffer, data_point_type,
                active_comm, active_rank,
                team_count, my_team.my_rank());

    int steps = P / (c * c);
    for (int s = 0; s < steps; ++s)
    {
        // Shift by c
        shift_right(exchange_buffer, data_point_type,
                    active_comm, active_rank,
                    team_count, c);
        compute_contribution(local_variogram, local_buffer, exchange_buffer);
    }

    variogram_data global_variogram(num_bins);
    global_variogram.max_distance = local_variogram.max_distance;

    MPI_Reduce(local_variogram.distance_averages.data(), global_variogram.distance_averages.data(),
               num_bins, MPI_DOUBLE, MPI_SUM, 0, active_comm);
    MPI_Reduce(local_variogram.gamma.data(), global_variogram.gamma.data(),
               num_bins, MPI_DOUBLE, MPI_SUM, 0, active_comm);
    MPI_Reduce(local_variogram.num_pairs.data(), global_variogram.num_pairs.data(),
               num_bins, MPI_UINT64_T, MPI_SUM, 0, active_comm);

    finalize_data(global_variogram);

    return global_variogram;
}
