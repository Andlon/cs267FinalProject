#include "variogram.h"
#include <pair_index_set.h>
#include <mpi_util.h>

#include <limits>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include <functional>
#include "team.h"

// Anonymous namespace to hold functions that should not be exported,
// i.e. functions that are local to this file only.
namespace {

/**
 * Calls function f for every pair k of data points in data_points,
 * where k is an integer in the closed interval [a, b] where
 * b >= a and both a, b are members of the implicit index set.
 *
 * The expected signature of f is:
 * (const data_point &, const data_point&)
 */
template <typename Func>
void for_each_pair(const std::vector<data_point> & data_points,
                   const pair_index_set & index_set,
                   index_type a,
                   index_type b,
                   Func f)
{
    for (index_type k = a; k <= b; ++k)
    {
        auto pair = index_set.map_to_pair(k);
        const data_point & p1 = data_points[pair.first - 1];
        const data_point & p2 = data_points[pair.second - 1];
        f(p1, p2);
    }
}

struct distance_range
{
    double max;
    double min;

    distance_range(double min, double max)
        :   min(min), max(max) {}

    distance_range()
        :   min(0.0), max(0.0) {}
};

/**
 * @brief compute_distance_range Computes the maximum and minimum distances in the subset of
 * the given data set defined by the indices in the closed interval [a, b], which is a subset
 * of the supplied index set.
 * @param data_points
 * @param index_set
 * @param a
 * @param b
 * @return
 */
distance_range compute_distance_range(const std::vector<data_point> & data_points,
                                      const pair_index_set & index_set,
                                      index_type a,
                                      index_type b)
{
    double d_max = 0;
    double d_min = std::numeric_limits<double>::max();

    for_each_pair(data_points, index_set, a, b,
                  [&d_min, &d_max] (const auto &p1, const auto &p2)
    {
        double dist = distance(p1, p2);
        if (dist > d_max)
            d_max = dist;
        if (dist < d_min)
            d_min = dist;
    });

    return distance_range(d_min, d_max);
}

double compute_gamma_contribution(const data_point &p1, const data_point &p2)
{
    return pow(p2.value - p1.value, 2);
}

/**
 * @brief compute_partial_variogram Computes a partial empirical variogram over the
 * indices [a, b] from the supplied index set.
 *
 * Note that returned data is not averaged.
 *
 * @param data_points The data points in the data set.
 * @param index_set The index set from which to iterate over.
 * @param a Index of first interaction, inclusive.
 * @param b Index of last interaction, inclusive.
 * @param range The distance range information for the particle data.
 * @param num_bins Number of bins to use when placing interactions in bins based on distance.
 * @return Variogram data for the specified interval of interactions, without averaging.
 */
variogram_data compute_partial_variogram(const std::vector<data_point> & data_points,
                                         const pair_index_set & index_set,
                                         index_type a,
                                         index_type b,
                                         distance_range range,
                                         size_t num_bins)
{
    variogram_data data(num_bins);
    data.min_distance = range.min;
    data.max_distance = range.max;

    // Note that we perturb the size of the interval slightly to make sure that
    // the maximum distance falls within a valid interval. In this case we use some
    // arbitrary factor of machine epsilon.
    const double interval = (range.max - range.min) / num_bins + 10.0 * std::numeric_limits<double>::epsilon();
    auto compute_bin = [interval, range] (double distance) -> size_t {
        return floor((distance - range.min) / interval);
    };

    // Loop over all interaction pairs
    for_each_pair(data_points, index_set, a, b,
                  [&data, &compute_bin] (const auto &p1, const auto &p2)
    {
        double dist = distance(p1, p2);
        size_t bin = compute_bin(dist);

        assert(bin >= 0);
        assert(bin < data.num_bins);

        data.gamma[bin] += compute_gamma_contribution(p1, p2);
        data.num_pairs[bin] += 1;
        data.distance_averages[bin] += dist;
    });

    return data;
}

void average_data(variogram_data & data)
{
    for (size_t i = 0; i < data.num_bins; ++i) {
        data.gamma[i] /= (2 * data.num_pairs[i]);
        data.distance_averages[i] /= data.num_pairs[i];
    }
}

} // End anonymous namespace

variogram_data empirical_variogram(const std::vector<data_point> &data_points, size_t num_bins)
{
    // The index set lets us easily iterate over interactions
    // in a single loop, as opposed to a nested loop over data points
    pair_index_set index_set(data_points.size());
    distance_range range = compute_distance_range(data_points, index_set, 1, index_set.count());

    // In this case, we compute the entire variogram by computing a partial variogram over
    // the entire index set
    variogram_data data = compute_partial_variogram(data_points, index_set, 1, index_set.count(), range, num_bins);
    average_data(data);
    return data;
}

variogram_data::variogram_data()
    : num_bins(0), max_distance(0.0), min_distance(0.0)
{

}

variogram_data::variogram_data(size_t num_bins)
    : num_bins(num_bins), max_distance(0.0), min_distance(0.0)
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

    std::vector<int> active_ranks(P);
    std::iota(active_ranks.begin(), active_ranks.end(), 0);

    // Create a group that is the active subset of the complete set of available processors,
    // and associate a communicator with this group
    MPI_Group global_group, active_group;
    MPI_Comm_group(MPI_COMM_WORLD, &global_group);
    MPI_Group_incl(global_group, active_ranks.size(), active_ranks.data(), &active_group);

    MPI_Comm active_comm;
    MPI_Comm_create(MPI_COMM_WORLD, active_group, &active_comm);

    // Kick any inactive processors out
    if (active_comm == MPI_COMM_NULL)
        return variogram_data(num_bins);

    int active_rank;
    MPI_Comm_rank(active_comm, &active_rank);

    // Create my team
    int team_count = P / c;
    int my_team_index = active_rank % (team_count);
    team my_team(active_comm, options, my_team_index);

    parallel_read_result read_result;

    if (my_team.my_node_is_leader())
        read_result = read_file_chunk_parallel(input_file, team_count, my_team_index);

    my_team.broadcast(read_result, data_point_type);

    std::cout << "Processor " << active_rank << " has " << read_result.data.size() << " of "
              << read_result.global_point_count << " data points." << std::endl;

    MPI_Barrier(active_comm);

    variogram_data data(num_bins);
    return data;
}
