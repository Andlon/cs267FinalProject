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
#include "grid.h"

// Anonymous namespace to hold functions that should not be exported,
// i.e. functions that are local to this file only.
namespace {

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
    const auto interval = variogram.max_distance / variogram.bin_count + eps;
    auto compute_bin = [interval] (double distance) -> size_t {
        return floor(distance / interval);
    };

    for (size_t i = 0; i < local_buffer.size(); ++i)
    {
        const auto local_point = local_buffer[i];
        for (size_t j = 0; j < exchange_buffer.size(); ++j)
        {
            const auto exchange_point = exchange_buffer[j];
            double dist = distance(local_point, exchange_point);
            double gamma = compute_gamma_contribution(local_point, exchange_point);

            auto bin = compute_bin(dist);
            assert(bin < variogram.bin_count);
            variogram.distance_averages[bin] += dist;
            variogram.num_pairs[bin] += 1;
            variogram.gamma[bin] += gamma;
        }
    }
}

void finalize_data(variogram_data & data)
{
    // Adjust for the fact that we've added self-interactions,
    // which would be n interactions of distance zero (and contribution zero),
    // where n is the number of points in the data set.
    data.num_pairs[0] -= data.point_count;

    for (size_t i = 0; i < data.bin_count; ++i) {
        data.gamma[i] /= (2 * data.num_pairs[i]);
        data.distance_averages[i] /= data.num_pairs[i];
        data.num_pairs[i] /= 2;
    }
}

} // End anonymous namespace

variogram_data::variogram_data()
    : bin_count(0), point_count(0), max_distance(0.0)
{

}

variogram_data::variogram_data(size_t num_bins)
    : bin_count(num_bins), point_count(0), max_distance(0.0)
{
    this->gamma = std::vector<double>(num_bins, 0);
    this->distance_averages = std::vector<double>(num_bins, 0);
    this->num_pairs = std::vector<size_t>(num_bins, 0u);
    this->bin_count = num_bins;
}

variogram_data empirical_variogram_parallel(const std::string &input_file, parallel_options options, size_t num_bins)
{
    node_grid grid(options);
    int P = options.active_processor_count();
    int c = options.replication_factor();

    // Kick out inactive processors
    if (!grid.local_is_active())
        return variogram_data();

    // Read and distribute data within each team
    column_team my_team = grid.create_team();
    parallel_read_result result;
    if (my_team.my_node_is_leader())
    {
        result = read_file_chunk_parallel(
                    input_file,
                    grid.team_count(),
                    grid.local_column());

        // Determine global bounds and use diagonal for global max distance
        auto bounds = bounding_rectangle(result.data);
        auto global_bounds = grid.reduce_bounds(bounds);
        result.max_distance = global_bounds.diagonal();
    }

    my_team.broadcast(result);

    const std::vector<data_point> & local_buffer = result.data;
    std::vector<data_point> exchange_buffer = local_buffer;
    variogram_data local_variogram(num_bins);
    local_variogram.max_distance = result.max_distance;
    local_variogram.point_count = result.global_point_count;

    // Given kth row processor, shift exchange_buffer by k along row
    int k = my_team.my_rank();
    exchange_buffer = grid.shift_along_row(std::move(exchange_buffer), k);

    int steps = P / (c * c);
    for (int s = 0; s < steps; ++s)
    {
        // Shift by c and compute interactions between the two buffers
        exchange_buffer = grid.shift_along_row(std::move(exchange_buffer), c);
        compute_contribution(local_variogram, local_buffer, exchange_buffer);
    }

    auto global_variogram = grid.reduce_variogram(local_variogram);
    finalize_data(global_variogram);
    return global_variogram;
}
