#include "variogram.h"
#include <pair_index_set.h>
#include <mpi_util.h>
#include <mpi.h>
#include <rect.h>
#include "column_team.h"
#include "node_grid.h"

#include <limits>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>

namespace pev {

// Anonymous namespace to hold functions that should not be exported,
// i.e. functions that are local to this file only.
namespace {

double compute_gamma_contribution(const data_point &p1, const data_point &p2)
{
    return pow(p2.value - p1.value, 2);
}

class compute_helper
{
public:
    explicit compute_helper(variogram_data & variogram)
        :   variogram(variogram)
    {
        // Note that we perturb the size of the maximum distance when computing the
        // interval for handling cases where the largest distance might not
        // "floor down", yielding an incorrect bin.
        const auto eps = 1e-6 * variogram.max_distance;
        const auto interval = (variogram.max_distance + eps) / variogram.bin_count;
        bin_func = [interval] (double distance) -> size_t {
            return floor(distance / interval);
        };
    }

    size_t compute_bin(double dist) { return bin_func(dist); }

    void interact(const data_point &p1, const data_point &p2)
    {
        auto dist = distance(p1, p2);
        auto gamma = compute_gamma_contribution(p1, p2);

        auto bin = compute_bin(dist);
        assert(bin < variogram.bin_count);
        variogram.distance_averages[bin] += dist;
        variogram.num_pairs[bin] += 1;
        variogram.gamma[bin] += gamma;
    }

private:
    variogram_data & variogram;
    std::function<size_t(double)> bin_func;
};


variogram_data compute_partial_contribution(variogram_data variogram,
                                            const std::vector<data_point> &local_buffer,
                                            const std::vector<data_point> &exchange_buffer,
                                            size_t a, size_t b)
{
    // Distribute the interactions as N buckets of M interactions
    size_t N = local_buffer.size();
    size_t M = exchange_buffer.size();

    // By using closed intervals for buckets (local particles)
    // and half-open intervals for indices within buckets (exchange particles),
    // we make sure that each processor computes a disjoint subset,
    // with the union of all subsets equal to the entire set.

    // Bucket interval [start_i, end_i] is closed
    size_t start_i = a / M;
    size_t end_i = (b - 1) / M;

    // Index interval within bucket [start_j, end_j) is half-open.
    size_t start_j = a - start_i * M;
    size_t end_j = b - end_i * M;

    if (b <= a)
        throw std::logic_error("Interval must be non-empty.");

    compute_helper comp(variogram);

    for (size_t i = start_i; i <= end_i; ++i)
    {
        auto local_point = local_buffer[i];
        size_t begin = i == start_i ? start_j : 0;
        size_t end = i == end_i ? end_j : M;
        for (size_t j = begin; j < end; ++j)
            comp.interact(local_point, exchange_buffer[j]);
    }

    return variogram;
}

variogram_data compute_self_contribution(variogram_data variogram,
                                         const std::vector<data_point> & data)
{
    // Note that we perturb the size of the maximum distance when computing the
    // interval for handling cases where the largest distance might not
    // "floor down", yielding an incorrect bin.
    const auto eps = 1e-6 * variogram.max_distance;
    const auto interval = (variogram.max_distance + eps) / variogram.bin_count;
    auto compute_bin = [interval] (double distance) -> size_t {
        return floor(distance / interval);
    };

    for (size_t i = 0; i < data.size(); ++i)
    {
        const auto p1 = data[i];
        for (size_t j = i + 1; j < data.size(); ++j)
        {
            const auto p2 = data[j];

            double dist = distance(p1, p2);
            double gamma = compute_gamma_contribution(p1, p2);

            auto bin = compute_bin(dist);
            assert(bin < variogram.bin_count);
            variogram.distance_averages[bin] += dist;
            variogram.num_pairs[bin] += 1;
            variogram.gamma[bin] += gamma;
        }
    }

    return variogram;
}

/**
 * @brief compute_contribution
 * @param variogram
 * @param local_buffer
 * @param exchange_buffer
 */
variogram_data compute_contribution(variogram_data variogram,
                                    const std::vector<data_point> & local_buffer,
                                    const std::vector<data_point> & exchange_buffer)
{
    size_t interaction_count = local_buffer.size() * exchange_buffer.size();
    return compute_partial_contribution(variogram, local_buffer,
                                        exchange_buffer, 0, interaction_count);
}

void finalize_data(variogram_data & data)
{
    // Adjust for the fact that we've added self-interactions,
    // which would be n interactions of distance zero (and contribution zero),
    // where n is the number of points in the data set.

    for (size_t i = 0; i < data.bin_count; ++i) {
        data.gamma[i] /= (2 * data.num_pairs[i]);
        data.distance_averages[i] /= data.num_pairs[i];
    }
}

} // End anonymous namespace

variogram_data::variogram_data()
    : bin_count(0), point_count(0), max_distance(0.0),
      contains_timing(false),
      options(1, 1)
{

}

variogram_data::variogram_data(size_t num_bins)
    : variogram_data()
{
    this->bin_count = num_bins;
    this->gamma = std::vector<double>(num_bins, 0);
    this->distance_averages = std::vector<double>(num_bins, 0);
    this->num_pairs = std::vector<size_t>(num_bins, 0u);
    this->bin_count = num_bins;
}

variogram_data empirical_variogram_parallel(const std::string &input_file, parallel_options options, size_t num_bins)
{
    using std::chrono::steady_clock;
    typedef std::chrono::duration<double> fsecs;

    node_grid grid(options);
    int P = options.active_processor_count();
    int c = options.replication_factor();

    // Kick out inactive processors
    if (!grid.node_is_active())
        return variogram_data();

    auto input_read_start = steady_clock::now();
    // Read and distribute data within each team
    column_team my_team = grid.create_team();
    chunked_read_result result;
    if (my_team.my_node_is_leader())
    {
        result = read_file_chunk_parallel(
                    input_file,
                    grid.team_count(),
                    grid.node_column());

        // Determine global bounds and use diagonal for global max distance
        auto bounds = bounding_rectangle(result.data);
        auto global_bounds = grid.reduce_bounds(bounds);
        result.max_distance = global_bounds.diagonal();
    }
    auto input_read_time = steady_clock::now() - input_read_start;

    auto input_broadcast_start = steady_clock::now();
    my_team.broadcast(result);
    auto input_broadcast_time = steady_clock::now() - input_broadcast_start;

    const std::vector<data_point> & local_buffer = result.data;
    std::vector<data_point> exchange_buffer = local_buffer;
    variogram_data local_variogram(num_bins);
    local_variogram.max_distance = result.max_distance;
    local_variogram.point_count = result.global_point_count;

    // Given kth row processor, shift exchange_buffer by k along row
    int k = my_team.my_rank();
    auto shift_start = steady_clock::now();
    exchange_buffer = grid.shift_along_row(exchange_buffer, k);
    auto shift_time = steady_clock::now() - shift_start;

    int steps = P / (c * c);
    steady_clock::duration computation_time = steady_clock::duration::zero();
    for (int s = 0; s < steps; ++s)
    {
        // Shift by c and compute interactions between the two buffers
        shift_start = steady_clock::now();
        exchange_buffer = grid.shift_along_row(exchange_buffer, c);
        shift_time += steady_clock::now() - shift_start;

        auto computation_start = steady_clock::now();
        local_variogram = compute_contribution(local_variogram, local_buffer, exchange_buffer);
        computation_time += steady_clock::now() - computation_start;
    }

    auto reduction_start = steady_clock::now();
    auto global_variogram = grid.reduce_variogram(local_variogram);
    auto reduction_time = steady_clock::now() - reduction_start;

    finalize_data(global_variogram);

    // Append timing and globally reduce timing information
    // TODO: Add switch for whether timing is necessary
    timing_info timing;
    timing.input_read_time = fsecs(input_read_time).count();
    timing.input_broadcast_time = fsecs(input_broadcast_time).count();
    timing.shifting_time = fsecs(shift_time).count();
    timing.computation_time = fsecs(computation_time).count();
    timing.reduction_time = fsecs(reduction_time).count();
    global_variogram.timing = grid.reduce_timing(timing);
    global_variogram.contains_timing = true;
    global_variogram.options = options;
    return global_variogram;
}

std::string timing_info::format_description()
{
    return "Appends one line to the given timing output file: \n\n"
           "<N> <P> <c> <input read> <input broadcast> <shift> <computation> <reduction>\n\n"
           "Here N is the global number of points, P is the number of \n"
           "ACTIVE processors and c is the replication factor, \n"
           "with the remaining parameters being time measurements. \n"
           "All timing values are measured in fractional seconds. The sum of \n"
           "the values constitutes roughly the total time. Values are delimited \n"
           "by TAB characters.";
}

void print_timing_info(std::ostream &out, const variogram_data &variogram)
{
    // NOTE: This must be kept up-to-date with timing_info::format_description
    out << variogram.point_count << "\t"
        << variogram.options.active_processor_count() << "\t"
        << variogram.options.replication_factor() << "\t"
        << variogram.timing
        << std::endl;
}

std::ostream &operator <<(std::ostream &out, const timing_info &timing_info)
{
    out << timing_info.input_read_time << "\t"
        << timing_info.input_broadcast_time << "\t"
        << timing_info.shifting_time << "\t"
        << timing_info.computation_time << "\t"
        << timing_info.reduction_time;
    return out;
}

std::string variogram_data::format_description()
{
    return "Truncates the given output file and writes one line per \n"
           "bin in the variogram. Each line has the form: \n\n"
           "<pairs in bin> <distance average for bin> <gamma for bin>\n\n"
           "Values are delimited by TAB characters.";
}

void print_variogram(std::ostream &out, const variogram_data &variogram)
{
    // NOTE: This must be kept up-to-date with variogram_data::format_description
    for(int i = 0; i < variogram.bin_count; i++)
    {
        out << variogram.num_pairs[i] << "\t"
            << variogram.distance_averages[i] << "\t"
            << variogram.gamma[i] << std::endl;
    }
    out << std::endl;
}

variogram_data empirical_variogram_parallel(const std::string &input_file, size_t num_bins)
{
    int P;
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    // We can re-use the grid functionality by setting c = 1
    node_grid grid(parallel_options(P, 1));
    int p = grid.node_column();

    // Read data and determine global diagonal
    chunked_read_result result = read_file_chunk_parallel(input_file, P, p);
    auto bounds = bounding_rectangle(result.data);
    result.max_distance = grid.reduce_bounds(bounds).diagonal();

    const std::vector<data_point> & local_buffer = result.data;
    auto exchange_buffer = local_buffer;
    variogram_data local_variogram(num_bins);
    local_variogram.max_distance = result.max_distance;
    local_variogram.point_count = result.global_point_count;

    // First interactions between particles in the same chunk
    local_variogram = compute_self_contribution(local_variogram, local_buffer);

    // Do floor(P / 2) - 1 steps of shifting and computing contributions of disjoint data sets
    size_t steps = P % 2 == 0 ? P / 2 - 1 : P / 2;
    for (size_t i = 0; i < steps; ++i)
    {
        exchange_buffer = grid.shift_along_row(exchange_buffer, 1);
        local_variogram = compute_contribution(local_variogram, local_buffer, exchange_buffer);
    }

    if (P % 2 == 0)
    {
        // If P % 2 the two processors which share the same set of data points
        // need to compute half the number of interactions each. So we split the number of computations
        // to be performed roughly equally across the two processors by
        // dividing the interactions instead of the data. We create a
        // lower interval and an upper interval, corresponding to
        // the first half and second half of the interaction set respectively.

        exchange_buffer = grid.shift_along_row(exchange_buffer, 1);
        size_t interaction_count = local_buffer.size() * exchange_buffer.size();

        size_t lower_interval_start = 0;
        size_t lower_interval_end = interaction_count / 2;
        size_t upper_interval_start = interaction_count / 2;
        size_t upper_interval_end = interaction_count;

        if (p < P / 2)
        {
            local_variogram = compute_partial_contribution(local_variogram,
                                                           local_buffer,
                                                           exchange_buffer,
                                                           lower_interval_start,
                                                           lower_interval_end);
        }
        else
        {
            local_variogram = compute_partial_contribution(local_variogram,
                                                           exchange_buffer,
                                                           local_buffer,
                                                           upper_interval_start,
                                                           upper_interval_end);
        }
    }

    auto global_variogram = grid.reduce_variogram(local_variogram);
    finalize_data(global_variogram);
    return global_variogram;
}

}
