#pragma once

#include <cstdlib>
#include <vector>
#include <cstdint>
#include "parallel_options.h"
#include "data.h"

namespace pev {

struct variogram_data
{
    variogram_data();
    explicit variogram_data(size_t bin_count);

    size_t bin_count;
    size_t point_count;
    std::vector<u_int64_t> num_pairs;
    std::vector<double> distance_averages;
    std::vector<double> gamma;
    double max_distance;
};

/**
 * @brief empirical_variogram_parallel Computes an empirical variogram of the data in the
 * given input file by distributing the workload in the manner described by the supplied
 * parallel_options, using the given number of bins to categorize the data points.
 *
 * The returned data is only valid on the root node (rank 0), and undefined for
 * all other processes.
 *
 * @param input_file The path to the input file in question.
 * @param options The parallel strategy to adhere to.
 * @param num_bins The number of distance bins to place data points in.
 * @return Valid variogram data if local node is the root node, otherwise undefined.
 */
variogram_data empirical_variogram_parallel(const std::string & input_file, parallel_options options, size_t num_bins);

}
