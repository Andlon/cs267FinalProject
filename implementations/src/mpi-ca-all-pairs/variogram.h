#pragma once

#include <cstdlib>
#include <vector>
#include <cstdint>
#include "parallel_options.h"
#include "data.h"

struct variogram_data
{
    variogram_data();
    explicit variogram_data(size_t num_bins);

    size_t num_bins;
    std::vector<u_int64_t> num_pairs;
    std::vector<double> distance_averages;
    std::vector<double> gamma;
    double max_distance;
};

variogram_data empirical_variogram_parallel(const std::string & input_file, parallel_options options, size_t num_bins);
