#pragma once

#include <cstdlib>
#include <vector>
#include "data.h"

struct variogram_data
{
    size_t num_bins;
    std::vector<size_t> num_pairs;
    std::vector<double> distance_averages;
    std::vector<double> gamma;
    double max_distance;
    double min_distance;
};

/**
 * @brief empirical_variogram Computes the empirical variogram associated with
 * the given data points and the number of distance bins.
 * @param data_points A vector of data points to use in the computation.
 * @param num_bins The number of distance bins in the resulting variogram.
 * @return
 */
variogram_data empirical_variogram(const std::vector<data_point> & data_points,
                                   size_t num_bins);
