#pragma once

#include <cstdlib>
#include <vector>
#include <cstdint>
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

/**
 * @brief empirical_variogram_parallel Computes the empirical variogram associated with
 * the given data points and the number of distance bins in parallel.
 *
 * Note that this function uses MPI behind the scenes, and so it is necessary that
 * MPI_Init has been called before calling this function.
 *
 * Also note that the result of this function call is undefined on any other ranks
 * than the root node.
 *
 * @param data_points A vector of data points to use in the computation.
 * @param num_bins The number of distance bins in the resulting variogram.
 * @param root The rank of the processor where the resulting variogram will be stored.
 * The result is undefined if the rank is invalid.
 * @return
 */
variogram_data empirical_variogram_parallel(const std::vector<data_point> & data_points,
                                            size_t num_bins, int root = 0);
