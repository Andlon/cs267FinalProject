#include "variogram.h"
#include "pair_index_set.h"

#include <limits>
#include <cassert>
#include <iostream>

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
};

/**
 * @brief compute_distance_range Computes the maximum and minimum distanes in the subset of
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

} // End anonymous namespace

variogram_data empirical_variogram(const std::vector<data_point> &data_points, size_t num_bins)
{
    variogram_data data;
    data.gamma = std::vector<double>(num_bins, 0);
    data.distance_averages = std::vector<double>(num_bins, 0);
    data.num_pairs = std::vector<size_t>(num_bins, 0u);
    data.num_bins = num_bins;

    // The index set lets us easily iterate over interactions
    // in a single loop, as opposed to a nested loop over data points
    pair_index_set index_set(data_points.size());

    // Note that we add a very tiny number, a fraction of the difference between the maximum and minimum
    // distance, to the interval calculation to avoid points at the max boundary to be placed
    // in a non-existent bin.
    distance_range range = compute_distance_range(data_points, index_set, 1, index_set.count());
    const double difference = range.max - range.min;
    const double eps = 1e-12 * difference;
    const double interval = (range.max - range.min + eps) / num_bins;

    auto compute_bin = [interval, range] (double distance) -> size_t {
        return floor((distance - range.min) / interval);
    };

    // Loop over all interaction pairs
    for_each_pair(data_points, index_set, 1, index_set.count(),
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

    // Average out numbers
    for (size_t i = 0; i < num_bins; ++i) {
        data.gamma[i] /= (2 * data.num_pairs[i]);
        data.distance_averages[i] /= data.num_pairs[i];
    }

    return data;
}
