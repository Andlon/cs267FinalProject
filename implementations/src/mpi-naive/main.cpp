#include <iostream>
#include <string>
#include <algorithm>
#include "pair_index_set.h"
#include "data.h"

#include <fstream>

void print_gamma(const std::vector<double> & values,
                 const std::vector<size_t> & counts,
                 const std::vector<double> & dists,
                 size_t num_bins, const std::string & filename)
{
    std::ofstream ofile(filename.c_str());
    ofile << "#from \tto  \tn_pairs \tav_dist \tsemivariogram \n";
    for(int i = 0; i < num_bins; i++)
    {
        ofile << "?" << "\t" << "?" << "\t" << counts[i] << "\t" <<
                 dists[i] << "\t" << values[i] << "\n";
    }
    ofile.flush();
    std::cout << std::endl;
}

int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        std::cerr << "Missing input file argument. Aborting...\n";
        return 1;
    }

    // Hardcode these for now...
    const size_t num_bins = 15;
    const double min_dist = 18.45;
    const double max_dist = 20.000;
    const double interval = (max_dist - min_dist) / num_bins;

    std::vector<double> gamma(num_bins, 0);
    std::vector<double> average_distances(num_bins, 0);
    std::vector<size_t> num_pairs(num_bins, 0u);

    auto compute_bin = [interval] (double distance) -> size_t {
        return floor(distance / interval);
    };

    auto compute_gamma_contribution = [] (const data_point &p1, const data_point &p2) {
        return pow(p2.value - p1.value, 2);
    };

    auto compute_distance = [] (const data_point &p1, const data_point &p2) {
        return sqrt(pow((p1.x - p2.x), 2) + pow(p1.y - p2.y, 2));
    };

    const std::string input_path = argv[1];
    const std::vector<data_point> data_points = read_file_data(input_path);

    const pair_index_set index_set(data_points.size());

    // Loop over all interaction pairs
    for (index_type k = 1; k <= index_set.count(); ++k) {
        index_pair pair = index_set.map_to_pair(k);
        index_type i = pair.first;
        index_type j = pair.second;

        data_point p1 = data_points[i - 1];
        data_point p2 = data_points[j - 1];

        double distance = compute_distance(p1, p2);
        size_t bin = compute_bin(distance);

        gamma[bin] += compute_gamma_contribution(p1, p2);
        num_pairs[bin] += 1;
        average_distances[bin] += distance;
    }

    for (size_t i = 0; i < num_bins; ++i) {
        gamma[i] /= (2 * num_pairs[i]);
        average_distances[i] /= num_pairs[i];
    }

    print_gamma(gamma, num_pairs, average_distances, num_bins, "output.txt");

    return 0;
}
