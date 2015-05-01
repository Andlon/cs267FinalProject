#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <mpi.h>
#include <ios>
#include "data.h"
#include "variogram.h"

void print_gamma(const std::vector<double> & values,
                 const std::vector<size_t> & counts,
                 const std::vector<double> & dists,
                 size_t num_bins,
                 const std::string & filename)
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

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Compile-time constant for now
    const size_t NUM_BINS = 15;

    const std::string input_path = argv[1];
    parallel_read_result read_result = read_file_chunk_parallel(input_path);

    print_points(std::cout, read_result.data, true);

//    const variogram_data data = empirical_variogram_parallel(data_points, NUM_BINS);

//    if (rank == 0)
//        print_gamma(data.gamma, data.num_pairs, data.distance_averages, data.num_bins, "output.txt");

    MPI_Finalize();

    return 0;
}
