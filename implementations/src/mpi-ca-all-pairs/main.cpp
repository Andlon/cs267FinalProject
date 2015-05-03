#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <mpi.h>
#include <ios>
#include <iomanip>
#include "data.h"
#include "variogram.h"

void print_gamma(std::ostream & out,
                 const std::vector<double> & values,
                 const std::vector<size_t> & counts,
                 const std::vector<double> & dists,
                 size_t num_bins)
{
    using std::setw;
    using std::left;

    for(int i = 0; i < num_bins; i++)
    {
        out << left << setw(10) << counts[i] << "\t"
            << left << setw(10) << dists[i] << "\t"
            << left << setw(10) << values[i] << std::endl;
    }
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
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    // Compile-time constant for now
    const size_t NUM_BINS = 15;

    const std::string input_path = argv[1];

    parallel_options options(ranks, 16);

    if (rank == 0)
    {
        std::cout << "Using " << options.active_processor_count() << " of "
                  << ranks << " processors with replication factor " << options.replication_factor()
                  << std::endl;
    }

    variogram_data variogram = empirical_variogram_parallel(input_path, options, NUM_BINS);

    if (rank == 0)
    {
        print_gamma(std::cout,
                    variogram.gamma,
                    variogram.num_pairs,
                    variogram.distance_averages,
                    variogram.bin_count);
        int count = std::accumulate(variogram.num_pairs.begin(), variogram.num_pairs.end(),
                        0);
        std::cout << "Pair count: " << count << std::endl;
    }

    MPI_Finalize();

    return 0;
}
