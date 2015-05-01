#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <mpi.h>
#include <ios>
#include "data.h"
#include "variogram.h"

void print_all_points(const parallel_read_result & read_result, MPI_Comm communicator)
{
    int rank, ranks;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &ranks);

    // Sequentally print all points
    for (int i = 0; i < ranks; ++i)
    {
        if (rank == i)
            print_points(std::cout, read_result.data, true, read_result.start);
        MPI_Barrier(communicator);
    }
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
    parallel_read_result read_result = read_file_chunk_parallel(input_path);

    print_all_points(read_result, MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
