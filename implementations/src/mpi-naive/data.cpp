#include "data.h"
#include <fstream>
#include <sstream>
#include <mpi.h>

std::vector<data_point> read_file_data(const std::string &path)
{
    std::vector<data_point> points;
    std::ifstream file(path.c_str());

    if (!file)
        throw std::ios_base::failure("Failed to open input file.");

    for (std::string line; std::getline(file, line);)
    {
        std::istringstream stream(line);
        data_point point;
        stream >> point.value >> point.x >> point.y;
        points.push_back(point);
    }

    // Note that the move constructor makes returning
    // the potentially large vector a cheap operation
    return points;
}


std::vector<data_point> read_file_data_parallel(const std::string &path, int root)
{
    int rank;
    int num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    std::vector<data_point> data_points;

    // Only root node reads file data directly
    if (rank == root)
        data_points = read_file_data(path);

    // Broadcast number of data points to all ranks
    u_int64_t num_points = data_points.size();
    MPI_Bcast(&num_points, 1, MPI_UINT64_T, root, MPI_COMM_WORLD);

    MPI_Datatype type = create_mpi_data_point_type();

    // Allocate space for points and broadcast data points from root to all ranks
    data_points.resize(num_points);
    MPI_Bcast(data_points.data(), num_points, type, root, MPI_COMM_WORLD);
    return data_points;
}


MPI_Datatype create_mpi_data_point_type()
{
    MPI_Datatype type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);
    return type;
}
