#include "data.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <mpi.h>
#include <cassert>
#include <iomanip>

namespace {

data_point read_point(const std::string & line)
{
    std::istringstream stream(line);
    data_point point;
    stream >> point.value >> point.x >> point.y;
    return point;
}

std::streamsize file_size(const std::string & path)
{
    std::ifstream file(path, std::ios_base::ate | std::ios_base::binary);
    return file.tellg();
}

size_t file_line_width(const std::string & path)
{
    std::ifstream file(path);
    if (!file)
        throw std::ios_base::failure("Failed to open input file.");

    std::string line;
    std::getline(file, line);

    return line.size();
}

/**
 * @brief estimate_data_point_count Estimates the number of data points
 * in the file, assuming fixed line width. If the file does not have a fixed line width,
 * the returned value will be incorrect.
 *
 * The fixed line width is determined from the length of the first line.
 *
 * @param path
 * @return
 */
size_t estimate_data_point_count(const std::string & path)
{
    // Add one to account for the newline character.
    // This might cause problems on Windows platforms, as there you have
    // an extra carriage return character (\r), but our solution
    // is UNIX-only at this point.
    return file_size(path) / (file_line_width(path) + 1);
}

/**
 * @brief read_local_points Reads the points from the given file that are associated
 * with the given rank.
 * @param path The path to the file to read.
 * @param point_count The number of points in the file.
 * @param ranks The number of ranks the file should be distributed over.
 * @param rank The rank of the local processor.
 * @return A parallel_read_result containing the result.
 */
parallel_read_result read_local_points(const std::string & path, size_t point_count, int ranks, int rank) {
    parallel_read_result result(point_count, ranks);
    result.data.reserve(point_count);

    // Define local interval [a, b) of point indices to work on
    interval<size_t> interval = result.distribution.interval_for_bucket(rank);
    result.start = interval.a();
    result.end = interval.b();

    // Add 1 to account for newline. See comment in estimate_data_point_count()
    size_t line_width = file_line_width(path) + 1;

    // Read file data
    std::ifstream file(path);
    if (!file)
        throw std::ios_base::failure("Failed to open input file.");

    std::streampos startpos = interval.a() * line_width;

    file.seekg(startpos);
    std::string line;
    for (size_t i = result.start; i < result.end; ++i)
    {
        bool success = std::getline(file, line);
        if (success)
        {
            data_point point = read_point(line);
            result.data.push_back(point);
        } else {
            std::stringstream error_stream;
            error_stream << "Failed to read point " << i << " from the interval ["
                         << result.start << ", " << result.end << ")";
            throw new std::ios_base::failure(error_stream.str());
        }
    }

    assert(result.data.size() == result.distribution.objects_in_bucket(rank));

    return result;
}

}

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

parallel_read_result read_file_chunk_parallel(const std::string &path, int chunk_count, int chunk_index)
{
    size_t point_count = estimate_data_point_count(path);
    return read_local_points(path, point_count, chunk_count, chunk_index);
}


// TODO: Fix this function. Output looks really bad.
void print_points(std::ostream &out, const std::vector<data_point> &points, bool print_indices, size_t start_index)
{
    auto print = [&out] (auto object)
    {
        out << std::left << std::setw(10) << object << "\t";
    };

    for (size_t i = 0; i < points.size(); ++i)
    {
        const data_point & point = points[i];
        if (print_indices)
            print(i + start_index);

        print(point.value);
        print(point.x);
        print(point.y);
        out << std::endl;
    }
}


parallel_read_result::parallel_read_result()
    :   distribution(0, 1)
{

}

parallel_read_result::parallel_read_result(size_t point_count, size_t node_count)
    :   distribution(point_count, node_count)
{

}
