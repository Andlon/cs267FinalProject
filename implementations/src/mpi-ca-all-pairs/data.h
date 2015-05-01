#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include <ostream>
#include "uniform_distribution.h"

struct data_point {
    double x;
    double y;
    double value;
};

MPI_Datatype create_mpi_data_point_type();

/**
 * @brief distance Computes the distance between two data points.
 * @param p1
 * @param p2
 * @return The distance between the two points.
 */
inline double distance (const data_point &p1, const data_point &p2)
{
    return sqrt(pow((p1.x - p2.x), 2) + pow(p1.y - p2.y, 2));
}

/**
 * @brief print_points Prints the given data points to the given output stream.
 * @param out The output stream to print to.
 * @param points The data set from which to print points.
 * @param print_indices Toggles whether or not to print indices for the points
 */
void print_points(std::ostream & out, const std::vector<data_point> & points, bool print_indices = false, size_t start_index = 0);

/**
 * @brief read_file_data Reads points in <value> <x> <y> format from the specified file.
 * @param path The path to the file to read from.
 * @return A vector of data points contained in the file.
 */
std::vector<data_point> read_file_data(const std::string &path);

/**
 * @brief read_file_data_parallel Reads input file on root node and broadcasts to local nodes.
 * Returned vector contains all data points on all ranks.
 *
 * Note that since MPI is used behind the scenes, MPI must be initialized prior to calling the function.
 * @param path Path to input file.
 * @return Vector of data points.
 */
std::vector<data_point> read_file_data_parallel(const std::string &path, int root = 0);

struct parallel_read_result
{
    parallel_read_result();
    parallel_read_result(size_t point_count, size_t node_count);

    // Holds the data points local to this node
    std::vector<data_point> data;

    // Represents the distribution of data points (objects) across nodes (buckets)
    uniform_distribution distribution;

    // [start, end) is the interval of data point indices
    // in the larger set { 0, ..., n - 1 } that belongs
    // to this processor.
    size_t start;
    size_t end;
};

/**
 * @brief read_file_chunk_parallel Reads roughly equal sized chunks of a file distributed across
 * processors. The returned struct holds information about how many points each processor was
 * awarded, as well as the local interval of indices and the actual data points read for this
 * particular processor.
 *
 * You can specify a communicator to only split the file across the ranks in the given communicator.
 * Note that you should only call this function for the processors that are part of the given
 * communicator.
 *
 * @param path The path to the file to read.
 * @param communicator The communicator to split the file across.
 * @return A parallel_read_result struct that contains the aforementioned information.
 */
parallel_read_result read_file_chunk_parallel(const std::string & path, int chunk_count, int chunk_index);


