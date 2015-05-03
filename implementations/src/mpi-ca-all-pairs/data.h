#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include <ostream>
#include "uniform_distribution.h"
#include <rect.h>

namespace pev {

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

pev::rect<double> bounding_rectangle(const std::vector<data_point> & data_points);

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

struct chunked_read_result
{
    chunked_read_result();

    // Holds the data points for this chunk
    std::vector<data_point> data;

    // Holds the total number of points contained in the file
    size_t global_point_count;

    // Optionally set to the maximum distance of the data points
    double max_distance;
};

/**
 * @brief read_file_chunk_parallel Reads a chunk of a file, in the sense that it reads
 * a certain amount of data points from the file, where each line holds a data point.
 * The result includes the local data and the global number of data points. Note that the
 * max_distance field is not set by this function.
 *
 * @param path The path to the file to read.
 * @param chunk_count The number of chunks to split the file into.
 * @param chunk_index The chunk index in { 0, ..., chunk_count - 1 } that this node should read.
 * @return A parallel_read_result struct that contains the aforementioned information. Note that
 * max_distance is not set.
 */
chunked_read_result read_file_chunk_parallel(const std::string & path, int chunk_count, int chunk_index);

}


