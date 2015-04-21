#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>

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
