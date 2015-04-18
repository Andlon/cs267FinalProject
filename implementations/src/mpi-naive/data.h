#pragma once

#include <vector>
#include <string>

struct data_point {
    double x;
    double y;
    double value;
};

/**
 * @brief read_file_data Reads points in <value> <x> <y> format from the specified file.
 * @param path The path to the file to read from.
 * @return A vector of data points contained in the file.
 */
std::vector<data_point> read_file_data(const std::string &path);
