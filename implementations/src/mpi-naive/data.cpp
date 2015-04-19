#include "data.h"
#include <fstream>
#include <sstream>

std::vector<data_point> read_file_data(const std::string &path)
{
    std::vector<data_point> points;
    std::ifstream file(path.c_str());

    if (!file)
        throw std::ios_base::failure("Failed to open input file");

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
