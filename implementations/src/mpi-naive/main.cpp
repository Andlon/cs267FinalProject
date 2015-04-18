#include <iostream>
#include <string>
#include "pair_index_set.h"
#include "data.h"

int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        std::cerr << "Missing input file argument. Aborting...\n";
        return 1;
    }

    const std::string input_path = argv[1];
    const std::vector<data_point> data_points = read_file_data(input_path);

    std::cout << "Read " << data_points.size() << " points from input." << std::endl;
    return 0;
}
