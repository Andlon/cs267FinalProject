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

    const pair_index_set index_set(data_points.size());

    // Loop over all interaction pairs
    for (index_type k = 1; k <= index_set.count(); ++k) {
        index_pair pair = index_set.map_to_pair(k);
        index_type i = pair.first;
        index_type j = pair.second;

        std::cout << "k = " << k << "- (" << i << ", " << j << ")" << std::endl;
    }

    return 0;
}
