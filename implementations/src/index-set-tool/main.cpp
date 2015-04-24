#include <iostream>
#include <iomanip>
#include <pair_index_set.h>
#include <cmath>
#include <sstream>
#include <string>

unsigned int count_digits_base10(u_int64_t max_number)
{
    // Note that this will only work up to a certain precision
    return static_cast<u_int64_t>(floor( log10( max_number ) )) + 1;
}

void explain_usage()
{
    std::cout << "Given a number of points and the number of processors, prints triplets (k, i, j) "
                 "which relate an indices k in an index set { 1, ..., K} to unique pairs (i, j)."
              << std::endl
              << "Usage: " << std::endl
              << "index-set-tool <num-points>" << std::endl
              << std::endl;
}

int handle_invalid_point_count()
{
    std::cerr << "ERROR: Invalid number of points supplied. Must be an integer greater than zero." << std::endl;
    return -1;
}

std::string string_from_pair(index_type i, index_type j)
{
    std::stringstream ss;
    ss << "(" << i << ", " << j << ")";
    return ss.str();
}

void print_pairs(size_t point_count)
{
    pair_index_set index_set(point_count);
    auto interaction_count = index_set.count();

    // Determine column_width based on what takes more space
    // to write out of the largest k and the pairs (i, j)
    size_t column_width = 2 + std::max(
                count_digits_base10(interaction_count),
                4 + 2 * count_digits_base10(point_count));

    std::streamsize old_width = std::cout.width();

    // Lambda function to return a string of column width given data to contain
    auto col = [column_width] (auto data) {
        std::stringstream ss;
        ss << std::left << std::setw(column_width) << data;
        return ss.str();
    };

    // Print header
    std::cout << "N = " << point_count << " points" << std::endl
              << "K = " << interaction_count << " interactions" << std::endl
              << col("k") << col("(i, j)") << std::endl;

    // Print separator
    std::cout << std::string(2 * column_width, '-') << std::endl;

    // Print indices and pairs
    for (index_type k = 1; k <= interaction_count; ++k)
    {
        auto pair = index_set.map_to_pair(k);
        auto i = pair.first;
        auto j = pair.second;

        std::cout << col(k) << col(string_from_pair(i, j)) << std::endl;
    }

    std::cout.width(old_width);
}

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        explain_usage();
    }
    else
    {
        size_t point_count = strtoul(argv[1], 0, 0);

        if (point_count == 0u)
            return handle_invalid_point_count();

        print_pairs(point_count);
    }

    return 0;
}
