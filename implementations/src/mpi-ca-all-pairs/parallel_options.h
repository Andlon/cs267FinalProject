#pragma once

#include <cstdlib>
#include <stdexcept>
#include <iostream>

namespace pev {

/**
 * @brief Given a total count of processors and a desired replication factor,
 * the parallel_options ensures that
 *
 * P / c^2 = k for some integer k,
 *
 * where P is the number of active processors and c is the replication factor.
 * For predictable results, it is recommended to deliberately choose
 * processor_count and desired_replication_factor such that they
 * fulfill the above criterion.
 */
class parallel_options {
public:
    parallel_options(int processor_count, int desired_replication_factor);

    int replication_factor() const { return repl_factor; }
    int active_processor_count() const { return active_proc_count; }

private:
    int repl_factor;
    int active_proc_count;
};

inline parallel_options::parallel_options(int processor_count, int desired_replication_factor)
{
    if (processor_count < 1)
        throw new std::logic_error("processor_count must be greater than zero.");
    if (desired_replication_factor < 1)
        throw new std::logic_error("desired_replication_factor must be greater than zero.");

    int c = desired_replication_factor;
    while (processor_count / (c * c) == 0)
        --c;

    int k = processor_count / (c * c);

    // We want p / c^2 = k, where p, c and k are all integers.
    repl_factor = c;
    active_proc_count = (c * c) * k;
}

}
