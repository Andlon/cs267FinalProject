#pragma once

class parallel_options {
public:
    parallel_options(int processor_count, int desired_replication_factor);

    int replication_factor() const;
    int processor_count() const;

private:

};
