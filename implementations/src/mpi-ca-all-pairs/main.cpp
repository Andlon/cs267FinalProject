#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <mpi.h>
#include <ios>
#include <iomanip>
#include <tclap/CmdLine.h>
#include "data.h"
#include "variogram.h"

void print_gamma(std::ostream & out,
                 const std::vector<double> & values,
                 const std::vector<size_t> & counts,
                 const std::vector<double> & dists,
                 size_t num_bins)
{
    using std::setw;
    using std::left;

    for(int i = 0; i < num_bins; i++)
    {
        out << left << setw(10) << counts[i] << "\t"
            << left << setw(10) << dists[i] << "\t"
            << left << setw(10) << values[i] << std::endl;
    }
    std::cout << std::endl;
}

struct positive_constraint : TCLAP::Constraint<int> {
    std::string description() const { return "Verifies that the number is positive."; }
    std::string shortID() const { return "positive integer"; }
    bool check(const int &value) const { return value > 0; }
} positive_constraint;

int main(int argc, char ** argv)
{
    try {
        const std::string description =
                "Parallel Empirical Variogram: A massively parallel solution "
                "for computing empirical variograms of large-scale data.";

        const std::string input_desc = "Path to input file";
        const std::string output_desc = "Path to output file. Defaults to stdout.";
        const std::string timing_desc = "Path to timing file. "
                                        "Defaults to not writing timing information.";
        const std::string repl_desc = "Desired replication factor. Defaults to 1.";
        const std::string bins_desc = "Number of distance bins to use for empirical variogram. "
                                      "Defaults to 15.";

        // TODO: Improve input file descriptions to describe output formats
        TCLAP::CmdLine cmd(description, ' ', "0.1");
        TCLAP::ValueArg<std::string> input_file_arg("i", "input", input_desc, true, "", "path");
        TCLAP::ValueArg<std::string> output_file_arg("o", "output", output_desc, false, "", "path");
        TCLAP::ValueArg<std::string> timing_file_arg("t", "timing", timing_desc, false, "", "path");
        TCLAP::ValueArg<int> repl_factor_arg("c", "replication", repl_desc,
                                             false, 1, &positive_constraint);
        TCLAP::ValueArg<int> num_bins_arg("b", "bins", bins_desc, false, 15, &positive_constraint);
        cmd.add(input_file_arg);
        cmd.add(output_file_arg);
        cmd.add(timing_file_arg);
        cmd.add(repl_factor_arg);
        cmd.add(num_bins_arg);
        cmd.parse(argc, argv);

        const int num_bins = num_bins_arg.getValue();
        const int replication_factor = repl_factor_arg.getValue();
        const std::string input_path = input_file_arg.getValue();
        const std::string output_path = output_file_arg.getValue();

        MPI_Init(&argc, &argv);
        int rank, ranks;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);

        pev::parallel_options options(ranks, replication_factor);
        if (rank == 0)
        {
            std::cout << "Using " << options.active_processor_count() << " of "
                      << ranks << " processors with replication factor " << options.replication_factor()
                      << std::endl;
        }

        pev::variogram_data variogram = pev::empirical_variogram_parallel(input_path, options, num_bins);

        if (rank == 0)
        {
            std::ofstream outfile;
            if (output_file_arg.isSet())
                outfile.open(output_path);

            // Print to stdout if output file not supplied
            std::ostream & out = output_file_arg.isSet()
                    ? outfile
                    : std::cout;

            print_gamma(out,
                        variogram.gamma,
                        variogram.num_pairs,
                        variogram.distance_averages,
                        variogram.bin_count);
        }

        MPI_Finalize();

        return 0;

    } catch (const TCLAP::ArgException & e)
    {
        std::cerr << "ERROR: " << e.what() << " for argument " << e.argId() << std::endl;
        return -1;
    }
}
