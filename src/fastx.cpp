#include <cstring>
#include <iostream>
#include "getopt.h"

#include "version.hpp"
#include "fastx_sample.hpp"


static
void Usage() {
    std::cerr << "fastx " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  fasta/fastq tool kit.\n"
              << std::endl;
    std::cerr << "Usage: fastx <command> <arguments>\n" << std::endl;
    std::cerr
            << "Commands:\n"
            << "  sample                      subsample sequences.\n"
            << std::endl;
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        Usage();
        return 0;
    }

    if ( strcmp(argv[1], "sample") == 0 )
    {
        FastxSampleMain(argc - 1, argv + 1);
    }

    return 0;
}
