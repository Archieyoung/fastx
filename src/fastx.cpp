#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include "fastx_subseq.hpp"
#include "getopt.h"

#include "version.hpp"
#include "fastx_sample.hpp"
#include "fastx_head.hpp"
#include "fastx_split.hpp"


static
void Usage() {
    std::cerr << "fastx " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  fasta/fastq tool kit.\n"
              << std::endl;
    std::cerr << "Usage: fastx <command> <arguments>\n" << std::endl;
    std::cerr
            << "Commands:\n"
            << "  head           head sequences\n"
            << "  sample         subsample sequences\n"
            << "  split          split fasta/fastq files.\n"
            << "  subseq         extract subsequences of fasta/fastq"
            << std::endl;
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        Usage();
        return 0;
    }

    std::map<std::string, bool> registered_commands = {
        {"head", true},
        {"sample", true},
        {"split", true},
        {"subseq", true}
        };
    
    if (registered_commands.find(argv[1]) != registered_commands.end())
    {
        if(!registered_commands.at(argv[1])) {
            std::cerr << "Error! " << argv[1] << " is deprecated!" << std::endl;
            std::exit(1);
        }
    } else {
        std::cerr << "Error! Unrecognized commond " << argv[1] << std::endl;
        std::exit(1);
    }

    if ( strcmp(argv[1], "head") == 0 ) {
        return FastxHeadMain(argc - 1, argv + 1);
    } else if ( strcmp(argv[1], "sample") == 0 )
    {
        return FastxSampleMain(argc - 1, argv + 1);
    } else if ( strcmp(argv[1], "split") == 0 )
    {
        return FastxSplitMain(argc - 1, argv + 1);
    } else if ( strcmp(argv[1], "subseq") == 0 )
    {
        return FastxSubseqMain(argc - 1, argv + 1);
    }

    return 0;
}
