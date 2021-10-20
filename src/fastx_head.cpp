#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <string>
#include <cstdint>
#include <random>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <thread>

#include <getopt.h>
#include "zlib.h"
#include "htslib/kseq.h"
#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/process.hpp"
#include "common.hpp"
#include "version.hpp"



void FastxHeadBases(
    const std::string &ifilename1, const std::string &ifilename2,
    const std::string &ofilename1, const std::string &ofilename2,
    int64_t bases, const std::string &pigz)
{
    gzFile fp1 = gzopen(ifilename1.c_str(), "r");
    gzFile fp2 = gzopen(ifilename2.c_str(), "r");

    if (fp1 == nullptr)
    {
        std::perror(("Error! Can not open " + ifilename1).c_str());
        std::exit(1);
    }

    if (fp2 == nullptr)
    {
        std::perror(("Error! Can not open " + ifilename2).c_str());
        std::exit(1);
    }

    kseq_t *read1 = kseq_init(fp1);
    kseq_t *read2 = kseq_init(fp2);

    std::ostringstream pigz_command1;
    std::ostringstream pigz_command2;
    pigz_command1 << pigz << " - 1> " << ofilename1;
    pigz_command2 << pigz << " - 1> " << ofilename2; 

    FILE *stream1 = popen(pigz_command1.str().c_str(), "w");
    FILE *stream2 = popen(pigz_command2.str().c_str(), "w");

    int64_t base_count = 0;

    int ret1, ret2;
    while ((ret1 = kseq_read(read1) >= 0) &&
        (ret2 = kseq_read(read2) >= 0))
    {
        base_count += read1->seq.l;
        base_count += read2->seq.l;
        if (base_count <= bases) {
            std::string seq_str1 = kseqToStr(read1);
            std::string seq_str2 = kseqToStr(read2);
            fprintf(stream1, "%s\n", seq_str1.c_str());
            fprintf(stream2, "%s\n", seq_str2.c_str());
        } else {
            break;
        }
    }

    if (ret1 < -1)
    {
        std::cerr << "Error! Input fastq1 truncated! Last read name was "
            << read1->name.s << std::endl;
        std::exit(1);
    }

    if (ret2 < -1)
    {
        std::cerr << "Error! Input fastq2 truncated! Last read name was "
            << read2->name.s << std::endl;
        std::exit(1);
    }

    kseq_destroy(read1);
    kseq_destroy(read2);
    gzclose(fp1);
    gzclose(fp2);

    fflush(stream1);
    fflush(stream2);

    ret1 = pclose(stream1);
    ret2 = pclose(stream2);

    if (ret1 != 0)
    {
        std::cerr << "Error! Can not run pigz compression for read1 ! "
            << "Error code: " << ret1 << std::endl;
        std::exit(1);
    }

    if (ret2 != 0)
    {
        std::cerr << "Error! Can not run pigz compression for read2 ! "
            << "Error code: " << ret2 << std::endl;
        std::exit(1);
    }
}


void FastxHeadReads(
    const std::string &ifilename, const std::string &ofilename,
    int64_t reads, const std::string &pigz)
{
    gzFile fp = gzopen(ifilename.c_str(), "r");

    if (fp == nullptr)
    {
        std::perror(("Error! Can not open " + ifilename).c_str());
        std::exit(1);
    }

    kseq_t *read = kseq_init(fp);

    std::ostringstream pigz_command;
    pigz_command << pigz << " - 1> " << ofilename;

    FILE *stream = popen(pigz_command.str().c_str(), "w");

    int64_t read_count = 0;

    int ret;
    while ((ret = kseq_read(read)) >= 0)
    {
        ++read_count;
        if (read_count <= reads) {
            std::string seq_str = kseqToStr(read);
            fprintf(stream, "%s\n", seq_str.c_str());
        } else {
            break;
        }
    }

    if (ret < -1)
    {
        std::cerr << "Error! Input fastq truncated! File was "
            << ifilename << " Last read name was "
            << read->name.s  << std::endl;
        std::exit(1);
    }

    kseq_destroy(read);
    gzclose(fp);

    fflush(stream);

    ret = pclose(stream);

    if (ret != 0)
    {
        std::cerr << "Error! Can not run pigz compression! "
            << "Error code: " << ret << std::endl;
        std::exit(1);
    }
}


static
void Usage() {
    std::cerr << "fastx sample " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  subsample sequences.\n"
              << std::endl;
    std::cerr
            << "Options:\n"
            << "  -i, --in1, FILE             input fasta/fastq file name for read1.\n"
            << "  -I, --in2, FILE             input fasta/fastq file name for read2.\n"
            << "  -o, --out1, FILE            output fasta/fastq file name for read1.\n"
            << "  -O, --out2, FILE            output fasta/fastq file name for read2.\n"
            << "  -b, --bases, STR            get this value of bases(K/M/G).\n"
            << "  -n, --number, STR           get this value of read pairs(K/M/G).\n"
            << "  -p, --pigz, STR             path to pigz program.[pigz]\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -t, --thread, INT           number of threads for running pigz.[4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}


int FastxHeadMain(int argc, char **argv)
{
    if (argc == 1)
    {
        Usage();
        return 0;
    }

    static const struct option long_options[] = {
            {"in1", required_argument, 0, 'i'},
            {"in2", required_argument, 0, 'I'},
            {"out1", required_argument, 0, 'o'},
            {"out2", required_argument, 0, 'O'},
            {"bases", required_argument, 0, 'b'},
            {"number", required_argument, 0, 'n'},
            {"fraction", required_argument, 0, 'f'},
            {"level", required_argument, 0, 'l'},
            {"pigz", required_argument, 0, 'p'},
            {"seed", required_argument, 0, 's'},
            {"thread", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "i:I:o:O:b:n:p:l:t:hV";

    std::string input1;
    std::string input2;
    std::string output1;
    std::string output2;
    int64_t bases = -1;
    int64_t reads = -1;
    std::string pigz = "pigz";
    int compress_level = 6;
    int num_threads = 4;

    while ((c = getopt_long(
        argc, argv, opt_str, long_options, &long_idx)) != -1)
    {
        switch (c) {
            case 'i':
                input1 = optarg;
                break;
            case 'I':
                input2 = optarg;
                break;
            case 'o':
                output1 = optarg;
                break;
            case 'O':
                output2 = optarg;
                break;
            case 'b':
                bases = BasesStrToInt(optarg);
                if (bases <= 0)
                {
                    std::cerr << "Error! input bases must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'n':
                reads = BasesStrToInt(optarg);
                if (reads <= 0)
                {
                    std::cerr << "Error! input reads must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'p':
                pigz = optarg;
                break;
            case 'l':
                compress_level = boost::lexical_cast<int>(optarg);
                break;
            case 't':
                num_threads = boost::lexical_cast<int>(optarg);
                break;
            case 'h':
                Usage();
                return 0;
            case 'V':
                std::cerr << FASTX_VERSION << std::endl;
                return 0;
            default:
                Usage();
                return 1;
        }
    }

    if (bases < 0 && reads < 0)
    {
        std::cerr << "Error! must input bases(-b, --bases) "
            << "or number(-n, --number)." << std::endl;
        std::exit(1);
    }

    if (bases > 0 && reads > 0)
    {
        std::cerr << "Error! -b(--bases) and -n(--number) can not be "
            "used together!" << std::endl;
        std::exit(1);
    }

    std::ostringstream pigz_command;
    pigz_command << pigz << " -" << compress_level << " -p" << num_threads;

    if (bases > 0)
    {
        FastxHeadBases(input1, input2, output1, output2, bases,
            pigz_command.str());
    } else {
        std::thread th(FastxHeadReads,
            input1, output1, reads, pigz_command.str());
        FastxHeadReads(input2, output2, reads, pigz_command.str());
        th.join();
    }

    return 0;
}
