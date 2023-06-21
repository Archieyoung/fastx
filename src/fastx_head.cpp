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
#include "htslib/thread_pool.h"
#include "utils.hpp"
#include "kseq_utils.hpp"
#include "seq_reader.hpp"
#include "version.hpp"


void FastxHeadBasesSingle(
    const std::string &ifilename1, const std::string &ofilename1,
    int64_t bases, int threads, int compress_level)
{
    gzFile fp1 = gzopen(ifilename1.c_str(), "r");

    if (fp1 == nullptr)
    {
        std::perror(("Error! Can not open " + ifilename1).c_str());
        std::exit(1);
    }

    kseq_t *read1 = kseq_init(fp1);

    hts_tpool *pool = hts_tpool_init(threads);
    if (pool == NULL) {
        std::cerr << "Error! hts_tpool_init can not init thread pool "
            << std::endl;
        std::exit(1);
    }
    
    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp1 = bgzf_open(ofilename1.c_str(), mode_str.str().c_str());
    if (bgzfp1 == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename1 << " for writing" << std::endl;
        std::exit(1);
    }

    bgzf_thread_pool(bgzfp1, pool, 0);

    int64_t base_count = 0;

    int ret1, ret2;
    while ((ret1 = kseq_read(read1)) >= 0)
    {
        base_count += read1->seq.l;
        if (base_count <= bases) {
            ret2 = BgzfWriteKseq(bgzfp1, read1);
            if (ret2 < 0) {
                std::cerr << "Error! Failed to write read: "
                    << read1->name.s << std::endl;
                std::exit(1);
            }
        } else {
            break;
        }
    }

    if (ret1 < -1)
    {
        std::cerr << "Error! Input fastq truncated! Last read name was "
            << read1->name.s << std::endl;
        std::exit(1);
    }

    kseq_destroy(read1);
    bgzf_close(bgzfp1);
    hts_tpool_destroy(pool);
    gzclose(fp1);
}


void FastxHeadBasesPair(
    const std::string &ifilename1, const std::string &ifilename2,
    const std::string &ofilename1, const std::string &ofilename2,
    int64_t bases, int threads, int compress_level)
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

    SeqReader reader1 = SeqReader(fp1);
    SeqReader reader2 = SeqReader(fp2);

    hts_tpool *pool = hts_tpool_init(threads);
    if (pool == NULL) {
        std::cerr << "Error! hts_tpool_init can not init thread pool "
            << std::endl;
        std::exit(1);
    }

    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp1 = bgzf_open(ofilename1.c_str(), mode_str.str().c_str());
    if (bgzfp1 == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename1 << " for writing" << std::endl;
        std::exit(1);
    }

    bgzf_thread_pool(bgzfp1, pool, 0);
    BGZF* bgzfp2 = bgzf_open(ofilename2.c_str(), mode_str.str().c_str());
    if (bgzfp2 == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename2 << " for writing" << std::endl;
        std::exit(1);
    }

    bgzf_thread_pool(bgzfp2, pool, 0);

    int64_t base_count = 0;

    int ret;
    kseq_t *read1 = nullptr;
    kseq_t *read2 = nullptr;
    while ((read1 = reader1.read()) != nullptr &&
        (read2 = reader2.read()) != nullptr)
    {
        base_count += read1->seq.l;
        base_count += read2->seq.l;
        if (base_count <= bases) {
            ret = BgzfWriteKseq(bgzfp1, read1);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read1: "
                    << read1->name.s << std::endl;
                std::exit(1);
            }
            ret = BgzfWriteKseq(bgzfp2, read2);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read2: "
                    << read2->name.s << std::endl;
                std::exit(1);
            }
        } else {
            break;
        }
    }

    bgzf_close(bgzfp1);
    bgzf_close(bgzfp2);
    hts_tpool_destroy(pool);
    gzclose(fp1);
    gzclose(fp2);
}


void FastxHeadReadsSingle(
    const std::string &ifilename, const std::string &ofilename,
    int64_t reads, int threads, int compress_level)
{
    gzFile fp = gzopen(ifilename.c_str(), "r");

    if (fp == nullptr)
    {
        std::perror(("Error! Can not open " + ifilename).c_str());
        std::exit(1);
    }

    kseq_t *read = kseq_init(fp);

    hts_tpool *pool = hts_tpool_init(threads);
    if (pool == NULL) {
        std::cerr << "Error! hts_tpool_init can not init thread pool "
            << std::endl;
        std::exit(1);
    }

    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp = bgzf_open(ofilename.c_str(), mode_str.str().c_str());
    if (bgzfp == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename << " for writing" << std::endl;
        std::exit(1);
    }

    bgzf_thread_pool(bgzfp, pool, 0);

    int64_t read_count = 0;

    int ret, ret1;
    while ((ret = kseq_read(read)) >= 0)
    {
        ++read_count;
        if (read_count <= reads) {
            ret1 = BgzfWriteKseq(bgzfp, read);
            if (ret1 < 0) {
                std::cerr << "Error! Failed to write read2: "
                    << read->name.s << std::endl;
                std::exit(1);
            }
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
    bgzf_close(bgzfp);
    hts_tpool_destroy(pool);
    gzclose(fp);
}


static
void Usage() {
    std::cerr << "fastx head " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  get sequences from front of fasta/q files.\n"
              << std::endl;
    std::cerr
            << "Options:\n"
            << "  -i, --in1, FILE             input fasta/fastq file name for read1.\n"
            << "  -I, --in2, FILE             input fasta/fastq file name for read2.\n"
            << "  -o, --out1, FILE            output fasta/fastq file name for read1.\n"
            << "  -O, --out2, FILE            output fasta/fastq file name for read2.\n"
            << "  -b, --bases, STR            get this value of bases(K/M/G).\n"
            << "  -n, --number, STR           get this value of read pairs(K/M/G).\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -t, --thread, INT           number of threads.[4]\n"
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
                bases = KmgStrToInt(optarg);
                if (bases <= 0)
                {
                    std::cerr << "Error! input bases must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'n':
                reads = KmgStrToInt(optarg);
                if (reads <= 0)
                {
                    std::cerr << "Error! input reads must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'l':
                compress_level = SafeStrtol(optarg, 10);
                break;
            case 't':
                num_threads = SafeStrtol(optarg, 10);
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

    if (input1.empty()) {
        std::cerr << "Error! Must set at least one input fasta/fastq file "
            << "using -i(--in1)." << std::endl;
        std::exit(1);
    }

    if (output1.empty()) {
        std::cerr << "Error! Must set at least one output fasta/fastq file "
            << "using -o(--out1)." << std::endl;
        std::exit(1);
    }

    if (!input2.empty() && output2.empty()) {
        std::cerr << "Error! Must set at the second output fasta/fastq file "
            << "using -O(--out2) When inputting 2 fasta/fastq files."
            << std::endl;
        std::exit(1);
    }

    if (bases < 0 && reads < 0)
    {
        std::cerr << "Error! Must input bases(-b, --bases) "
            << "or number(-n, --number)." << std::endl;
        std::exit(1);
    }

    if (bases > 0 && reads > 0)
    {
        std::cerr << "Error! -b(--bases) and -n(--number) can not be "
            "used together!" << std::endl;
        std::exit(1);
    }

    if (compress_level < 0) {
        std::cerr << "Error! Compression level must be greater than or equal to"
            << " 0" << std::endl;
        std::exit(1);
    }

    if (num_threads < 1) {
        std::cerr << "Error! Number of threads -t(--threads) must greater"
            << " than 0" << std::endl;
        std::exit(1);
    }

    if (input2.empty()) {
        // single read
        if (bases > 0)
        {
            FastxHeadBasesSingle(
                input1, output1, bases, num_threads, compress_level);
        } else {
            FastxHeadReadsSingle(
                input1, output1, reads, num_threads, compress_level);
        }
    } else {
        // paired reads
        if (bases > 0)
        {
            FastxHeadBasesPair(input1, input2, output1, output2, bases,
                num_threads, compress_level);
        } else {
            FastxHeadReadsSingle(
                input1, output1, reads, num_threads, compress_level);
            FastxHeadReadsSingle(
                input2, output2, reads, num_threads, compress_level);
        }
    }

    return 0;
}
