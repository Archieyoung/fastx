#include <boost/lexical_cast.hpp>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <zlib.h>
#include <getopt.h>
#include "fastx_split.hpp"
#include "common.hpp"
#include "htslib/bgzf.h"
#include "version.hpp"
#include "htslib/thread_pool.h"

// static
// int NumberOfDigits(int64_t number)
// {
//     int n = 1;
//     int ir = static_cast<double>(number) / 10.0;
//     while (ir >= 1) {
//         ir = ir / 10.0;
//         ++n;
//     }
//     return n;
// }

void FastxSplit(const std::string &ifilename, int64_t reads,
    int64_t bases, const std::string &prefix, const std::string &suffix,
    int threads, int compress_level)
{
    if (reads <= 0 && bases <= 0) return;
    if (reads > 0 && bases > 0) {
        std::cerr << "Error! Must input one of reads and bases "
            << "per output file not both!" << std::endl;
        std::exit(1);
    }

    bool split_by_reads = false;
    if (reads > 0)
    {
        split_by_reads = true;
    }

    gzFile fp = gzopen(ifilename.c_str(), "r");
    kseq_t *read = kseq_init(fp);

    int ret;
    int64_t read_count = 0;
    int64_t base_count = 0;
    bool open_new = false;
    int64_t n = 0;

    std::ostringstream ofilename;
    ofilename << prefix << "." << n << "." << suffix;
    hts_tpool *pool = hts_tpool_init(threads);
    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp = bgzf_open(ofilename.str().c_str(), mode_str.str().c_str());
    bgzf_thread_pool(bgzfp, pool, 0);

    while ((ret = kseq_read(read)) >= 0)
    {
        if (open_new) {
            if (split_by_reads)
            {
                ++read_count;
            } else {
                base_count += read->seq.l;
            }
            BgzfWriteKseq(bgzfp, read);
            if (split_by_reads && read_count >= reads)
            {
                open_new = true;
            }
            if (!split_by_reads && base_count >= bases)
            {
                open_new = true;
            }
        } else {
            ++n;
            ofilename.str("");   // clear
            bgzf_close(bgzfp);
            ofilename << prefix << "." << n << "." << suffix;
            bgzfp = bgzf_open(ofilename.str().c_str(), mode_str.str().c_str());
            bgzf_thread_pool(bgzfp, pool, 0);
            base_count = 0;
            read_count = 0;
            if (split_by_reads)
            {
                ++read_count;
            } else {
                base_count += read->seq.l;
            }
            BgzfWriteKseq(bgzfp, read);
        }
    }

    bgzf_close(bgzfp);
    hts_tpool_destroy(pool);
    kseq_destroy(read);
    gzclose(fp);
}

static
void Usage() {
    std::cerr << "fastx split " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  split fasta/q file.\n"
              << std::endl;
    std::cerr
            << "Options:\n"
            << "  -i, --in1, FILE             input fasta/fastq file name for read1.\n"
            << "  -I, --in2, FILE             input fasta/fastq file name for read2(optional).\n"
            << "  -o, --out1, FILE            output fasta/fastq file name for read1.\n"
            << "  -O, --out2, FILE            output fasta/fastq file name for read2(with -I).\n"
            << "  -b, --bases, STR            put this value of bases per output file(K/M/G).\n"
            << "  -r, --reads, STR            put this value of reads per output file(K/M/G).\n"
            << "  -n, --number, int           generate this value of output files"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -t, --thread, INT           number of threads for running pigz.[4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}


int FastxSplit(int argc, char **argv)
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
            {"reads", required_argument, 0, 'r'},
            {"number", required_argument, 0, 'n'},
            {"level", required_argument, 0, 'l'},
            {"thread", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "i:I:o:O:b:f:p:l:s:t:hV";

    std::string input1 = "";
    std::string input2 = "";
    std::string output1 = "";
    std::string output2 = "";
    int64_t bases = -1;
    int64_t reads = -1;
    int number = -1;
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
            case 'r':
                reads = BasesStrToInt(optarg);
                if (reads <= 0)
                {
                    std::cerr << "Error! input reads must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'n':
                number = boost::lexical_cast<int>(opt_str);
                if (number <= 0)
                {
                    std::cerr << "Error! input number must be positive!"
                        << std::endl;
                    std::exit(1);
                }
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

    if (input1.empty())
    {
        std::cerr << "Error! must input the first fasta/q file." << std::endl;
        std::exit(1);
    }

    if (bases < 0 && reads < 0 && number < 0)
    {
        std::cerr << "Error! must input bases(-b, --bases) or "
            << "reads(-r, --reads) or number(-n, --number)"
            << std::endl;
        std::exit(1);
    }

    // single end
    if (!input1.empty() && input2.empty())
    {
        if (number > 0)
        {
            int64_t total_reads;
            int64_t total_bases;
            FastxCount(input1, total_reads, total_bases);
            bases = static_cast<double>(total_bases) / number;
            FastxSplit(input1, -1, bases, "T1", "fastq.gz", num_threads, compress_level);
        }
    }

    return 0;
}


