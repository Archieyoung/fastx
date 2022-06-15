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
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"
#include "zlib.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/process.hpp"
#include "common.hpp"
#include "version.hpp"


struct SubsampleSummary {
    int64_t total_bases;
    int64_t expected_subsample_bases;
    int64_t real_subsample_bases;
    double expected_subsample_fraction;
    double real_subsample_fraction;
};


int BgzfWriteKseq(BGZF *fp, const kseq_t *seq) {
    int ret;
    ret = bgzf_write(fp, (seq->qual.l ? "@" : ">"), 1);
    if (ret < 0) return ret;
    ret = bgzf_write(fp, seq->name.s, seq->name.l);
    if (ret < 0) return ret;
    if (seq->comment.l) {
        ret = bgzf_write(fp, " ", 1);
        if (ret < 0) return ret;
        ret = bgzf_write(fp, seq->comment.s, seq->comment.l);
        if (ret < 0) return ret;
    }
    ret = bgzf_write(fp, "\n", 1);
    if (ret < 0) return ret;
    ret = bgzf_write(fp, seq->seq.s, seq->seq.l);
    if (ret < 0) return ret;
    ret = bgzf_write(fp, "\n", 1);
    if (ret < 0) return ret;
    if (seq->qual.l) {
        ret = bgzf_write(fp, "+\n", 2);
        if (ret < 0) return ret;
        ret = bgzf_write(fp, seq->qual.s, seq->qual.l);
        if (ret < 0) return ret;
    }
    ret = bgzf_write(fp, "\n", 1);
    return ret;
}


void FastxSample(
    const std::string &ifilename1, const std::string &ifilename2,
    const std::string &ofilename1, const std::string &ofilename2,
    double fraction, int64_t bases, double mean_length, int seed,
    int compress_level, int threads, SubsampleSummary &summary)
{
    if (fraction >= 1.0)
    {
        boost::filesystem::copy_file(ifilename1, ofilename1,
            boost::filesystem::copy_option::overwrite_if_exists);
        boost::filesystem::copy_file(ifilename2, ofilename2,
            boost::filesystem::copy_option::overwrite_if_exists);
        summary.real_subsample_bases = summary.total_bases;
        summary.real_subsample_fraction = 1.0;
        return;
    }

    // subsample

    // scale fraction to 'avoid' subsample less bases than expected, 
    // side effect is that reads in the front of the files are more likely to be 
    // sampled than reads in the tail of the files
    fraction *= 1.05;

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

    std::random_device rd;
    std::mt19937 g(rd());
    g.seed(seed);
    std::uniform_real_distribution<double> random_u(0.0, 1.0);

    hts_tpool *pool = hts_tpool_init(threads);
    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp1 = bgzf_open(ofilename1.c_str(), mode_str.str().c_str());
    bgzf_thread_pool(bgzfp1, pool, 0);
    BGZF* bgzfp2 = bgzf_open(ofilename2.c_str(), mode_str.str().c_str());
    bgzf_thread_pool(bgzfp2, pool, 0);

    int64_t total_bases = 0;
    int64_t subsample_bases = 0;

    int ret1, ret2, ret3;
    while ((ret1 = kseq_read(read1) >= 0) &&
        (ret2 = kseq_read(read2) >= 0))
    {
        total_bases += read1->seq.l;
        total_bases += read2->seq.l;
        // weight p by read length
        double p = random_u(g) * (read1->seq.l + read2->seq.l) / mean_length;
        if (p <= fraction && subsample_bases < bases) {
            subsample_bases += read1->seq.l;
            subsample_bases += read2->seq.l;
            ret3 = BgzfWriteKseq(bgzfp1, read1);
            if (ret3 < 0) {
                std::cerr << "Error! Failed to write read1: "
                    << read1->name.s << std::endl;
                std::exit(1);
            }
            ret3 = BgzfWriteKseq(bgzfp2, read2);
            if (ret3 < 0) {
                std::cerr << "Error! Failed to write read2: "
                    << read2->name.s << std::endl;
                std::exit(1);
            }
        }
        if (subsample_bases >= bases)
        {
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
    bgzf_close(bgzfp1);
    bgzf_close(bgzfp2);
    hts_tpool_destroy(pool);
    gzclose(fp1);
    gzclose(fp2);

    summary.real_subsample_bases = subsample_bases;
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
            << "  -b, --bases, STR            expected bases to subsample(K/M/G).\n"
            << "  -f, --fraction, FLOAT       expected fraction of bases to subsample.\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -s, --seed, INT             random seed.[11]\n"
            << "  -t, --thread, INT           number of threads.[4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}


int FastxSampleMain(int argc, char **argv)
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
            {"fraction", required_argument, 0, 'f'},
            {"level", required_argument, 0, 'l'},
            {"seed", required_argument, 0, 's'},
            {"thread", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "i:I:o:O:b:f:p:l:s:t:hV";

    std::string input1;
    std::string input2;
    std::string output1;
    std::string output2;
    int64_t bases = -1;
    double fraction = -1.0;
    int compress_level = 6;
    int seed = 11;
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
            case 'f':
                fraction = boost::lexical_cast<double>(optarg);
                if (fraction <= 0)
                {
                    std::cerr << "Error! input fraction must be positive!"
                        << std::endl;
                    std::exit(1);
                }
                break;
            case 'l':
                compress_level = boost::lexical_cast<int>(optarg);
                break;
            case 's':
                seed = boost::lexical_cast<int>(optarg);
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

    if (bases < 0 && fraction < 0)
    {
        std::cerr << "Error! must input expected bases(-b, --bases) "
            << "or expected fraction(-f, --fraction)." << std::endl;
        std::exit(1);
    }

    if (bases > 0 && fraction > 0)
    {
        std::cerr << "Error! -b(--bases) and -f(--fraction) can not be "
            "used together!" << std::endl;
        std::exit(1);
    }

    SubsampleSummary summary;

    int64_t total_reads, total_bases;
    FastxCountPair(input1, input2, total_reads, total_bases);

    summary.total_bases = total_bases;

    // mean length of read1 + read2
    double mean_length = static_cast<double>(total_bases) / total_reads * 2.0;

    if (fraction < 0)
    {
        fraction = static_cast<double>(bases) / total_bases;
    }

    summary.expected_subsample_fraction = fraction;

    if (bases < 0)
    {
        bases = static_cast<int64_t>(std::round(fraction * total_bases));
    }

    FastxSample(input1, input2, output1, output2,
        fraction, bases, mean_length, seed, compress_level,
        num_threads, summary);

    summary.expected_subsample_bases = bases;

    summary.real_subsample_fraction =
        static_cast<double>(summary.real_subsample_bases) / summary.total_bases;
    
    std::cerr << "#Subsample Summary" << std::endl;
    std::cerr << "Total bases: "
        << summary.total_bases << std::endl;
    std::cerr << "Expected bases: "
        << summary.expected_subsample_bases << std::endl;
    std::cerr << "Real bases: "
        << summary.real_subsample_bases << std::endl;
    std::cerr << "Expected fraction: "
        << summary.expected_subsample_fraction << std::endl;
    std::cerr << "Real fraction: "
        << summary.real_subsample_fraction << std::endl;

    return 0;
}
