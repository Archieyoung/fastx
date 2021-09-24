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
#include "version.hpp"

KSEQ_INIT(gzFile, gzread)


void FastxCount(const std::string &filename, int64_t &reads, int64_t &bases)
{
    reads = 0;
    bases = 0;
    gzFile fp = gzopen(filename.c_str(), "r");

    kseq_t *read = kseq_init(fp);

    int ret;

    while ((ret = kseq_read(read) >= 0))
    {
        reads += 1;
        bases += read->seq.l;
    }

    if (ret < -1)
    {
        std::cerr << "Error! Input fastq truncated! Last read name was "
            << read->name.s << std::endl;
        std::exit(1);
    }
    
    kseq_destroy(read);
    gzclose(fp);
}


void FastxCountPair(const std::string &ifilename1,
    const std::string &ifilename2, int64_t &reads, int64_t &bases)
{
    int64_t read1_counts, read1_bases;
    int64_t read2_counts, read2_bases;

    std::thread th(FastxCount, ifilename1,
        std::ref(read1_counts), std::ref(read1_bases));

    FastxCount(ifilename2, read2_counts, read2_bases);
    
    th.join();

    reads = read1_counts + read2_counts;
    bases = read1_bases + read2_bases; 
}


static inline std::string kseqToStr(const kseq_t *seq) {
    std::ostringstream seq_str;
    seq_str << (seq->qual.l ? "@" : ">");
    seq_str << seq->name.s;
    if (seq->comment.l) {
        seq_str << " ";
        seq_str << seq->comment.s;
    }
    seq_str << "\n";
    seq_str << seq->seq.s;
    seq_str << "\n";
    if (seq->qual.l) {
        seq_str << "+\n";
        seq_str << seq->qual.s;
    }
    return seq_str.str();
}


struct SubsampleSummary {
    int64_t total_bases;
    int64_t expected_subsample_bases;
    int64_t real_subsample_bases;
    double expected_subsample_fraction;
    double real_subsample_fraction;
};


void FastxSample(
    const std::string &ifilename1, const std::string &ifilename2,
    const std::string &ofilename1, const std::string &ofilename2,
    double fraction, int64_t bases, double mean_length, int seed,
    const std::string &pigz, SubsampleSummary &summary)
{
    if (fraction >= 1.0)
    {
        boost::filesystem::copy_file(ifilename1, ofilename1,
            boost::filesystem::copy_option::overwrite_if_exists);
        boost::filesystem::copy_file(ifilename2, ofilename2,
            boost::filesystem::copy_option::overwrite_if_exists);
        summary.real_subsample_fraction = 1.0;
        return;
    }

    // subsample
    fraction *= 1.05;

    gzFile fp1 = gzopen(ifilename1.c_str(), "r");
    gzFile fp2 = gzopen(ifilename2.c_str(), "r");

    kseq_t *read1 = kseq_init(fp1);
    kseq_t *read2 = kseq_init(fp2);

    std::random_device rd;
    std::mt19937 g(rd());
    g.seed(seed);
    std::uniform_real_distribution<double> random_u(0.0, 1.0);

    std::ostringstream pigz_command1;
    std::ostringstream pigz_command2;
    pigz_command1 << pigz << " - 1> " << ofilename1;
    pigz_command2 << pigz << " - 1> " << ofilename2; 

    FILE *stream1 = popen(pigz_command1.str().c_str(), "w");
    FILE *stream2 = popen(pigz_command2.str().c_str(), "w");


    int64_t total_bases = 0;
    int64_t subsample_bases = 0;

    int ret1, ret2;
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
            std::string seq_str1 = kseqToStr(read1);
            std::string seq_str2 = kseqToStr(read2);
            fprintf(stream1, "%s\n", seq_str1.c_str());
            fprintf(stream2, "%s\n", seq_str2.c_str());
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

    summary.total_bases = total_bases;
    summary.real_subsample_bases = subsample_bases;
    summary.real_subsample_fraction =
        static_cast<double>(subsample_bases) / total_bases;
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
            << "  -p, --pigz, STR             path to pigz program.[pigz]\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -s, --seed, INT             random seed.[11]\n"
            << "  -t, --thread, INT           number of threads for running pigz.[4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}

static
int64_t BasesStrToInt(const std::string &bases_str) {
    if (bases_str.empty())
    {
        std::cerr << "Error! Input bases string is empty!" << std::endl;
        std::exit(1);
    }

    int64_t bases = 0;
    switch (bases_str.back()) {
        case 'G':
        case 'g':
            bases = 1000000000 * boost::lexical_cast<int64_t>(
                bases_str.substr(0, bases_str.size() - 1));
            break;
        case 'M':
        case 'm':
            bases = 1000000 * boost::lexical_cast<int64_t>(
                bases_str.substr(0, bases_str.size() - 1));
            break;
        case 'K':
        case 'k':
            bases = 1000 * boost::lexical_cast<int64_t>(
                bases_str.substr(0, bases_str.size() - 1));
            break;
        default:
            bases = boost::lexical_cast<int64_t>(
                bases_str.substr(0, bases_str.size() - 1));
            break;
    }

    return bases;
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
            {"pigz", required_argument, 0, 'p'},
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
    std::string pigz = "pigz";
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
                break;
            case 'f':
                fraction = boost::lexical_cast<double>(optarg);
                break;
            case 'p':
                pigz = optarg;
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

    SubsampleSummary summary;

    int64_t total_reads, total_bases;
    FastxCountPair(input1, input2, total_reads, total_bases);

    // mean length of read1 + read2
    double mean_length = static_cast<double>(total_bases) / total_reads * 2.0;

    if (fraction < 0)
    {
        fraction = static_cast<double>(bases) / total_bases;
    }

    summary.expected_subsample_fraction = fraction;

    if (bases < 0)
    {
        bases = fraction * total_bases;
    }

    std::ostringstream pigz_command;
    pigz_command << pigz << " -" << compress_level << " -p" << num_threads;

    FastxSample(input1, input2, output1, output2,
        fraction, bases, mean_length, seed, pigz_command.str(), summary);

    if (bases > 0) {
        summary.expected_subsample_bases = bases;
    } else {
        summary.expected_subsample_bases =
            static_cast<int>(std::round(
                summary.expected_subsample_fraction * summary.total_bases));
    }
    
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
