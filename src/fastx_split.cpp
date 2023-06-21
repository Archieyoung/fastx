#include <cstdint>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <zlib.h>
#include <getopt.h>
#include <thread>
#include "utils.hpp"
#include "kseq_utils.hpp"
#include "htslib/bgzf.h"
#include "seq_reader.hpp"
#include "version.hpp"
#include "htslib/thread_pool.h"


void FastxSplitReads(const std::string &ifilename, int64_t n_read_per_chunk,
    int64_t n_base_per_chunk, const std::string &prefix,
    const std::string &suffix, int threads, int compress_level)
{
    bool split_by_base = false;
    if (n_read_per_chunk <= 0) {
        if (n_base_per_chunk > 0) {
            split_by_base = true;
        } else {
            std::cerr << "Error! Must input a least one of "
                << "bases(-b, --bases) or reads(-r, --reads)"
                << std::endl;
            std::exit(1);
        }
    }

    gzFile fp = gzopen(ifilename.c_str(), "r");

    kseq_t *read = kseq_init(fp);

    int ret1, ret2;
    int64_t read_count = 0;
    int64_t base_count = 0;
    bool open_new = false;
    int64_t n = 0;

    std::ostringstream ofilename;
    ofilename << prefix << "." << n << "." << suffix;
    hts_tpool *pool = hts_tpool_init(threads);
    if (pool == NULL) {
        std::cerr << "Error! hts_tpool_init can not init thread pool "
            << std::endl;
        std::exit(1);
    }

    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp = bgzf_open(ofilename.str().c_str(), mode_str.str().c_str());
    if (bgzfp == NULL) {
        std::cerr << "Error! Can not open " << ofilename.str()
            << " for writing" << std::endl;
        std::exit(1);
    }
    bgzf_thread_pool(bgzfp, pool, 0);

    while ((ret1 = kseq_read(read)) >= 0)
    {
        if (!open_new) {
            ret2 = BgzfWriteKseq(bgzfp, read);
            if (ret2 < 0) {
                std::cerr << "Error! Failed to write read: "
                    << read->name.s << " to " << ofilename.str() << std::endl;
                std::exit(1);
            }
            if (split_by_base) {
                base_count += read->seq.l;
                if (base_count >= n_base_per_chunk) open_new = true;
            } else {
                ++read_count;
                if (read_count >= n_read_per_chunk) open_new = true;
            }
        } else {
            ++n;
            ofilename.str("");   // clear
            bgzf_close(bgzfp);
            ofilename << prefix << "." << n << "." << suffix;
            bgzfp = bgzf_open(ofilename.str().c_str(), mode_str.str().c_str());
            if (bgzfp == NULL) {
                std::cerr << "Error! Can not open "
                    << ofilename.str() << " for writing" << std::endl;
                std::exit(1);
            }

            bgzf_thread_pool(bgzfp, pool, 0);

            ret2 = BgzfWriteKseq(bgzfp, read);
            if (ret2 < 0) {
                std::cerr << "Error! Failed to write read: "
                    << read->name.s << " to " << ofilename.str() << std::endl;
                std::exit(1);
            }

            if (split_by_base) {
                base_count = 0;
                base_count += read->seq.l;
                if (base_count < n_base_per_chunk) open_new = false;
            } else {
                read_count = 0;
                ++read_count;
                if (read_count < n_read_per_chunk) open_new = false;
            }
        }
    }

    if (ret1 == -2) {
        std::cerr << "Error! Input fastq file quality is truncated "
            << ifilename << std::endl;
        std::exit(1);
    } else if (ret1 == -3) {
        std::cerr << "Error! Input fastq file stream error!"
            << ifilename << std::endl;
        std::exit(1);
    }

    bgzf_close(bgzfp);
    hts_tpool_destroy(pool);
    kseq_destroy(read);
    gzclose(fp);
}

/**
 * @brief split paired end fastq/fasta by number of bases
 * 
 * @param input1 first input fastq/fasta file path
 * @param input2 second input fastq/fasta file path
 * @param n_base_per_chunk number of bases per chunk
 * @param prefix output prefix
 * @param suffix output suffix
 * @param threads number of threads
 * @param compress_level compress level
 */
void FastxSplitReadsByBasesPair(
    const std::string &input1, const std::string &input2,
    int64_t n_base_per_chunk, const std::string &prefix,
    const std::string &suffix, int threads, int compress_level)
{
    gzFile fp1 = gzopen(input1.c_str(), "r");
    gzFile fp2 = gzopen(input2.c_str(), "r");

    SeqReader reader1 = SeqReader(fp1);
    SeqReader reader2 = SeqReader(fp2);

    int64_t base_count = 0;
    bool open_new = false;
    int64_t n = 0;

    std::ostringstream ofilename1;
    ofilename1 << prefix << "." << n << ".R1." << suffix;
    std::ostringstream ofilename2;
    ofilename2 << prefix << "." << n << ".R2." << suffix;

    hts_tpool *pool = hts_tpool_init(threads);
    if (pool == NULL) {
        std::cerr << "Error! hts_tpool_init can not init thread pool "
            << std::endl;
        std::exit(1);
    }

    std::ostringstream mode_str;
    mode_str << "w" << compress_level;
    BGZF* bgzfp1 = bgzf_open(ofilename1.str().c_str(), mode_str.str().c_str());
    if (bgzfp1 == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename1.str() << " for writing" << std::endl;
        std::exit(1);
    }

    BGZF* bgzfp2 = bgzf_open(ofilename2.str().c_str(), mode_str.str().c_str());
    if (bgzfp2 == NULL) {
        std::cerr << "Error! Can not open "
            << ofilename2.str() << " for writing" << std::endl;
        std::exit(1);
    }


    bgzf_thread_pool(bgzfp1, pool, 0);
    bgzf_thread_pool(bgzfp2, pool, 0);

    int ret;
    kseq_t *read1 = nullptr;
    kseq_t *read2 = nullptr;
    while ((read1 = reader1.read()) != nullptr &&
        (read2 = reader2.read()) != nullptr)
    {
        if (!open_new) {
            ret = BgzfWriteKseq(bgzfp1, read1);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read1: "
                    << read1->name.s << " to " << ofilename1.str() << std::endl;
                std::exit(1);
            }
            ret = BgzfWriteKseq(bgzfp2, read2);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read2: "
                    << read2->name.s << " to " << ofilename2.str() << std::endl;
                std::exit(1);
            }
            base_count += read1->seq.l;
            base_count += read2->seq.l;
            if (base_count >= n_base_per_chunk) {
                open_new = true;
            }
        } else {
            ++n;
            ofilename1.str("");   // clear
            ofilename2.str("");   // clear
            bgzf_close(bgzfp1);
            bgzf_close(bgzfp2);
            ofilename1 << prefix << "." << n << ".R1." << suffix;
            ofilename2 << prefix << "." << n << ".R2." << suffix;
            
            bgzfp1 = bgzf_open(ofilename1.str().c_str(),
                mode_str.str().c_str());
            if (bgzfp1 == NULL) {
                std::cerr << "Error! Can not open "
                    << ofilename1.str() << " for writing" << std::endl;
                std::exit(1);
            }

            bgzfp2 = bgzf_open(ofilename2.str().c_str(),
                mode_str.str().c_str());
            if (bgzfp2 == NULL) {
                std::cerr << "Error! Can not open "
                    << ofilename2.str() << " for writing" << std::endl;
                std::exit(1);
            }
            
            bgzf_thread_pool(bgzfp1, pool, 0);
            bgzf_thread_pool(bgzfp2, pool, 0);

            ret = BgzfWriteKseq(bgzfp1, read1);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read1: "
                    << read1->name.s << " to " << ofilename1.str() << std::endl;
                std::exit(1);
            }
            ret = BgzfWriteKseq(bgzfp2, read2);
            if (ret < 0) {
                std::cerr << "Error! Failed to write read2: "
                    << read2->name.s << " to " << ofilename2.str() << std::endl;
                std::exit(1);
            }
            
            base_count = 0;
            base_count += read1->seq.l;
            base_count += read2->seq.l;
            if (base_count < n_base_per_chunk) {
                open_new = false;
            }
        }
    }


    bgzf_close(bgzfp1);
    bgzf_close(bgzfp2);
    hts_tpool_destroy(pool);
    gzclose(fp1);
    gzclose(fp2);
}


static
void Usage() {
    std::cerr << "fastx split " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  split fasta/fastq files.\n"
              << std::endl;
    std::cerr
            << "Options:\n"
            << "  -i, --in1, FILE             input fasta/fastq file name for read1.\n"
            << "  -I, --in2, FILE             input fasta/fastq file name for read2(optional).\n"
            << "  -p, --prefix, FILE          output fasta/fastq file name prefix.\n"
            << "  -b, --bases, STR            put this value of bases per output file(K/M/G).\n"
            << "  -n, --reads, STR            put this value of reads per output file(K/M/G).\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11).[6]\n"
            << "  -t, --thread, INT           number of threads.[4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}


int FastxSplitMain(int argc, char **argv)
{
    if (argc == 1)
    {
        Usage();
        return 0;
    }

    static const struct option long_options[] = {
            {"in1", required_argument, 0, 'i'},
            {"in2", required_argument, 0, 'I'},
            {"prefix", required_argument, 0, 'p'},
            {"bases", required_argument, 0, 'b'},
            {"reads", required_argument, 0, 'n'},
            {"level", required_argument, 0, 'l'},
            {"thread", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "i:I:p:b:n:l:t:hV";

    std::string input1 = "";
    std::string input2 = "";
    std::string prefix = "";
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
            case 'p':
                prefix = optarg;
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

    if (prefix.empty()) {
        std::cerr << "Error! Must set the output prefix using -p(--prefix). "
            << std::endl;
        std::exit(1);
    }

    if (bases <= 0 && reads <= 0)
    {
        std::cerr << "Error! Must input a least one of "
            << "bases(-b, --bases) or reads(-r, --reads)"
            << std::endl;
        std::exit(1);
    }

    if (bases > 0 && reads > 0) {
        std::cerr << "Error! bases(-b, --bases) is conflict with "
            << "reads(-r, --reads)." << std::endl;
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
        // single end reads
        bool input1_is_fq = IsFastq(input1.c_str());
        std::string input1_suffix = "fastq.gz";
        if (!input1_is_fq) {
            input1_suffix = "fasta.gz";
        }
        FastxSplitReads(input1, reads, bases, prefix, input1_suffix,
            num_threads, compress_level);
    } else {
        // paired end reads
        bool input1_is_fq = IsFastq(input1.c_str());
        std::string input1_suffix = "fastq.gz";
        if (!input1_is_fq) {
            input1_suffix = "fasta.gz";
        }
        bool input2_is_fq = IsFastq(input2.c_str());
        std::string input2_suffix = "fastq.gz";
        if (!input2_is_fq) {
            input2_suffix = "fasta.gz";
        }

        if (input1_is_fq != input2_is_fq) {
            std::cerr << "Error! fastx do not support mixed input of fasta and "
                << "fastq for input1 and input2." << std::endl;
            std::exit(1);
        }

        if (reads > 0) {
            FastxSplitReads(input1, reads, bases, prefix, "R1."+input1_suffix,
                num_threads, compress_level);
            FastxSplitReads(input2, reads, bases, prefix, "R2."+input2_suffix,
                num_threads, compress_level);
        } else {
            FastxSplitReadsByBasesPair(input1, input2, bases, prefix,
                input1_suffix, num_threads, compress_level);
        }
    }

    return 0;
}
