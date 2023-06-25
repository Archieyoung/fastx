#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <getopt.h>

#include "fastx_subseq.hpp"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/thread_pool.h"
#include "version.hpp"
#include "utils.hpp"


// copy from faidx.c of htslib, expose hidden struct faidx_t
typedef struct {
    int id; // faidx_t->name[id] is for this struct.
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

// copy from faidx.c of htslib, expose hidden struct faidx_t
struct faidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
    enum fai_format_options format;
};


const int64_t FASTX_POS_MAX = std::numeric_limits<int64_t>::max();

#ifdef FASTX_MAX_NAME_LINE
const int64_t FASTX_MAX_NAME_LINE_LEN = FASTX_MAX_NAME_LINE;
#else
const int64_t FASTX_MAX_NAME_LINE_LEN = 4096;
#endif


static
std::vector<std::string> SplitString(std::string str, char delimiter) {
    std::vector<std::string> result;
    size_t i = 0;
    size_t j = 0;
    while ((j = str.find(delimiter, i)) != std::string::npos) {
        result.push_back(str.substr(i, j-i));
        i = j + 1;
    }
    if (i < str.size()) result.push_back(str.substr(i));
    return result;
}


struct Interval {
    Interval(std::string name_, int64_t start_, int64_t end_):
        name(std::move(name_)), start(start_), end(end_) {}
    
    explicit Interval(const std::string &region) {
        if (region.find('\t') != std::string::npos) {
            size_t i = 0;
            i = region.find('\t', i);
            name = region.substr(0, i);
            ++i;
            size_t j = 0;
            j = region.find('\t', i);
            if (j == std::string::npos) {
                std::cerr << "Error! Malformed bed line "
                    << region << std::endl;
                std::exit(1);
            }
            start = SafeStrtol(region.substr(i, j-i).c_str(), 10);
            i = j + 1;
            j = region.find('\t', i);
            if (j == std::string::npos) {
                end = SafeStrtol(region.substr(i).c_str(), 10) - 1;
            } else {
                end = SafeStrtol(region.substr(i, j-i).c_str(), 10) - 1;
            }
        } else if (region.find(':') != std::string::npos) {
            size_t i = 0;
            i = region.find(':', i);
            name = region.substr(0, i);
            ++i;
            size_t j = 0;
            j = region.find('-', i);
            if (j == std::string::npos) {
                std::cerr << "Error! Malformed bed line "
                    << region << std::endl;
                std::exit(1);
            }
            start = SafeStrtol(region.substr(i, j-i).c_str(), 10) - 1;
            i = j + 1;
            end = SafeStrtol(region.substr(i).c_str(), 10) - 1;
        } else {
            name = region;
            start = 0;
            end = FASTX_POS_MAX;
        }
    }

    std::string name;
    // 0-based start position
    int64_t start;
    // 0-based end position
    int64_t end;
};


/**
 * @brief get the name line(include name and comment)
 * 
 * @param fai faidx
 * @param name target sequence name
 * @return std::string name line
 */
std::string FaiGetNameLine(const faidx_t *fai, const char *name) {
    khiter_t iter;

    iter = kh_get(s, fai->hash, name);

    if (iter == kh_end(fai->hash)) {
        std::cerr << "[FaiGetNameLine] Error! The sequence "
            << name << " was not found!" << std::endl;
        std::exit(1);
    }

    faidx1_t val = kh_value(fai->hash, iter);
    int64_t offset = val.seq_offset;
    int ret;

    int64_t offset1 = offset - FASTX_MAX_NAME_LINE_LEN;
    offset1 = offset1 >= 0 ? offset1 : 0;
    ret = bgzf_useek(fai->bgzf, offset1, SEEK_SET);

    if (ret < 0) {
        std::cerr << "[FaiGetNameLine] Error! Failed to get name line! "
            << "bgzf_useek cab not seek to " << offset1
            << ". name=" << name << std::endl;
        std::exit(1);
    }

    std::vector<int> new_line_indices;
    char c;
    int i = 0;
    while ((c=bgzf_getc(fai->bgzf)) >= 0) {
        if (c == '\n') {
            new_line_indices.push_back(i);
        }
        ++i;
        if (i >= offset - offset1) {
            // search to the position before offset
            break;
        }
    }

    if (new_line_indices.empty()) {
        std::cerr << "[FaiGetNameLine] Error! Failed to get name line! "
            << "New line seperator was not found!"
            << " name=" << name << std::endl;
        std::exit(1);
    }

    char name_line[FASTX_MAX_NAME_LINE_LEN];
    if (new_line_indices.size() >= 2 || offset1 == 0) {
        int64_t start;
        if (new_line_indices.size() >= 2) {
            start = new_line_indices[new_line_indices.size()-2] + 1;
        } else {
            start = 0;
        }
        start += offset1;

        ret = bgzf_useek(fai->bgzf, start, SEEK_SET);
        if (ret < 0) {
            std::cerr << "[FaiGetNameLine] Error! Failed to get name line! "
                << "bgzf_useek can not seek to start (" << start
                << ") of name line. name=" << name << std::endl;
            std::exit(1);
        }

        i = 0;
        while ((c=bgzf_getc(fai->bgzf)) >= 0) {
            if (c != '\n' && c != '\r') {
                name_line[i++] = c;
            } else {
                break;
            }
        }
    }

    name_line[i] = '\0';

    return name_line;
}


int FastxSubseq(const faidx_t *fai, const std::vector<Interval> &intervals,
    bool input_name_list, const std::string &output,
    int compress_level, int threads)
{
    BGZF *outfp = NULL;
    hts_tpool *pool = NULL;

    if (output == "-") {
        outfp = bgzf_open(output.c_str(), "wu");

        if (outfp == NULL) {
            std::cerr << "[FastxSubseq] Error! Can not open " << output
                << " for writing" << std::endl;
            perror("");
            std::exit(1);
        }
    }
    else if (output.size() >= 3 && output.substr(output.size()-3) == ".gz") {
        char mode[4] = {0};
        sprintf(mode, "w%d", compress_level);
        outfp = bgzf_open(output.c_str(), mode);

        if (outfp == NULL) {
            std::cerr << "[FastxSubseq] Error! Can not open " << output
                << " for writing" << std::endl;
            perror("");
            std::exit(1);
        }
        
        pool = hts_tpool_init(threads);
        if (pool == NULL) {
            std::cerr << "Error! hts_tpool_init can not init thread pool "
                << std::endl;
            std::exit(1);
        }
        bgzf_thread_pool(outfp, pool, 0);
    } else {
        outfp = bgzf_open(output.c_str(), "wu");

        if (outfp == NULL) {
            std::cerr << "[FastxSubseq] Error! Can not open " << output
                << " for writing" << std::endl;
            perror("");
            std::exit(1);
        }
    }

    int ret;

    for (auto &interval: intervals) {
        int64_t target_len = 0;
        char *seq = faidx_fetch_seq64(fai, interval.name.c_str(),
            interval.start, interval.end, &target_len);
        if (target_len < 0 || seq == nullptr) {
            std::cerr << "[FastxSubseq] Error! Can not fetch "
                "sequence of " << interval.name.c_str();
            if (!input_name_list) {
                std::cerr << ":" << interval.start << "-"
                    << interval.end;
            }
            std::cerr << std::endl;
            free(seq);
            std::exit(1);
        }

        std::string new_name;
        if (!input_name_list) {
            std::ostringstream new_name_stream;
            new_name_stream << ">"
                << interval.name.c_str()
                << ":" << interval.start + 1
                << "-" << interval.end + 1;
            new_name = new_name_stream.str();
        } else {
            new_name = FaiGetNameLine(fai, interval.name.c_str());
        }

        ret = bgzf_write(outfp, new_name.c_str(), new_name.size());
        if (ret < 0 || ret < (int)new_name.size()) {
            std::cerr << "[FastxSubseq] Error! failed to write fasta name line "
                << new_name << std::endl;
            std::exit(1);
        }

        ret = bgzf_write(outfp, "\n", 1);
        if (ret < 0) {
            std::cerr << "[FastxSubseq] Error! failed to write new line"
                << std::endl;
            std::exit(1);
        }

        ret = bgzf_write(outfp, seq, target_len);
        if (ret < 0) {
            std::cerr << "[FastxSubseq] Error! failed to write fasta seq for "
                << "name=" << interval.name << std::endl;
            std::exit(1);
        }

        ret = bgzf_write(outfp, "\n", 1);
        if (ret < 0) {
            std::cerr << "[FastxSubseq] Error! failed to write new line"
                << std::endl;
            std::exit(1);
        }

        free(seq);        
    }

    if (outfp) bgzf_close(outfp);
    if (pool) hts_tpool_destroy(pool);

    return 0;
}


static
void Usage() {
    std::cerr << "fastx subseq " << FASTX_VERSION << std::endl;
    std::cerr << std::endl;
    std::cerr << "  extract subsequences of fasta/fastq.\n"
              << std::endl;
    std::cerr
            << "Usage: fastx subseq [options] <file.fasta>\n\n"
            << "Options:\n"
            << "  -o, --output, FILE          output file name [stdout]\n"
            << "  -r, --region, STR           comma-separated list of regions\n"
            << "  -R, --region-file, FILE     regions list in file(can be bed file or target name list)\n"
            << "  -l, --level, INT            compression level(0 to 9, or 11), valid if output file type is gzip [6]\n"
            << "  -t, --thread, INT           number of threads for compression, valid if output file type is gzip [4]\n"
            << "  -h, --help                  print this message and exit.\n"
            << "  -V, --version               print version."
            << std::endl;
}


int FastxSubseqMain(int argc, char **argv)
{
    if (argc == 1)
    {
        Usage();
        return 0;
    }

    static const struct option long_options[] = {
            {"output", required_argument, 0, 'o'},
            {"region", required_argument, 0, 'r'},
            {"region-file", required_argument, 0, 'R'},
            {"level", required_argument, 0, 'l'},
            {"thread", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "o:r:R:l:t:hV";

    std::string input;
    std::string output = "-";
    std::string region;
    std::string region_file;
    int compress_level = 6;
    int num_threads = 4;

    while ((c = getopt_long(
        argc, argv, opt_str, long_options, &long_idx)) != -1)
    {
        switch (c) {
            case 'o':
                output = optarg;
                break;
            case 'r':
                region = optarg;
                break;
            case 'R':
                region_file = optarg;
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

    if (optind < argc) {
        input = argv[optind];
    } else {
        std::cerr << "Error! Missing input file" << std::endl;
        std::exit(1);
    }

    if (region.empty() && region_file.empty()) {
        std::cerr << "Error! Must set regions by -r(--region) or "
            << "-R(--region-file)" << std::endl;
        std::exit(1);
    }

    if (!region.empty() && !region_file.empty()) {
        std::cerr << "Error! -r(--region) is conflict with "
            << "-R(--region-file)" << std::endl;
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

    std::vector<Interval> intervals;
    if (!region.empty()) {
        std::vector<std::string> regions = SplitString(region, ',');
        for (auto &r: regions) {
            intervals.emplace_back(r);
        }
    } else {
        std::ifstream instream(region_file);
        std::string line;
        while (std::getline(instream, line)) {
            intervals.emplace_back(line);
        }
    }

    bool is_name_list = true;
    for (auto &i: intervals) {
        if (i.end != FASTX_POS_MAX) {
            is_name_list = false;
        }
    }

    faidx_t *fai = fai_load(input.c_str());
    if (fai == nullptr) {
        std::cerr << "Error! Fail to load fai index for fasta " << input
            << std::endl;
        std::exit(1);
    }

    FastxSubseq(fai, intervals, is_name_list, output,
        compress_level, num_threads);

    fai_destroy(fai);
    return 0;
}
