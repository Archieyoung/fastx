#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "fastx_subseq.hpp"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"


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

struct Interval {
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
        std::cerr << "Error! The sequence " << name << " was not found!"
            << std::endl;
        std::exit(1);
    }

    faidx1_t val = kh_value(fai->hash, iter);
    int64_t offset = val.seq_offset;
    int ret;

    int64_t offset1 = offset - FASTX_MAX_NAME_LINE_LEN;
    offset1 = offset1 >= 0 ? offset1 : 0;
    ret = bgzf_useek(fai->bgzf, offset1, SEEK_SET);

    if (ret < 0) {
        std::cerr << "Error! Failed to get name line!" << std::endl;
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
        std::cerr << "Error! Failed to get name line! "
            << "New line seperator was not found!" << std::endl;
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
            std::cerr << "Error! Failed to get name line!" << std::endl;
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
    bool input_name_list)
{
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
            new_name_stream << "<"
                << interval.name.c_str()
                << ":" << interval.start + 1
                << "-" << interval.end + 1
                << std::endl;
            new_name = new_name_stream.str();
        } else {
            new_name = FaiGetNameLine(fai, interval.name.c_str());
        }

        std::cout << new_name << "\n" << seq << std::endl;
        free(seq);        
    }
    return 0;
}


int main(int argc, char **argv) {
    faidx_t *fai = fai_load(argv[1]);
    if (fai == nullptr) {
        std::cerr << "Error! Fail to load fai index for fasta " << argv[1]
            << std::endl;
        std::exit(1);
    }
    
    // for (int i = 0; i < fai->n; ++i)
    //     std::cout << FaiGetNameLine(fai, fai->name[i]) << std::endl;

    std::vector<Interval> intervals;
    std::ifstream instream(argv[2]);
    std::string line;
    while (std::getline(instream, line)) {
        intervals.push_back(Interval{line, 0, FASTX_POS_MAX});
    }

    FastxSubseq(fai, intervals, true);

    fai_destroy(fai);
    return 0;
}

