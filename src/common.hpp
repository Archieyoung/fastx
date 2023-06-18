#ifndef FASTX_COMMON_HPP
#define FASTX_COMMON_HPP

#include <cstdint>
#include <string>
#include "zlib.h"
#include "kseq.h"
#include "htslib/bgzf.h"

KSEQ_INIT(gzFile, gzread)

// count fasta/q reads and bases
void FastxCount(const std::string &filename, int64_t &reads, int64_t &bases);

// count fasta/q reads and bases using two thread
void FastxCountPair(const std::string &ifilename1,
    const std::string &ifilename2, int64_t &reads, int64_t &bases);


std::string kseqToStr(const kseq_t *seq);

int64_t BasesStrToInt(const std::string &bases_str);


int BgzfWriteKseq(BGZF *fp, const kseq_t *seq);

#endif  // FASTX_COMMON_HPP
