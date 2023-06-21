#ifndef FASTX_KSEQ_UTILS_HPP
#define FASTX_KSEQ_UTILS_HPP


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


int BgzfWriteKseq(BGZF *fp, const kseq_t *seq);

/**
 * @brief detect if the input file is a fasta/fastq
 * 
 * @param path input file path
 * @return true if the input file is a fastq
 * @return false if the input file is a fastq
 */
bool IsFastq(const char *path);

#endif  // FASTX_KSEQ_UTILS_HPP
