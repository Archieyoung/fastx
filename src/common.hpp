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

/**
 * @brief safe wrapper of strtol.
 * 
 * @param str C-string beginning with the representation of an integral number.
 * @param base Numerical base (radix) that determines the valid characters and their interpretation.
 * @return long 
 */
long SafeStrtol(const char *str, int base);


/**
 * @brief safe wrapper of strtod.
 * 
 * @param str C-string beginning with the representation of a floating-point number.
 * @return double 
 */
double SafeStrtod(const char *str);


/**
 * @brief detect if the input file is a fasta/fastq
 * 
 * @param path input file path
 * @return true if the input file is a fastq
 * @return false if the input file is a fastq
 */
bool IsFastq(const char *path);



#endif  // FASTX_COMMON_HPP
