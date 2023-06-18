#ifndef FASTX_SEQ_HPP
#define FASTX_SEQ_HPP


#include "htslib/kstring.h"


/**
 * @brief fasta/fastq sequence
 * 
 */
struct FastxSeq {

    kstring_t *name;
    kstring_t *comment;
    kstring_t *seq;
    kstring_t *qual;
};


#endif  // FASTX_SEQ_HPP
