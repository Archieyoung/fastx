#include "common.hpp"
#include "htslib/bgzf.h"
#include "seq_reader.hpp"
#include "zlib.h"
#include <cstdio>
#include "htslib/thread_pool.h"

int main(int argc, char **argv) {
    char *input = argv[1];
    gzFile fp = gzopen(input, "r");
    if (fp == NULL) {
        std::cerr << "Error! Failed to open " << argv[1] << " for reading."
            << std::endl;
        perror("");
        std::exit(1);
    }

    SeqReader reader = SeqReader(fp);

    BGZF* bgzfp = bgzf_open(argv[2], "w6");
    hts_tpool *pool = hts_tpool_init(atoi(argv[3]));
    bgzf_thread_pool(bgzfp, pool, 0);

    // int N = 0;
    kseq_t *ks = nullptr;
    while ((ks = reader.read()) != nullptr) {
        BgzfWriteKseq(bgzfp, ks);
        // ++N;
        // if (N >= 1000) break;
    }

    bgzf_close(bgzfp);
    hts_tpool_destroy(pool);
    gzclose(fp);
}
