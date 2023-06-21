#include <climits>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <iostream>
#include <sstream>
#include <cmath>
#include "kseq_utils.hpp"

void FastxCount(const std::string &filename, int64_t &reads, int64_t &bases)
{
    reads = 0;
    bases = 0;
    gzFile fp = gzopen(filename.c_str(), "r");

    if (fp == nullptr)
    {
        std::perror(("Error! Can not open " + filename).c_str());
        std::exit(1);
    }

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

    if (read1_counts != read2_counts) {
        std::cerr << "Error! Record number not equal for paired inputs."
            << "Input1 has " << read1_counts << " records, but "
            << "Input1 has " << read2_counts << " records" << std::endl;
        std::exit(1);
    }

    reads = read1_counts + read2_counts;
    bases = read1_bases + read2_bases; 
}


std::string kseqToStr(const kseq_t *seq) {
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


bool IsFastq(const char *path)
{
    gzFile fp = gzopen(path, "r");
    kseq_t *read = kseq_init(fp);

    int ret = kseq_read(read);
    if (ret >= 0) {
        if (read->qual.l) {
            kseq_destroy(read);
            gzclose(fp);
            return true;
        }
    } else if (ret == -2) {
        std::cerr << "Error! Input fastq file quality is truncated "
            << path << std::endl;
        std::exit(1);
    } else if (ret == -3) {
        std::cerr << "Error! Input fastq file stream error!"
            << path << std::endl;
        std::exit(1);
    }

    kseq_destroy(read);
    gzclose(fp);
    return false;
}
