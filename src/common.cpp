#include <thread>
#include <iostream>
#include <sstream>
#include "common.hpp"
#include "boost/lexical_cast.hpp"


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
            bases = static_cast<int64_t>(std::round(1000000000 *
                boost::lexical_cast<double>(
                bases_str.substr(0, bases_str.size() - 1))));
            break;
        case 'M':
        case 'm':
            bases = static_cast<int64_t>(std::round(1000000 *
                boost::lexical_cast<double>(
                bases_str.substr(0, bases_str.size() - 1))));
            break;
        case 'K':
        case 'k':
            bases = static_cast<int64_t>(std::round(1000 *
                boost::lexical_cast<double>(
                bases_str.substr(0, bases_str.size() - 1))));
            break;
        default:
            bases = boost::lexical_cast<int64_t>(
                bases_str.substr(0, bases_str.size()));
            break;
    }

    return bases;
}
