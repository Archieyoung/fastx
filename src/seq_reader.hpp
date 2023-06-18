#ifndef FASTX_SEQ_READER_HPP
#define FASTX_SEQ_READER_HPP


#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <zlib.h>

#include "common.hpp"


/**
 * @brief a single producer single consumer sequence reader implemented using
 * ring buffer. It uses one extra thread for reading.
 * 
 */
class SeqReader {
public:
    SeqReader(gzFile fp, int size = 1024): size_(size), head_(0), tail_(0),
        is_eof_(false)
    {
        ks_ = kseq_init(fp);
        buffer_ = (kseq_t **)malloc(size*sizeof(kseq_t *));
        for (int i = 0; i < size_; ++i) {
            buffer_[i] = KseqInitWithKstream(ks_->f);
        }

        producer_ = std::thread(
            [this]() {
                while (!this->is_eof_) {
                    this->Push();
                }
            }
        );
    }

    ~SeqReader() {
        if (buffer_) {
            for (int i = 0; i < size_; ++i) {
                if (buffer_[i]) {
                    free(buffer_[i]->name.s);
                    free(buffer_[i]->comment.s);
                    free(buffer_[i]->seq.s);
                    free(buffer_[i]->qual.s);
                    free(buffer_[i]);
                }
            }
            free(buffer_);
        }

        kseq_destroy(ks_);
        if (producer_.joinable()) producer_.join();
    }

    kseq_t *read() {
        return Pop();
    }

    void join() {
        if (producer_.joinable()) producer_.join();
    }

private:

    bool IsEmpty() const {
        return head_ == tail_;
    }

    bool IsFull() const {
         return (tail_ + 1) % size_ == head_;
    }

    void Push() {
        while (IsFull()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        buffer_[tail_]->last_char = ks_->last_char;
        buffer_[tail_]->is_fastq = ks_->is_fastq;
        int64_t r = kseq_read(buffer_[tail_]);
        ks_->last_char = buffer_[tail_]->last_char;
        ks_->is_fastq = buffer_[tail_]->is_fastq;
        tail_ = (tail_ + 1) % size_;
        if (r == -1) {
            is_eof_ = true;
            return;
        }
    }

    kseq_t *Pop() {
        if (is_eof_ && IsEmpty()) {
            return nullptr;
        }
        while (IsEmpty()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        kseq_t *ks = buffer_[head_];
        head_ = (head_ + 1) % size_;
        return ks;
    }

    kseq_t *KseqInitWithKstream(kstream_t *f) {
        kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));	
        s->f = f;
        return s;
    }

    kseq_t *ks_;
    int size_;
    kseq_t **buffer_;
    std::atomic_int head_;
    std::atomic_int tail_;

    std::atomic_bool is_eof_;

    std::thread producer_;
};


#endif  // FASTX_SEQ_READER_HPP
