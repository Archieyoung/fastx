#ifndef FASTX_SEQ_READER_HPP
#define FASTX_SEQ_READER_HPP


#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <zlib.h>

#include "kseq_utils.hpp"


const int KSEQ_ARRAY_CAPACITY = 1024;
const int SEQ_READER_QUEUE_SIZE = 10;


/**
 * @brief fixed size array of kseq_t *
 * 
 */
class KseqArray {
public:
    KseqArray(kstream_t *kstream, size_t capacity = KSEQ_ARRAY_CAPACITY):
        cur_(0), size_(0), capacity_(capacity)
    {
        data_ = (kseq_t **)malloc(capacity_*sizeof(kseq_t *));
        for (int i = 0; i < capacity_; ++i) {
            data_[i] = (kseq_t*)calloc(1, sizeof(kseq_t));	
            data_[i]->f = kstream;
        }
    }

    ~KseqArray() {
        if (data_) {
            for (int i = 0; i < capacity_; ++i) {
                if (data_[i]) {
                    free(data_[i]->name.s);
                    free(data_[i]->comment.s);
                    free(data_[i]->seq.s);
                    free(data_[i]->qual.s);
                    free(data_[i]);
                }
            }
            free(data_);
        }
    }


    kseq_t *get(int pos) {
        if (pos > size_ - 1 || pos < 0) {
            std::cerr << "[KseqArray::get] Out of range error! pos: " << pos
                << "size: " << size_;
            std::exit(1);
        }
        return data_[pos];
    }

    kseq_t *get() {
        return data_[cur_++];
    }

    int set(int pos, kseq_t *ks_state) {
        if (pos > capacity_ - 1 || pos < 0) {
            std::cerr << "[KseqArray::set] Out of range error! pos: " << pos
                << "capacity: " << capacity_;
            std::exit(1);
        }
        
        data_[pos]->last_char = ks_state->last_char;
        data_[pos]->is_fastq = ks_state->is_fastq;
        int64_t r = kseq_read(data_[pos]);
        ks_state->last_char = data_[pos]->last_char;
        ks_state->is_fastq = data_[pos]->is_fastq;
        
        // reach end of file
        if (r == -1) return r;
        
        if (r < -1) {
            std::cerr << "[KseqArray::set] Failed to read seq from stream! "
                << "kseq_read return code " << r << std::endl;
            std::exit(1);
        }

        size_ = pos + 1;

        return 0;
    }

    int &cur() {
        return cur_;
    }

    int size() {
        return size_;
    }

    int capacity() {
        return capacity_;
    }

    void clear() {
        cur_ = 0;
        size_ = 0;
    }

    bool empty() {
        return cur_ == size_;
    }

private:
    kseq_t **data_ = nullptr;
    int cur_;
    int size_;
    int capacity_;
};



class SeqReader {
public:
    SeqReader(gzFile fp): stop_(false)
    {
        ks_ = kseq_init(fp);
        for (int i = 0; i < SEQ_READER_QUEUE_SIZE; ++i) {
            KseqArray *kseq_array = new KseqArray(ks_->f);
            empty_queue_.push(kseq_array);
        }

        // start a reading thread
        producer_ = std::thread([this]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mutex_);
                producer_cv_.wait(lock,
                    [this]{return !empty_queue_.empty() || stop_;});
                if (stop_) break;
                KseqArray *kseq_array = empty_queue_.front();
                empty_queue_.pop();
                int r = Fill(kseq_array);
                filled_queue_.push(kseq_array);
                if (r < 0) {
                    // reach end of file
                    stop_ = true;
                    consumer_cv_.notify_one();
                    break;
                }
                consumer_cv_.notify_one();
            }
        });
    }

    ~SeqReader() {
        stop_ = true;
        producer_cv_.notify_one();
        if (producer_.joinable()) producer_.join();

        if (filled_queue_.size() + empty_queue_.size() < SEQ_READER_QUEUE_SIZE)
        {
            delete reading_array_;
        }

        while (!filled_queue_.empty()) {
            delete filled_queue_.front();
            filled_queue_.pop();
        }

        while (!empty_queue_.empty()) {
            delete empty_queue_.front();
            empty_queue_.pop();
        }

        kseq_destroy(ks_);
    }

    kseq_t *read() {
        if (reading_array_ == nullptr) {
            std::unique_lock<std::mutex> lock(mutex_);
            consumer_cv_.wait(lock,
                [this](){return !filled_queue_.empty() || stop_;});
            if (!filled_queue_.empty()) {
                reading_array_ = filled_queue_.front();
                filled_queue_.pop();
                producer_cv_.notify_one();
            } else {
                std::cerr << "[SeqReader::read] Error! Failed to init reading"
                    << " array!" << std::endl;
                std::exit(1);
            }
        }

        if (reading_array_->empty()) {
            std::unique_lock<std::mutex> lock(mutex_);
            consumer_cv_.wait(lock,
                [this](){return !filled_queue_.empty() || stop_;});
            
            reading_array_->clear();
            empty_queue_.push(reading_array_);

            if (!filled_queue_.empty()) {
                reading_array_ = filled_queue_.front();
                filled_queue_.pop();
                producer_cv_.notify_one();
            }
        }

        if (reading_array_->empty()) {
            return nullptr;            
        } else {
            return reading_array_->get();
        }
    }

    void stop() {
        stop_ = true;
        producer_cv_.notify_one();
    }

private:

    int Fill(KseqArray * kseq_array) {
        for (int i = 0; i < kseq_array->capacity(); ++i) {
            int r = kseq_array->set(i, ks_);
            if (r < 0) return r;
        }
        return 0;
    }

    kseq_t *ks_;
    std::queue<KseqArray *> filled_queue_;
    std::queue<KseqArray *> empty_queue_;
    KseqArray *reading_array_ = nullptr;
    std::thread producer_;

    std::atomic_bool stop_;
    std::mutex mutex_;
    std::condition_variable producer_cv_;
    std::condition_variable consumer_cv_;
};


#endif  // FASTX_SEQ_READER_HPP
