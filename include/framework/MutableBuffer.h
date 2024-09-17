/*
 * include/framework/MutableBuffer.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Revised BSD License.
 *
 */
#pragma once

#include <cstdlib>
#include <atomic>
#include <cassert>
#include <numeric>
#include <algorithm>

#include "base/weighted_sampling.h"
#include "util/alignment.h"
#include "base/bf_config.h"
#include "ds/BloomFilter.h"
#include "base/record.h"
#include "ds/Alias.h"
#include "util/timer.h"

namespace extension {


class MutableBuffer {
public:
    MutableBuffer(size_t capacity, size_t max_tombstone_cap)
    : m_cap(capacity), m_tombstone_cap(max_tombstone_cap), m_reccnt(0)
    , m_tombstonecnt(0), m_weight(0), m_max_weight(0) {
        auto len = capacity * sizeof(record_t);
        size_t aligned_buffersize = len + (psudb::CACHELINE_SIZE - (len % psudb::CACHELINE_SIZE));
        m_data = (record_t*) std::aligned_alloc(psudb::CACHELINE_SIZE, aligned_buffersize);
        m_tombstone_filter = nullptr;
        if (max_tombstone_cap > 0) {
            m_tombstone_filter = new psudb::BloomFilter<skey_t>(BF_FPR, max_tombstone_cap, BF_HASH_FUNCS);
        }
    }

    ~MutableBuffer() {
        if (m_data) free(m_data);
        if (m_tombstone_filter) delete m_tombstone_filter;
    }

    int append(const key_t& key, const value_t& value, weight_t weight = 1, bool is_tombstone = false) {
        if (is_tombstone && m_tombstonecnt + 1 > m_tombstone_cap) return 0;

        int32_t pos = 0;
        if ((pos = try_advance_tail()) == -1) return 0;

        if (is_tombstone) {
            weight = 0;
        }

        m_data[pos].key = key;
        m_data[pos].value = value;
        m_data[pos].header = ((pos << 2) | (is_tombstone ? 1 : 0));
        m_data[pos].weight = weight;

        if (is_tombstone) {
            m_tombstonecnt.fetch_add(1);
            if (m_tombstone_filter) m_tombstone_filter->insert(key);
        }

        double old_val, new_val;
        do {
            old_val = m_weight.load();
            new_val = old_val + weight;
        } while (!m_weight.compare_exchange_strong(old_val, new_val));


        double old = m_max_weight.load();
        while (old < weight) {
            m_max_weight.compare_exchange_strong(old, weight);
            old = m_max_weight.load();
        }

        return 1;     
    }

    bool truncate() {
        m_tombstonecnt.store(0);
        m_reccnt.store(0);
        m_weight.store(0);
        m_max_weight.store(0);
        if (m_tombstone_filter) m_tombstone_filter->clear();

        return true;
    }

    record_t* sorted_output() {
        TIMER_INIT();
        TIMER_START();
        std::sort(m_data, m_data + m_reccnt.load(), buffer_record_cmp);
        TIMER_STOP();

        #ifdef INSTRUMENT_MERGING
        fprintf(stderr, "sort\t%ld\n", TIMER_RESULT());
        #endif
        return m_data;
    }
    
    size_t get_record_count() {
        return m_reccnt;
    }
    
    size_t get_capacity() {
        return m_cap;
    }

    bool is_full() {
        return m_reccnt == m_cap;
    }

    size_t get_tombstone_count() {
        return m_tombstonecnt.load();
    }

    bool delete_record(const record_t &rec) {
        auto offset = 0;
        while (offset < m_reccnt.load()) {
            if (m_data[offset] == rec) {
                m_data[offset].set_delete_status();
                return true;
            }
            offset++;
        }

        return false;
    }

    bool check_tombstone(const record_t &rec) {
        if (m_tombstone_filter && !m_tombstone_filter->lookup(rec.key)) return false;

        auto offset = 0;
        while (offset < m_reccnt.load()) {
            if (m_data[offset] == rec && m_data[offset].is_tombstone()) return true;
            offset++;;
        }
        return false;
    }

    const record_t* get_record_at(size_t idx) {
        return m_data + idx;
    }

    size_t get_memory_usage() {
        return m_cap * sizeof(record_t);
    }

    size_t get_aux_memory_usage() {
        return m_tombstone_filter->memory_usage();
    }

    double get_sample_range(const key_t &lower, const key_t &upper, psudb::Alias **alias, 
                            std::vector<record_t *> &records, size_t *cutoff) {
      *cutoff = std::atomic_load(&m_reccnt) - 1;
      std::vector<double> weights;

      records.clear();

      double tot_weight = 0.0;
      for (size_t i = 0; i < (*cutoff) + 1; i++) {
        if (!m_data[i].is_tombstone() && !m_data[i].get_delete_status()) {
            if (m_data[i].key >= lower && m_data[i].key <= upper) {
                tot_weight += m_data[i].weight;
                weights.push_back(m_data[i].weight);
                records.push_back(&m_data[i]);
            }
        } 
      }

      for (size_t i = 0; i < weights.size(); i++) {
        weights[i] = weights[i] / tot_weight;
      }

      *alias = new psudb::Alias(weights);

      return tot_weight;
    }

    // rejection sampling
    record_t* get_sample(const key_t &lower, const key_t &upper, gsl_rng *rng) {
        size_t reccnt = m_reccnt.load();
        if (reccnt == 0) {
            return nullptr;
        }

        auto idx = (reccnt == 1) ? 0 : gsl_rng_uniform_int(rng, reccnt - 1);
        auto test = gsl_rng_uniform(rng) * m_max_weight.load();

        // reject tombstones and deleted records automatically
        if (m_data[idx].is_tombstone() || m_data[idx].get_delete_status()) {
            return nullptr;
        }

        // next reject based on key range
        if (m_data[idx].key < lower || m_data[idx].key > upper) {
            return nullptr;
        }

        // finally, do the weighted rejection check
        return (test <= m_data[idx].weight) ? m_data + idx : nullptr;
    }

    size_t get_tombstone_capacity() {
        return m_tombstone_cap;
    }

    double get_total_weight() {
        return m_weight.load();
    }

private:
    int32_t try_advance_tail() {
        size_t new_tail = m_reccnt.fetch_add(1);

        if (new_tail < m_cap) return new_tail;
        else return -1;
    }

    size_t m_cap;
    size_t m_tombstone_cap;
    
    record_t* m_data;
    psudb::BloomFilter<skey_t>* m_tombstone_filter;

    alignas(64) std::atomic<size_t> m_tombstonecnt;
    alignas(64) std::atomic<uint32_t> m_reccnt;
    alignas(64) std::atomic<double> m_weight;
    alignas(64) std::atomic<double> m_max_weight;
};

}
