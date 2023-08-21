/*
 * include/framework/WSS.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#pragma once

#include <vector>
#include <cassert>
#include <queue>
#include <memory>

#include "framework/MutableBuffer.h"

#include "base/Cursor.h"
#include "util/timer.h"
#include "ds/Alias.h"

namespace extension {

struct sample_state;
bool check_deleted(const record_t *record, sample_state *state);

class ShardWSS {
public:

    ShardWSS(MutableBuffer* mbuffer, psudb::BloomFilter<skey_t>* bf, bool tagging)
    : m_reccnt(0), m_tombstone_cnt(0), m_rejection_cnt(0), m_ts_check_cnt(0), 
      m_deleted_cnt(0), m_total_weight(0), m_tagging(tagging) {
        TIMER_INIT();

        TIMER_START();
        std::vector<double> weights;
        weights.reserve(mbuffer->get_record_count());
        size_t alloc_len = mbuffer->get_record_count() * sizeof(record_t);
        m_data = (record_t*) psudb::sf_aligned_alloc(psudb::CACHELINE_SIZE, &alloc_len);

        auto base = mbuffer->sorted_output();
        auto stop = base + mbuffer->get_record_count();
        TIMER_STOP();

        auto setup_time = TIMER_RESULT();

        TIMER_START();
        for (auto ptr = base; ptr < stop; ptr++) {
            if (m_tagging && ptr->get_delete_status()) {
                continue;
            } else if (!(ptr->is_tombstone()) && (ptr + 1) < stop) {
                if (*ptr == *(ptr + 1) && (ptr + 1)->is_tombstone()) {
                    // we need to skip over both the record _and_ its tombstone,
                    // so we'll increment the pointer an extra time here.
                    ptr++;
                    continue;
                }
            }

            ptr->header &= 1;
            m_data[m_reccnt++] = *ptr;
            m_total_weight+= ptr->weight;
            weights.push_back(ptr->weight);

            if (bf && ptr->is_tombstone()) {
                m_tombstone_cnt++;
                bf->insert(ptr->key);
            }
        }
        TIMER_STOP();

        auto data_const_time = TIMER_RESULT();

        TIMER_START();
        // normalize the weights array
        for (size_t i=0; i<weights.size(); i++) {
            weights[i] = (double) weights[i] / (double) m_total_weight;
        }

        // build the alias structure
        m_alias = new psudb::Alias(weights);
        TIMER_STOP();

        auto ssi_const_time = TIMER_RESULT();

        #ifdef INSTRUMENT_MERGING
        fprintf(stderr, "buffer merge\t%ld\t%ld\t%ld\n", setup_time, data_const_time, ssi_const_time);
        #endif
    }

    ShardWSS(ShardWSS** shards, size_t len, psudb::BloomFilter<skey_t>* bf, bool tagging)
    : m_reccnt(0), m_tombstone_cnt(0), m_total_weight(0), m_rejection_cnt(0), m_ts_check_cnt(0), m_deleted_cnt(0), m_tagging(tagging) {

        TIMER_INIT();
        TIMER_START();
        std::vector<Cursor> cursors;
        std::vector<double> weights;
        cursors.reserve(len);

        size_t attempted_reccnt = 0;
        
        for (size_t i = 0; i < len; ++i) {
            if (shards[i]) {
                auto base = shards[i]->sorted_output();
                cursors.emplace_back(Cursor{base, base + shards[i]->get_record_count(), 0, shards[i]->get_record_count()});
                attempted_reccnt += shards[i]->get_record_count();
            } else {
                cursors.emplace_back(g_empty_cursor);
            }
        }

        auto alloc_len = attempted_reccnt * sizeof(record_t);
        m_data = (record_t*) psudb::sf_aligned_alloc(psudb::CACHELINE_SIZE, &alloc_len);

        weights.reserve(attempted_reccnt);
        TIMER_STOP();
        auto setup_time = TIMER_RESULT();

        TIMER_START();
        size_t reccnt = 0;
        Cursor *next = get_next(cursors);
        do {
            Cursor *current = next;
            next = get_next(cursors, current); 

            // Handle tombstone cancellation case if record tagging is not
            // enabled
            if (!m_tagging && !current->ptr->is_tombstone() && current != next &&
                next->ptr < next->end && *(current->ptr) == *(next->ptr) && 
                next->ptr->is_tombstone()) {
                
                reccnt += 2;
                advance_cursor(current);
                advance_cursor(next);
                next = get_next(cursors);
                continue;
            }

            // If the record is tagged as deleted, we can drop it
            if (current->ptr->get_delete_status()) {
                advance_cursor(current);
                reccnt += 1;
                continue;
            }

            // If the current record is a tombstone and tagging isn't in
            // use, add it to the Bloom filter and increment the tombstone
            // counter.
            if (!m_tagging && current->ptr->is_tombstone()) {
                ++m_tombstone_cnt;
                bf->insert(current->ptr->key);
            }

            // Copy the record into the new shard and increment the associated counters
            m_data[m_reccnt++] = *(current->ptr);
            m_total_weight += current->ptr->weight;
            weights.push_back((double) current->ptr->weight);
            advance_cursor(current);
            reccnt += 1;
        } while (reccnt < attempted_reccnt);
        TIMER_STOP();
        auto merge_time = TIMER_RESULT();

        TIMER_START();
        // normalize the weights array
        for (size_t i=0; i<weights.size(); i++) {
            weights[i] = weights[i] / (double) m_total_weight;
        }

        // build the alias structure
        m_alias = new psudb::Alias(weights);
        TIMER_STOP();

        auto alias_time = TIMER_RESULT();


        #ifdef INSTRUMENT_MERGING
        fprintf(stderr, "shard merge\t%ld\t%ld\t%ld\n", setup_time, merge_time, alias_time);
        #endif
   }

    ~ShardWSS() {
        if (m_data) free(m_data);
        if (m_alias) delete m_alias;
    }


    bool delete_record(const record_t &rec) {
        size_t idx = get_lower_bound(rec.key);
        if (idx >= m_reccnt) {
            return false;
        }

        while (idx < m_reccnt && m_data[idx] < rec) ++idx;

        if (idx < m_reccnt && m_data[idx] == rec) {
            m_data[idx].set_delete_status();
            m_deleted_cnt++;
            return true;
        }

        return false;
    }

    record_t* sorted_output() const {
        return m_data;
    }
    
    size_t get_record_count() const {
        return m_reccnt;
    }

    size_t get_tombstone_count() const {
        return m_tombstone_cnt;
    }

    const record_t* get_record_at(size_t idx) const {
        if (idx >= m_reccnt) return nullptr;
        return m_data + idx;
    }

    //
    // returns the number of records sampled
    // NOTE: This operation returns records strictly between the lower and upper bounds, not
    // including them.
    size_t get_samples(record_t *sample_set, size_t sample_sz, sample_state *state, gsl_rng *rng) {
        if (sample_sz == 0) {
            return 0;
        }

        size_t sampled_cnt=0;
        for (size_t i=0; i<sample_sz; i++) {
            size_t idx = m_alias->get(rng);

            if (m_data[idx].is_tombstone() || (state && check_deleted(m_data + idx, state))) {
                continue;
            }

            sample_set[sampled_cnt++] = m_data[idx];
        }

        return sampled_cnt;
    }

    size_t get_lower_bound(const skey_t& key) const {
        size_t min = 0;
        size_t max = m_reccnt - 1;

        while (min < max) {
            size_t mid = (min + max) / 2;

            if (key > m_data[mid].key) {
                min = mid + 1;
            } else {
                max = mid;
            }
        }

        return min;
    }

    bool check_tombstone(const record_t &rec) {
        m_ts_check_cnt += 1;

        size_t idx = get_lower_bound(rec.key);
        if (idx >= m_reccnt) {
            return false;
        }

        while (idx < m_reccnt && m_data[idx] < rec) ++idx;

        if (idx >= m_reccnt) {
            return false;
        }

        bool result = m_data[idx] == rec;

        if (result && m_data[idx].is_tombstone()) {
            m_rejection_cnt++;
        }
        return result;
    }


    size_t get_memory_usage() {
        return 0;
    }


    weight_t get_total_weight() {
        return m_total_weight;
    }

    size_t get_rejection_count() {
        return m_rejection_cnt;
    }

    size_t get_ts_check_count() {
        return m_ts_check_cnt;
    }

    size_t get_deleted_count() {
        assert(m_tagging);
        return m_deleted_cnt;
    }
    
private:
    record_t* m_data;
    psudb::Alias *m_alias;

    size_t m_reccnt;
    size_t m_tombstone_cnt;
    size_t m_deleted_cnt;
    weight_t m_total_weight;

    size_t m_rejection_cnt;
    size_t m_ts_check_cnt;

    bool m_tagging;
};

}
