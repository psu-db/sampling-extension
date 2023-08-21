/*
 * include/framework/WIRS.h
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
#include "ds/Alias.h"

namespace extension {

struct sample_state;
bool check_deleted(const record_t* record, sample_state *state);

struct wirs_node {
    struct wirs_node *left, *right;
    skey_t low, high;
    weight_t weight;
    psudb::Alias* alias;
};

struct ShardWIRSState {
    double tot_weight;
    std::vector<wirs_node*> nodes;
    psudb::Alias* top_level_alias;

    ~ShardWIRSState() {
        if (top_level_alias) delete top_level_alias;
    }
};

class ShardWIRS {
public:
    ShardWIRS(MutableBuffer* mbuffer, psudb::BloomFilter<skey_t>* bf, bool tagging)
    : m_reccnt(0), m_tombstone_cnt(0), m_deleted_cnt(0), m_total_weight(0), m_rejection_cnt(0), 
      m_ts_check_cnt(0), m_tagging(tagging), m_root(nullptr) {

        TIMER_INIT();
        TIMER_START();
        size_t len = mbuffer->get_record_count() * sizeof(record_t);
        m_data = (record_t*) psudb::sf_aligned_alloc(psudb::CACHELINE_SIZE, &len);

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

            if (bf && ptr->is_tombstone()) {
                m_tombstone_cnt++;
                bf->insert(ptr->key);
            }
        }
        TIMER_STOP();
        auto data_const_time = TIMER_RESULT();

        TIMER_START();
        if (m_reccnt > 0) {
            build_wirs_structure();
        }
        TIMER_STOP();
        auto ssi_const_time = TIMER_RESULT();

        #ifdef INSTRUMENT_MERGING
        fprintf(stderr, "buffer merge\t%ld\t%ld\t%ld\n", setup_time, data_const_time, ssi_const_time);
        #endif
    }

    ShardWIRS(ShardWIRS** shards, size_t len, psudb::BloomFilter<skey_t>* bf, bool tagging)
    : m_reccnt(0), m_tombstone_cnt(0), m_deleted_cnt(0), m_total_weight(0), m_rejection_cnt(0), m_ts_check_cnt(0), 
      m_tagging(tagging), m_root(nullptr) {

        TIMER_INIT();
        TIMER_START();
        std::vector<Cursor> cursors;
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

        auto attempted_len = attempted_reccnt * sizeof(record_t);
        m_data = (record_t *) psudb::sf_aligned_alloc(psudb::CACHELINE_SIZE, &attempted_len);
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
            advance_cursor(current);
            reccnt += 1;
        } while (reccnt < attempted_reccnt);

        TIMER_STOP();
        auto data_const_time = TIMER_RESULT();
        

        TIMER_START();
        if (m_reccnt > 0) {
            build_wirs_structure();
        }
        TIMER_STOP();
        auto ssi_const_time = TIMER_RESULT();
   }

    ~ShardWIRS() {
        if (m_data) free(m_data);
        for (size_t i=0; i<m_alias.size(); i++) {
            if (m_alias[i]) delete m_alias[i];
        }

        free_tree(m_root);
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

    void free_tree(struct wirs_node* node) {
        if (node) {
            delete node->alias;
            free_tree(node->left);
            free_tree(node->right);
            delete node;
        }
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

    // low - high -> decompose to a set of nodes.
    // Build psudb::Alias across the decomposed nodes.
    ShardWIRSState* get_shard_sample_state(const key_t& lower_key, const key_t& upper_key) {
        ShardWIRSState* res = new ShardWIRSState();
        //std::vector<struct wirs_node*> nodes;
        //double tot_weight = decompose_node(m_root, lower_key, upper_key, res->nodes);

        // Simulate a stack to unfold recursion.        
        double tot_weight = 0.0;
        struct wirs_node* st[64] = {0};
        st[0] = m_root;
        size_t top = 1;
        while(top > 0) {
            auto now = st[--top];
            if (covered_by(now, lower_key, upper_key) ||
                (now->left == nullptr && now->right == nullptr && intersects(now, lower_key, upper_key))) {
                res->nodes.emplace_back(now);
                tot_weight += now->weight;
            } else {
                if (now->left && intersects(now->left, lower_key, upper_key)) st[top++] = now->left;
                if (now->right && intersects(now->right, lower_key, upper_key)) st[top++] = now->right;
            }
        }
        
        //assert(tot_weight > 0.0);
        std::vector<double> weights;
        for (const auto& node: res->nodes) {
            weights.emplace_back(node->weight / tot_weight);
        }
        res->tot_weight = tot_weight;
        res->top_level_alias = new psudb::Alias(weights);

        return res;
    }

    // returns the number of records sampled
    // NOTE: This operation returns records strictly between the lower and upper bounds, not
    // including them.
    size_t get_samples(ShardWIRSState* shard_state, record_t *sample_set, const key_t& lower_key, const key_t& upper_key, size_t sample_sz, sample_state *state, gsl_rng *rng) {
        if (sample_sz == 0) {
            return 0;
        }

        // k -> sampling: three levels. 1. select a node -> select a fat point -> select a record.
        size_t cnt = 0;
        size_t attempts = 0;
        do {
            ++attempts;
            // first level....
            auto node = shard_state->nodes[shard_state->top_level_alias->get(rng)];
            // second level...
            auto fat_point = node->low + node->alias->get(rng);
            // third level...
            size_t rec_offset = fat_point * m_group_size + m_alias[fat_point]->get(rng);
            auto record = m_data + rec_offset;

            // bounds rejection
            if (lower_key > record->key || upper_key < record->key) {
                continue;
            // tombstone/delete rejection
            } else if (record->is_tombstone() || (state && check_deleted(record, state))) {
                continue;
            }

            sample_set[cnt++] = *record;
            
        } while (attempts < sample_sz);

        return cnt;
    }

    size_t get_lower_bound(const key_t& key) const {
        size_t min = 0;
        size_t max = m_reccnt - 1;

        const char * record_key;
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
    
private:
    bool covered_by(struct wirs_node* node, const key_t& lower_key, const key_t& upper_key) {
        auto low_index = node->low * m_group_size;
        auto high_index = std::min((node->high + 1) * m_group_size - 1, m_reccnt - 1);
        return lower_key < m_data[low_index].key && m_data[high_index].key < upper_key;
    }

    bool intersects(struct wirs_node* node, const key_t& lower_key, const key_t& upper_key) {
        auto low_index = node->low * m_group_size;
        auto high_index = std::min((node->high + 1) * m_group_size - 1, m_reccnt - 1);
        return lower_key < m_data[high_index].key && m_data[low_index].key < upper_key;
    }
    
    struct wirs_node* construct_wirs_node(const std::vector<weight_t>& weights, size_t low, size_t high) {
        if (low == high) {
            return new wirs_node{nullptr, nullptr, low, high, weights[low], new psudb::Alias({1.0})};
        } else if (low > high) return nullptr;

        std::vector<double> node_weights;
        weight_t sum = 0;
        for (size_t i = low; i < high; ++i) {
            node_weights.emplace_back(weights[i]);
            sum += weights[i];
        }

        for (auto& w: node_weights)
            if (sum) w /= sum;
            else w = 1.0 / node_weights.size();
        
        
        size_t mid = (low + high) / 2;
        return new wirs_node{construct_wirs_node(weights, low, mid),
                             construct_wirs_node(weights, mid + 1, high),
                             low, high, sum, new psudb::Alias(node_weights)};
    }

    void build_wirs_structure() {
        m_group_size = std::ceil(std::log(m_reccnt));
        size_t n_groups = std::ceil((double) m_reccnt / (double) m_group_size);
        
        // Fat point construction + low level alias....
        double sum_weight = 0.0;
        std::vector<weight_t> weights;
        std::vector<double> group_norm_weight;
        size_t i = 0;
        size_t group_no = 0;
        while (i < m_reccnt) {
            double group_weight = 0.0;
            group_norm_weight.clear();
            for (size_t k = 0; k < m_group_size && i < m_reccnt; ++k, ++i) {
                auto w = m_data[i].weight;
                group_norm_weight.emplace_back(w);
                group_weight += w;
                sum_weight += w;
            }

            for (auto& w: group_norm_weight)
                if (group_weight) w /= group_weight;
                else w = 1.0 / group_norm_weight.size();
            m_alias.emplace_back(new psudb::Alias(group_norm_weight));

            
            weights.emplace_back(group_weight);
        }

        assert(weights.size() == n_groups);

        m_root = construct_wirs_node(weights, 0, n_groups-1);
    }

    record_t* m_data;
    std::vector<psudb::Alias *> m_alias;
    wirs_node* m_root;

    size_t m_reccnt;
    size_t m_tombstone_cnt;
    size_t m_deleted_cnt;
    weight_t m_total_weight;
    size_t m_group_size;

    size_t m_rejection_cnt;
    size_t m_ts_check_cnt;

    bool m_tagging;
};

}
