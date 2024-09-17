/*
 * include/framework/Level.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Revised BSD License.
 *
 */
#pragma once

#include <vector>
#include <memory>

#include "base/types.h"
#include "base/bf_config.h"
#include "framework/WSS.h"
#include "ds/BloomFilter.h"

namespace extension {


class Level {

static const size_t REJECTION_TRIGGER_THRESHOLD = 1024;

private:
    struct InternalLevelStructure {
        InternalLevelStructure(size_t cap)
        : m_cap(cap)
        , m_shards(new ShardWSS*[cap]{nullptr})
        , m_bfs(new psudb::BloomFilter<skey_t>*[cap]{nullptr}) {} 

        ~InternalLevelStructure() {
            for (size_t i = 0; i < m_cap; ++i) {
                if (m_shards[i]) delete m_shards[i];
                if (m_bfs[i]) delete m_bfs[i];
            }

            delete[] m_shards;
            delete[] m_bfs;
        }

        size_t m_cap;
        ShardWSS** m_shards;
        psudb::BloomFilter<skey_t>** m_bfs;
    };

public:
    Level(ssize_t level_no, size_t shard_cap, bool tagging)
    : m_level_no(level_no), m_shard_cnt(0)
    , m_structure(new InternalLevelStructure(shard_cap))
    , m_tagging(tagging) {}

    // Create a new memory level sharing the shards and repurposing it as previous level_no + 1
    // WARNING: for leveling only.
    Level(Level* level, bool tagging)
    : m_level_no(level->m_level_no + 1), m_shard_cnt(level->m_shard_cnt)
    , m_structure(level->m_structure)
    , m_tagging(tagging) {
        assert(m_structure->m_cap == 1 && m_shard_cnt == 1);
    }


    ~Level() {}

    // WARNING: for leveling only.
    // assuming the base level is the level new level is merging into. (base_level is larger.)
    static Level* merge_levels(Level* base_level, Level* new_level, bool tagging) {
        assert(base_level->m_level_no > new_level->m_level_no || (base_level->m_level_no == 0 && new_level->m_level_no == 0));
        auto res = new Level(base_level->m_level_no, 1, tagging);
        res->m_shard_cnt = 1;
        res->m_structure->m_bfs[0] =
            new psudb::BloomFilter<skey_t>(BF_FPR,
                            new_level->get_tombstone_count() + base_level->get_tombstone_count(),
                            BF_HASH_FUNCS);
        ShardWSS* shards[2];
        shards[0] = base_level->m_structure->m_shards[0];
        shards[1] = new_level->m_structure->m_shards[0];

        res->m_structure->m_shards[0] = new ShardWSS(shards, 2, res->m_structure->m_bfs[0], tagging);
        return res;
    }

    void append_buffer(MutableBuffer* mbuffer) {
        assert(m_shard_cnt < m_structure->m_cap);
        m_structure->m_bfs[m_shard_cnt] = new psudb::BloomFilter<skey_t>(BF_FPR, mbuffer->get_tombstone_count(), BF_HASH_FUNCS);
        m_structure->m_shards[m_shard_cnt] = new ShardWSS(mbuffer, m_structure->m_bfs[m_shard_cnt], m_tagging);
        ++m_shard_cnt;
    }

    void append_merged_shards(Level* level) {
        assert(m_shard_cnt < m_structure->m_cap);
        m_structure->m_bfs[m_shard_cnt] = new psudb::BloomFilter<skey_t>(BF_FPR, level->get_tombstone_count(), BF_HASH_FUNCS);
        m_structure->m_shards[m_shard_cnt] = new ShardWSS(level->m_structure->m_shards, level->m_shard_cnt, m_structure->m_bfs[m_shard_cnt], m_tagging);
        ++m_shard_cnt;
    }

    ShardWSS *get_merged_shard() {
        ShardWSS *shards[m_shard_cnt];

        for (size_t i=0; i<m_shard_cnt; i++) {
            shards[i] = (m_structure->m_shards[i]) ? m_structure->m_shards[i] : nullptr;
        }

        return new ShardWSS(shards, m_shard_cnt, nullptr, m_tagging);
    }

    // Append the sample range in-order
    void get_shard_weights(std::vector<double>& weights, std::vector<std::pair<ShardId, ShardWSS *>> &shards) {
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_shards[i] && m_structure->m_shards[i]->get_total_weight() > 0) {
                shards.push_back({{m_level_no, (ssize_t) i}, m_structure->m_shards[i]});
                weights.push_back(m_structure->m_shards[i]->get_total_weight());
            }
        }
    }

    bool bf_rejection_check(size_t shard_stop, const key_t& key) {
        for (size_t i = 0; i < shard_stop; ++i) {
            if (m_structure->m_bfs[i] && m_structure->m_bfs[i]->lookup(key))
                return true;
        }
        return false;
    }

    bool check_tombstone(size_t shard_stop, const record_t &rec) {
        if (m_shard_cnt == 0) return false;

        for (int i = m_shard_cnt - 1; i >= (ssize_t) shard_stop;  i--) {
            if (m_structure->m_shards[i] && (m_structure->m_bfs[i]->lookup(rec.key))
                && m_structure->m_shards[i]->check_tombstone(rec))
                return true;
        }
        return false;
    }

    bool delete_record(const record_t &rec) {
        for (size_t i = 0; i < m_structure->m_cap;  ++i) {
            if (m_structure->m_shards[i] && m_structure->m_shards[i]->delete_record(rec)) {
                return true;
            }
        }

        return false;
    }

    const record_t* get_record_at(size_t shard_no, size_t idx) {
        return m_structure->m_shards[shard_no]->get_record_at(idx);
    }
    
    ShardWSS* get_shard(size_t idx) {
        return m_structure->m_shards[idx];
    }

    size_t get_shard_count() {
        return m_shard_cnt;
    }

    size_t get_record_cnt() {
        size_t cnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            cnt += m_structure->m_shards[i]->get_record_count();
        }

        return cnt;
    }
    
    size_t get_tombstone_count() {
        size_t res = 0;
        for (size_t i = 0; i < m_shard_cnt; ++i) {
            res += m_structure->m_shards[i]->get_tombstone_count();
        }
        return res;
    }

    size_t get_aux_memory_usage() {
        size_t cnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_bfs[i]) {
                cnt += m_structure->m_bfs[i]->memory_usage();
            }
        }

        return cnt;
    }

    size_t get_memory_usage() {
        size_t cnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_shards[i]) {
                cnt += m_structure->m_shards[i]->get_memory_usage();
            }
        }

        return cnt;
    }

    double get_tombstone_prop() {
        size_t tscnt = 0;
        size_t reccnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_shards[i]) {
                tscnt += m_structure->m_shards[i]->get_tombstone_count();
                reccnt += m_structure->m_shards[i]->get_record_count();
            }
        }

        return (double) tscnt / (double) (tscnt + reccnt);
    }

    size_t get_rejection_count() {
        size_t rej_cnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_shards[i]) {
                rej_cnt += m_structure->m_shards[i]->get_rejection_count();
            }
        }

        return rej_cnt;
    }

    double get_rejection_rate() {
        size_t rej_cnt = 0;
        size_t attempt_cnt = 0;
        for (size_t i=0; i<m_shard_cnt; i++) {
            if (m_structure->m_shards[i]) {
                attempt_cnt += m_structure->m_shards[i]->get_ts_check_count();
                rej_cnt += m_structure->m_shards[i]->get_rejection_count();
            }
        }

        if (attempt_cnt == 0) return 0;

        // the rejection rate is considered 0 until we exceed an
        // absolute threshold of rejections.
        if (rej_cnt <= REJECTION_TRIGGER_THRESHOLD) return 0;

        return (double) rej_cnt / (double) attempt_cnt;
    }

private:
    ssize_t m_level_no;
    size_t m_shard_cnt;
    size_t m_shard_size_cap;
    bool m_tagging;

    std::shared_ptr<InternalLevelStructure> m_structure;
};

}
