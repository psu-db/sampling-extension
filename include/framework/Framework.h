/*
 * include/framework/Framework.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#pragma once

#include <atomic>
#include <numeric>

#include "framework/MutableBuffer.h"
#include "framework/Level.h"
#include "ds/Alias.h"

#include "util/timer.h"

namespace extension {

thread_local size_t sampling_attempts = 0;
thread_local size_t sampling_rejections = 0;
thread_local size_t deletion_rejections = 0;
thread_local size_t bounds_rejections = 0;
thread_local size_t tombstone_rejections = 0;

/*
 * thread_local size_t various_sampling_times go here.
 */
thread_local size_t mbuffer_alias_time = 0;
thread_local size_t alias_time = 0;
thread_local size_t alias_query_time = 0;
thread_local size_t rejection_check_time = 0;
thread_local size_t mbuffer_sample_time = 0;
thread_local size_t level_sample_time = 0;


/*
 * LSM Tree configuration global variables
 */

// True for mbuffer rejection sampling
static constexpr bool LSM_REJ_SAMPLE = true;

// True for leveling, false for tiering
static constexpr bool LSM_LEVELING = false;

// True for tagging, false for tombstones
static constexpr bool DELETE_TAGGING = true;

typedef ssize_t level_index;

class SamplingFramework;
struct sample_state {
    SamplingFramework *tree; 
    ShardId shid;
    char *buff;
    MutableBuffer *mbuffer;
    size_t mbuffer_cutoff;
};


class SamplingFramework {
    friend bool check_deleted();

public:
    SamplingFramework(size_t mbuffer_cap, size_t scale_factor, double max_tombstone_prop, 
                      double max_rejection_rate, gsl_rng *rng) 
        : 
          m_scale_factor(scale_factor), 
          m_max_tombstone_prop(max_tombstone_prop),
          m_max_rejection_rate(max_rejection_rate),
          m_last_level_idx(-1),
          m_buffer(new MutableBuffer(mbuffer_cap, LSM_REJ_SAMPLE, max_tombstone_prop*mbuffer_cap, rng))
    {}

    ~SamplingFramework() {
        delete m_buffer;

        for (size_t i=0; i<m_levels.size(); i++) {
            delete m_levels[i];
        }
    }

    int append(const key_t& key, const value_t& val, double weight, bool tombstone, gsl_rng *rng) {
        if (m_buffer->is_full()) {
            merge_mbuffer(rng);
        }

        return m_buffer->append(key, val, weight, tombstone);
    }


    int delete_record(const key_t& key, const value_t& val, gsl_rng *rng) {
        assert(DELETE_TAGGING);

        // Check the levels first. This assumes there aren't 
        // any undeleted duplicate records.
        for (auto level : m_levels) {
            if (level && level->delete_record(key, val)) {
                return 1;
            }
        }

        // the mbuffer will take the longest amount of time, and 
        // probably has the lowest probability of having the record,
        // so we'll check it last.
        return m_buffer->delete_record(key, val);
    }

    void range_sample(record_t *sample_set, size_t sample_sz, gsl_rng *rng) {
        TIMER_INIT();

        ds::Alias *mbuffer_alias = nullptr;
        size_t mbuffer_cutoff = 0;

        double mbuffer_weight = m_buffer->get_total_weight();

        TIMER_START();
        // Get the shard weights for each level. Index 0 is the mbuffer,
        // represented by nullptr.
        std::vector<std::pair<ShardId, ShardWSS *>> shards;
        shards.push_back({{-1, -1}, nullptr});

        std::vector<double> shard_weights;
        shard_weights.push_back(mbuffer_weight);

        for (auto &level : m_levels) {
            level->get_shard_weights(shard_weights, shards);
        }

        if (shard_weights.size() == 1 && shard_weights[0] == 0) {
            delete mbuffer_alias;
            return; // no records in the sampling range
        }
        double tot_weight = std::accumulate(shard_weights.begin(), shard_weights.end(), 0);
        for (auto& w: shard_weights) w /= tot_weight;

        // Construct alias structure
        auto alias = ds::Alias(shard_weights);
        TIMER_STOP();

        alias_time += TIMER_RESULT();

        std::vector<size_t> shard_samples(shard_weights.size(), 0);

        size_t rejections = sample_sz;
        size_t sample_idx = 0;

        sample_state state;
        state.tree = this;
        state.mbuffer = m_buffer;
        state.mbuffer_cutoff = mbuffer_cutoff;

        size_t passes = 0;
        do {
            TIMER_START();
            for (size_t i=0; i<rejections; i++) {
                shard_samples[alias.get(rng)] += 1;
            }
            TIMER_STOP();

            alias_query_time += TIMER_RESULT();

            rejections = 0;

            TIMER_START();
            while (shard_samples[0] > 0) {
                if (!LSM_REJ_SAMPLE && !mbuffer_alias) {
                    TIMER_START();
                    m_buffer->get_sample_range(&mbuffer_alias, &mbuffer_cutoff);
                    TIMER_STOP();

                    mbuffer_alias_time += TIMER_RESULT();
                }

                const record_t* rec;
                if (LSM_REJ_SAMPLE) {
                    rec = m_buffer->get_sample(rng);
                } else {
                    rec = m_buffer->get_record_at(mbuffer_alias->get(rng));
                }

                if (DELETE_TAGGING) {
                    if (rec && !rec->get_delete_status()) {
                        sample_set[sample_idx++] = *rec;
                    } else {
                        rejections++;
                    }
                } else {
                    if (rec && !m_buffer->check_tombstone(rec->key, rec->value)) {
                        sample_set[sample_idx++] = *rec;
                    } else {
                        rejections++;
                    }
                }

                shard_samples[0]--;
            }
            TIMER_STOP();
            mbuffer_sample_time += TIMER_RESULT();

            TIMER_START();
            for (size_t i=1; i<shard_samples.size(); i++) {
                // sample from each WIRS level
                state.shid = shards[i].first;
                auto sampled = shards[i].second->get_samples(sample_set + sample_idx, shard_samples[i], &state, rng) ;
                assert(sampled <= shard_samples[i]);
                sample_idx += sampled;
                rejections += shard_samples[i] - sampled;
                shard_samples[i] = 0;
            }
            TIMER_STOP();
            level_sample_time += TIMER_RESULT();
            passes++;
        } while (sample_idx < sample_sz);

        delete mbuffer_alias;

        enforce_rejection_rate(rng);
    }

    // Checks the tree and mbuffer for a tombstone corresponding to
    // the provided record in any shard *above* the shid, which
    // should correspond to the shard containing the record in question
    // 
    // Passing INVALID_shid indicates that the record exists within the MutableBuffer
    bool is_deleted(const record_t *record, const ShardId &shid, size_t mbuffer_cutoff) {

        TIMER_INIT();

        TIMER_START();
        // If tagging is enabled, we just need to check if the record has the delete tag set
        if (DELETE_TAGGING) {
            TIMER_STOP();
            rejection_check_time += TIMER_RESULT();
            return record->get_delete_status();
        }

        // Otherwise, we need to look for a tombstone.

        // check for tombstone in the mbuffer. This will require accounting for the cutoff eventually.
        if (m_buffer->check_tombstone(record->key, record->value)) {
            TIMER_STOP();
            rejection_check_time += TIMER_RESULT();
            return true;
        }

        // if the record is in the mbuffer, then we're done.
        if (shid == INVALID_shid) {
            TIMER_STOP();
            rejection_check_time += TIMER_RESULT();
            return false;
        }

        for (size_t lvl=0; lvl<=shid.level_idx; lvl++) {
            if (m_levels[lvl]->check_tombstone(0, record->key, record->value)) {
                TIMER_STOP();
                rejection_check_time += TIMER_RESULT();
                return true;
            }
        }

        // check the level containing the shard
        auto res = m_levels[shid.level_idx]->check_tombstone(shid.shard_idx + 1, record->key, record->value);

        TIMER_STOP();
        rejection_check_time += TIMER_RESULT();
        return res;
    }

    ShardWSS *create_static_structure() {
        std::vector<ShardWSS *> shards;

        if (m_levels.size() > 0) {
            for (int i=m_levels.size() - 1; i>= 0; i--) {
                if (m_levels[i]) {
                    shards.emplace_back(m_levels[i]->get_merged_shard());
                }
            }
        }

        shards.emplace_back(new ShardWSS(m_buffer, nullptr, DELETE_TAGGING));

        ShardWSS *shards_array[shards.size()];

        size_t j = 0;
        for (size_t i=0; i<shards.size(); i++) {
            if (shards[i]) {
                shards_array[j++] = shards[i];
            }
        }

        ShardWSS *flattened = new ShardWSS(shards_array, j, nullptr, DELETE_TAGGING);

        for (auto shard : shards) {
            delete shard;
        }

        return flattened;
    }


    size_t get_record_cnt() {
        size_t cnt = m_buffer->get_record_count();

        for (size_t i=0; i<m_levels.size(); i++) {
            if (m_levels[i]) cnt += m_levels[i]->get_record_cnt();
        }

        return cnt;
    }


    size_t get_tombstone_cnt() {
        size_t cnt = m_buffer->get_tombstone_count();

        for (size_t i=0; i<m_levels.size(); i++) {
            if (m_levels[i]) cnt += m_levels[i]->get_tombstone_count();
        }

        return cnt;
    }

    size_t get_height() {
        return m_levels.size(); // + m_disk_levels.size();
    }

    size_t get_memory_usage() {
        size_t cnt = m_buffer->get_memory_usage();

        for (size_t i=0; i<m_levels.size(); i++) {
            if (m_levels[i]) cnt += m_levels[i]->get_memory_usage();
        }

        return cnt;
    }

    size_t get_aux_memory_usage() {
        size_t cnt = m_buffer->get_aux_memory_usage();

        for (size_t i=0; i<m_levels.size(); i++) {
            if (m_levels[i]) cnt += m_levels[i]->get_aux_memory_usage();
        }

        return cnt;
    }

   bool validate_tombstone_proportion() {
        long double ts_prop;
        for (size_t i=0; i<m_levels.size(); i++) {
            if (m_levels[i]) {
                ts_prop = (long double) m_levels[i]->get_tombstone_count() / (long double) calc_level_record_capacity(i);
                if (ts_prop > (long double) m_max_tombstone_prop) {
                    return false;
                }
            }
        }

        return true;
    }

    size_t get_mbuffer_capacity() {
        return m_buffer->get_capacity();
    }

private:
    MutableBuffer *m_buffer;

    size_t m_scale_factor;
    double m_max_tombstone_prop;
    double m_max_rejection_rate;

    std::vector<Level *> m_levels;

    level_index m_last_level_idx;

    inline bool rejection(const record_t *record, ShardId shid, const key_t& lower_bound, const key_t& upper_bound, size_t mbuffer_cutoff) {
        if (record->is_tombstone()) {
            tombstone_rejections++;
            return true;
        } else if (record->key < lower_bound|| record->key > upper_bound) {
            bounds_rejections++;
            return true;
        } else if (is_deleted(record, shid, mbuffer_cutoff)) {
            deletion_rejections++;
            return true;
        }

        return false;
    }

    inline size_t shid_to_disk(ShardId shid) {
        return shid.level_idx - m_levels.size();
    }

    /*
     * Add a new level to the LSM Tree and return that level's index. Will
     * automatically determine whether the level should be on memory or on disk,
     * and act appropriately.
     */
    inline level_index grow() {
        level_index new_idx;

        size_t new_shard_cnt = (LSM_LEVELING) ? 1 : m_scale_factor;
        new_idx = m_levels.size();
        if (new_idx > 0) {
            assert(m_levels[new_idx - 1]->get_shard(0)->get_tombstone_count() == 0);
        }
        m_levels.emplace_back(new Level(new_idx, new_shard_cnt, DELETE_TAGGING));

        m_last_level_idx++;
        return new_idx;
    }


    // Merge the mutable buffer down into the tree, completing any required other
    // merges to make room for it.
    inline void merge_mbuffer(gsl_rng *rng) {
        TIMER_INIT();
        TIMER_START();

        if (!can_merge_with(0, m_buffer->get_record_count())) {
            merge_down(0, rng);
        }

        merge_mbuffer_into_l0(m_buffer, rng);
        TIMER_STOP();

        enforce_tombstone_maximum(0, rng);

        m_buffer->truncate();

        #ifdef INSTRUMENT_MERGING
            fprintf(stderr, "merge\t%ld\t%ld\n", TIMER_RESULT(), get_record_cnt());
        #endif
        return;
    }

    /*
     * Merge the specified level down into the tree. The level index must be
     * non-negative (i.e., this function cannot be used to merge the mbuffer). This
     * routine will recursively perform any necessary merges to make room for the 
     * specified level.
     */
    inline void merge_down(level_index idx, gsl_rng *rng) {
        level_index merge_base_level = find_mergable_level(idx);
        if (merge_base_level == -1) {
            merge_base_level = grow();
        }

        for (level_index i=merge_base_level; i>idx; i--) {
            merge_levels(i, i-1, rng);
            enforce_tombstone_maximum(i, rng);
        }

        return;
    }

    /*
     * Find the first level below the level indicated by idx that
     * is capable of sustaining a merge operation and return its
     * level index. If no such level exists, returns -1. Also
     * returns -1 if idx==0, and no such level exists, to simplify
     * the logic of the first merge.
     */
    inline level_index find_mergable_level(level_index idx) {

        if (idx == 0 && m_levels.size() == 0) return -1;

        bool level_found = false;
        bool disk_level;
        level_index merge_level_idx;

        size_t incoming_rec_cnt = get_level_record_count(idx);
        for (level_index i=idx+1; i<=m_last_level_idx; i++) {
            if (can_merge_with(i, incoming_rec_cnt)) {
                return i;
            }

            incoming_rec_cnt = get_level_record_count(i);
        }

        return -1;
    }

    /*
     * Merge the level specified by incoming level into the level specified
     * by base level. The two levels should be sequential--i.e. no levels
     * are skipped in the merge process--otherwise the tombstone ordering
     * invariant may be violated by the merge operation.
     */
    inline void merge_levels(level_index base_level, level_index incoming_level, gsl_rng *rng) {
        bool base_disk_level;
        bool incoming_disk_level;

        size_t base_idx = base_level; 
        size_t incoming_idx = incoming_level;

        // If the base level is a memory level, then the incoming level
        // cannot be a disk level.
        assert(!(!base_disk_level && incoming_disk_level));

        // merging two memory levels
        if (LSM_LEVELING) {
            auto tmp = m_levels[base_idx];
            m_levels[base_idx] = Level::merge_levels(m_levels[base_idx], m_levels[incoming_idx], DELETE_TAGGING, rng);
            mark_as_unused(tmp);
        } else {
            m_levels[base_idx]->append_merged_shards(m_levels[incoming_idx], rng);
        }

        mark_as_unused(m_levels[incoming_idx]);
        m_levels[incoming_idx] = new Level(incoming_level, (LSM_LEVELING) ? 1 : m_scale_factor, DELETE_TAGGING);
    }

    inline void merge_mbuffer_into_l0(MutableBuffer *mbuffer, gsl_rng *rng) {
        assert(m_levels[0]);
        if (LSM_LEVELING) {
            // FIXME: Kludgey implementation due to interface constraints.
            auto old_level = m_levels[0];
            auto temp_level = new Level(0, 1, DELETE_TAGGING);
            temp_level->append_buffer(mbuffer, rng);
            auto new_level = Level::merge_levels(old_level, temp_level, DELETE_TAGGING, rng);

            m_levels[0] = new_level;
            delete temp_level;
            mark_as_unused(old_level);
        } else {
            m_levels[0]->append_buffer(mbuffer, rng);
        }
    }

    /*
     * Mark a given memory level as no-longer in use by the tree. For now this
     * will just free the level. In future, this will be more complex as the
     * level may not be able to immediately be deleted, depending upon who
     * else is using it.
     */ 
    inline void mark_as_unused(Level *level) {
        delete level;
    }

    /*
     * Check the tombstone proportion for the specified level and
     * if the limit is exceeded, forcibly merge levels until all
     * levels below idx are below the limit.
     */
    inline void enforce_tombstone_maximum(level_index idx, gsl_rng *rng) {
        assert(m_levels[idx]);

        long double ts_prop = (long double) m_levels[idx]->get_tombstone_count() / (long double) calc_level_record_capacity(idx);

        if (ts_prop > (long double) m_max_tombstone_prop) {
            merge_down(idx, rng);
        }

        return;
    }

    inline void enforce_rejection_rate(gsl_rng *rng) {
        if (m_levels.size() == 0) {
            return;
        }

        for (size_t i=0; i<m_last_level_idx; i++) {
            if (m_levels[i]) {
                if (m_levels[i]->get_rejection_rate() > m_max_rejection_rate) {
                    merge_down(i, rng);
                }
            }
        } 
    }

    /*
     * Assume that level "0" should be larger than the mbuffer. The mbuffer
     * itself is index -1, which should return simply the mbuffer capacity.
     */
    inline size_t calc_level_record_capacity(level_index idx) {
        return m_buffer->get_capacity() * pow(m_scale_factor, idx+1);
    }

    /*
     * Returns the actual number of records present on a specified level. An
     * index value of -1 indicates the mutable buffer. Can optionally pass in
     * a pointer to the mutable buffer to use, if desired. Otherwise, there are
     * no guarantees about which mbuffer will be accessed if level_index is -1.
     */
    inline size_t get_level_record_count(level_index idx) {
        assert(idx >= -1);
        if (idx == -1) {
            return m_buffer->get_record_count();
        }

        return (m_levels[idx]) ? m_levels[idx]->get_record_cnt() : 0;
    }

    /*
     * Determines if the specific level can merge with another record containing
     * incoming_rec_cnt number of records. The provided level index should be 
     * non-negative (i.e., not refer to the mbuffer) and will be automatically
     * translated into the appropriate index into either the disk or memory level
     * vector.
     */
    inline bool can_merge_with(level_index idx, size_t incoming_rec_cnt) {
        bool disk_level;
        assert(idx >= 0);

        if (m_levels.size() <= idx || !m_levels[idx]) {
            return false;
        }

        if (LSM_LEVELING) {
            return m_levels[idx]->get_record_cnt() + incoming_rec_cnt <= calc_level_record_capacity(idx);
        } else {
            return m_levels[idx]->get_shard_count() < m_scale_factor;
        }

        // unreachable
        assert(true);
    }
};


bool check_deleted(const record_t* record, sample_state *state) {
    return state->tree->is_deleted(record, state->shid, state->mbuffer_cutoff);
}

}

