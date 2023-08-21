/*
 * include/util/types.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu>
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 * A centralized header file for various datatypes used throughout the
 * code base. There are a few very specific types, such as header formats,
 * that are defined within the header files that make direct use of them,
 * but all generally usable, simple types are defined here.
 */
#pragma once

#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <string>

#include "util/alignment.h"

namespace extension {

using std::byte;

// Represents a page offset within a specific file (physical or virtual)
typedef uint32_t PageNum;
const size_t PAGE_MAX = UINT32_MAX;

// Byte offset within a page. Also used for lengths of records, etc.,
// within the codebase. size_t isn't necessary, as the maximum offset
// is only parm::PAGE_SIZE 
typedef uint16_t PageOffset;
const size_t PAGE_OFFSET_MAX = UINT16_MAX;

// A unique identifier for a frame within a buffer or cache.
typedef int32_t FrameId;
const size_t FRAMEID_MAX = UINT32_MAX;

// A unique timestamp for use in MVCC concurrency control.
typedef uint32_t Timestamp;
const Timestamp TIMESTAMP_MIN = 0;
const Timestamp TIMESTAMP_MAX = UINT32_MAX;

// Invalid values for various IDs. Used throughout the code base to indicate
// uninitialized values and error conditions.
const PageNum INVALID_PNUM = 0;
const FrameId INVALID_FRID = -1;

// An ID for a given shard within the tree. The level_idx is the index
// in the memory_levels and disk_levels vectors corresponding to the
// shard, and the shard_idx is the index with the level (always 0 in the
// case of leveling) Note that the two vectors of levels are treated
// as a contiguous index space, so index 0-memory_levels.size() corresponds
// to a memory level, and memory_levels.size()-memory_levels.size() +
// disk_levels.size() corresponds to a disk level
struct ShardId {
    ssize_t level_idx;
    ssize_t shard_idx;

    friend bool operator==(const ShardId &shid1, const ShardId &shid2) {
        return shid1.level_idx == shid2.level_idx && shid1.shard_idx == shid2.shard_idx;
    }
};

const ShardId INVALID_SHID = {-1, -1};

/*
 * Record component types--adjust these to control data sizes, etc.
 */

/* 
 * Record header information--used primarily for tombstone/delete bits
 * right now.
 */
typedef uint32_t hdr_t;    

/*
 * Sampling key type. key_t is not used anymore due to confusing
 * conflict with the POSIX type of the same name.
 */
typedef uint64_t skey_t;    

/*
 * Record value type.
 */
typedef uint32_t value_t;

/*
 * Weight type for weighted sampling applications.
 */
typedef uint64_t weight_t;

struct SampleRange {
    ShardId shid;
    skey_t low;
    skey_t high;
};



}
