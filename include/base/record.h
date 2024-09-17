/*
 * include/util/record.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Revised BSD License.
 *
 */
#pragma once
#pragma once

#include "util/alignment.h"
#include "base/types.h"

#include <cstring>

namespace extension {

struct record_t {
    skey_t key;
    value_t value;
    hdr_t header;

    // the weight field can be omitted if unnecessary
    #ifdef WEIGHTED_SAMPLING
    weight_t weight;
    #endif

    inline void set_delete_status() {
        header |= 2;
    }

    inline bool get_delete_status() const {
        return header & 2;
    }

    inline bool is_tombstone() const {
        return header & 1;
    }

    inline bool operator<(const record_t& other) const {
        return key < other.key || (key == other.key && value < other.value);
    }

    inline bool operator==(const record_t& other) const {
        return key == other.key && value == other.value;
    }
};

static bool buffer_record_cmp(const record_t& a, const record_t& b) {
    return (a.key < b.key) || (a.key == b.key && a.value < b.value)
        || (a.key == b.key && a.value == b.value && a.header < b.header);
}

}
