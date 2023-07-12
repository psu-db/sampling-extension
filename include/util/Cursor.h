/*
 * include/util/Cursor.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#pragma once

#include <vector>

#include "util/base.h"
#include "util/record.h"

namespace extension {
struct Cursor {
    record_t *ptr;
    const record_t *end;
    size_t cur_rec_idx;
    size_t rec_cnt;

    friend bool operator==(const Cursor &a, const Cursor &b) {
        return a.ptr == b.ptr && a.end == b.end;
    }
};

static Cursor g_empty_cursor = {0};

/*
 * Advance the cursor to the next record. If the cursor is backed by an
 * iterator, will attempt to advance the iterator once the cursor reaches its
 * end and reset the cursor to the beginning of the read page.
 *
 * If the advance succeeds, ptr will be updated to point to the new record
 * and true will be returned. If the advance reaches the end, then ptr will
 * be updated to be equal to end, and false will be returned. Iterators will
 * not be closed.
 */
inline static bool advance_cursor(Cursor *cur) {
    cur->ptr++;
    cur->cur_rec_idx++;

    if (cur->cur_rec_idx >= cur->rec_cnt) return false;

    return true;
}

/*
 *   Process the list of cursors to return the cursor containing the next
 *   largest element. Does not advance any of the cursors. If current is
 *   specified, then skip the current head of that cursor during checking. 
 *   This allows for "peaking" at the next largest element after the current 
 *   largest is processed.
 */
inline static Cursor *get_next(std::vector<Cursor> &cursors, Cursor *current=&g_empty_cursor) {
    const record_t *min_rec = nullptr;
    Cursor *result = &g_empty_cursor;
    for (size_t i=0; i< cursors.size(); i++) {
        if (cursors[i] == g_empty_cursor) continue;

        const record_t *rec = (&cursors[i] == current) ? cursors[i].ptr + 1 : cursors[i].ptr;
        if (rec >= cursors[i].end) continue;

        if (min_rec == nullptr) {
            result = &cursors[i];
            min_rec = rec;
            continue;
        }

        if (*rec < *min_rec) {
            result = &cursors[i];
            min_rec = rec;
        }
    }

    return result;
} 

}
