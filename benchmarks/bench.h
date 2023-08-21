/*
 * benchmarks/bench.h
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#pragma once

#include "framework/Framework.h"

#include <ds/BTree.h>

#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <set>
#include <string>
#include <random>


typedef std::pair<extension::skey_t, extension::value_t> btree_record;
struct btree_key_extract {
    static const extension::skey_t &get(const btree_record &v) {
        return v.first;
    }
};

typedef psudb::BTree<extension::skey_t, btree_record, btree_key_extract> TreeMap;

typedef struct record {
    extension::skey_t key;
    extension::value_t value;
    extension::weight_t weight;

    friend bool operator<(const struct record &first, const struct record &other) {
        return (first.key < other.key) || (first.key == other.key && first.value < other.value);
    }

    friend bool operator==(const struct record &first, const struct record &other) {
        return (first.key == other.key) && (first.value == other.value);
    }
} record;

typedef std::pair<extension::skey_t, extension::skey_t> key_range;

static gsl_rng *g_rng;
static std::set<record> *g_to_delete;
static bool g_osm_data;

static extension::skey_t g_min_key = UINT64_MAX;
static extension::skey_t g_max_key = 0;

static size_t g_max_record_cnt = 0;
static size_t g_reccnt = 0;

static constexpr unsigned int DEFAULT_SEED = 0;

typedef enum Operation {
    READ,
    WRITE,
    DELETE
} Operation;

static unsigned int get_random_seed()
{
    unsigned int seed = 0;
    std::fstream urandom;
    urandom.open("/dev/urandom", std::ios::in|std::ios::binary);
    urandom.read((char *) &seed, sizeof(seed));
    urandom.close();

    return seed;
}


static extension::skey_t osm_to_key(const char *key_field) {
    double tmp_key = (atof(key_field) + 180) * 10e6;
    return (extension::skey_t) tmp_key;
}


static void init_bench_rng(unsigned int seed, const gsl_rng_type *type) 
{
    g_rng = gsl_rng_alloc(type);
    gsl_rng_set(g_rng, seed);
}


static void init_bench_env(size_t max_reccnt, bool random_seed, bool osm_correction=true)
{
    unsigned int seed = (random_seed) ? get_random_seed() : DEFAULT_SEED;
    init_bench_rng(seed, gsl_rng_mt19937);
    g_to_delete = new std::set<record>();
    g_osm_data = osm_correction;
    g_max_record_cnt = max_reccnt;
    g_reccnt = 0;
}


static void delete_bench_env()
{
    gsl_rng_free(g_rng);
    delete g_to_delete;
}

static bool next_record(std::fstream *file, extension::skey_t* key, extension::value_t* val, extension::weight_t *weight)
{
    if (g_reccnt >= g_max_record_cnt) return false;

    std::string line;
    if (std::getline(*file, line, '\n')) {
        std::stringstream line_stream(line);
        std::string key_field;
        std::string value_field;
        std::string weight_field;

        std::getline(line_stream, value_field, '\t');
        std::getline(line_stream, key_field, '\t');
        std::getline(line_stream, weight_field, '\t');

        *key = (g_osm_data) ? osm_to_key(key_field.c_str()) : atol(key_field.c_str());
        *val = atol(value_field.c_str());
        *weight = atof(weight_field.c_str());

        if (*key < g_min_key) g_min_key = *key;

        if (*key > g_max_key) g_max_key = *key;

        g_reccnt++;

        return true;
    }

    return false;
}


static bool build_insert_vec(std::fstream *file, std::vector<record> &vec, size_t n) {
    vec.clear();
    for (size_t i=0; i<n; i++) {
        record rec;
        if (!next_record(file, &rec.key, &rec.value, &rec.weight)) {
            if (i == 0) {
                return false;
            }

            break;
        }

        vec.emplace_back((record){rec.key, rec.value, rec.weight});
    }

    return true;
}


static bool build_btree_insert_vec(std::fstream *file, std::vector<std::pair<btree_record, extension::weight_t>> &vec, size_t n)
{
    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    vec.clear();
    for (size_t i=0; i<n; i++) {
        if (!next_record(file, &key, &val, &weight)) {
            if (i == 0) {
                return false;
            }

            break;
        }

        btree_record rec = {key, val};
        vec.push_back({rec, weight});
    }

    return true;

}

static bool build_avl_insert_vec(std::fstream *file, std::vector<std::pair<extension::skey_t, extension::weight_t>> &vec, size_t n)
{
    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    vec.clear();
    for (size_t i=0; i<n; i++) {
        if (!next_record(file, &key, &val, &weight)) {
            if (i == 0) {
                return false;
            }

            break;
        }

        vec.push_back({key, weight});
    }

    return true;

}
/*
 * helper routines for displaying progress bars to stderr
 */
static const char *g_prog_bar = "======================================================================";
static const size_t g_prog_width = 50;

static void progress_update(double percentage, std::string prompt) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * g_prog_width);
    int rpad = (int) (g_prog_width - lpad);
    fprintf(stderr, "\r(%3d%%) %20s [%.*s%*s]", val, prompt.c_str(), lpad, g_prog_bar, rpad, "");
    fflush(stderr);   

    if (percentage >= 1) fprintf(stderr, "\n");
}



static bool warmup(std::fstream *file, extension::SamplingFramework *lsmtree, size_t count, double delete_prop, bool progress=true)
{
    std::string line;

    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    size_t del_buf_size = 10000;
    size_t del_buf_ptr = del_buf_size;

    extension::record_t delbuf[del_buf_size];
    std::set<record> deleted_keys;

    size_t inserted = 0;
    
    double last_percent = 0;
    for (size_t i=0; i<count; i++) {
        if (!next_record(file, &key, &val, &weight)) {
            return false;
        }

        inserted++;
        lsmtree->append(key, val, weight, false, g_rng);

        if (i > lsmtree->get_buffer_capacity() && del_buf_ptr >= del_buf_size) {
            lsmtree->range_sample(delbuf, del_buf_size, g_rng);
            del_buf_ptr = 0;
            deleted_keys.clear();
        }

        if (i > lsmtree->get_buffer_capacity() && gsl_rng_uniform(g_rng) < delete_prop) {
            auto key = delbuf[del_buf_ptr].key;
            auto val = delbuf[del_buf_ptr].value;
            del_buf_ptr++;

            if (deleted_keys.find({key, val, weight}) == deleted_keys.end()) {
                if (extension::DELETE_POLICY) {
                    lsmtree->delete_record(key, val, g_rng);
                } else {
                    lsmtree->append(key, val, weight, true, g_rng);
                }
                deleted_keys.insert({key, val, weight});
            }
        }

        if (progress && ((double) i / (double) count) - last_percent > .01) {
            progress_update((double) i / (double) count, "warming up:");
            last_percent = (double) i / (double) count;
        }
    }

    if (progress) {
        progress_update(1, "warming up:");
    }

    return true;
}

static bool warmup(std::fstream *file, TreeMap *tree, size_t count, double delete_prop, bool progress=true)
{
    std::string line;

    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    size_t del_buf_size = 100;
    size_t del_buf_idx = del_buf_size;
    std::vector<extension::skey_t> delbuf(del_buf_size);

    
    double last_percent = 0;
    for (size_t i=0; i<count; i++) {
        if (!next_record(file, &key, &val, &weight)) {
            return false;
        }

        tree->insert({key, val}, weight);

        if ( i > 10*del_buf_size && del_buf_idx == del_buf_size) {
            tree->range_sample(g_min_key, g_max_key, del_buf_size, delbuf, g_rng);
            del_buf_idx = 0;
        }

        if (del_buf_idx != del_buf_size && gsl_rng_uniform(g_rng) < delete_prop) {
            tree->erase_one(delbuf[del_buf_idx++]);
        }

        if (progress && ((double) i / (double) count) - last_percent > .01) {
            progress_update((double) i / (double) count, "warming up: ");
            last_percent = (double) i / (double) count;
        }
    }

    if (progress) {
        progress_update(1, "warming up:");
    }

    return true;
}


static key_range get_key_range(extension::skey_t min, extension::skey_t max, double selectivity)
{
    size_t range_length = (max - min) * selectivity;

    extension::skey_t max_bottom = max - range_length;
    extension::skey_t bottom;

    while ((bottom = gsl_rng_get(g_rng)) > range_length) 
        ;

    return {min + bottom, min + bottom + range_length};
}


static void reset_extension_perf_metrics() {
    /*
     * rejection counters are zeroed automatically by the
     * sampling function itself.
     */

    extension::buffer_alias_time = 0;
    extension::alias_time = 0;
    extension::alias_query_time = 0;
    extension::buffer_sample_time = 0;
    extension::level_sample_time = 0;
    extension::rejection_check_time = 0;
}


static void build_lsm_tree(extension::SamplingFramework *tree, std::fstream *file) {

    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    size_t i=0;
    while (next_record(file, &key, &val, &weight)) {
        auto res = tree->append(key, val, weight, false, g_rng);
        assert(res);
    }
}

static void build_btree(TreeMap *tree, std::fstream *file) {
    // for looking at the insert time distribution
    extension::skey_t key;
    extension::value_t val;
    extension::weight_t weight;

    size_t i=0;
    while (next_record(file, &key, &val, &weight)) {
        auto res = tree->insert({key, val}, weight);
        assert(res.second);
    }
}
