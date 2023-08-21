/*
 * benchmarks/btree_benchmark.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#include "bench.h"

size_t g_insert_batch_size = 1000;

static bool insert_benchmark(TreeMap *tree, std::fstream *file, 
                      size_t insert_cnt, double delete_prop) {
    size_t deletes = insert_cnt * delete_prop;
    size_t delete_batch_size = g_insert_batch_size * delete_prop * 15;
    size_t delete_idx = delete_batch_size;

    std::vector<extension::skey_t> delbuf;
    tree->range_sample(g_min_key, g_max_key, delete_batch_size, delbuf, g_rng);
    size_t applied_deletes = 0;

    std::set<extension::skey_t> deleted;

    size_t applied_inserts = 0;
    std::vector<std::pair<btree_record, extension::weight_t>> insert_vec;
    insert_vec.reserve(g_insert_batch_size);
    bool continue_benchmark = true;

    size_t total_time = 0;

    while (applied_inserts < insert_cnt && continue_benchmark) { 
        continue_benchmark = build_btree_insert_vec(file, insert_vec, g_insert_batch_size);
        if (insert_vec.size() == 0) {
            break;
        }

        // if we've fully processed the delete vector, sample a new
        // set of records to delete.
        if (delete_idx >= delete_batch_size) {
            tree->range_sample(g_min_key, g_max_key, delete_batch_size, delbuf, g_rng);
            deleted.clear();
        }

        progress_update((double) applied_inserts / (double) insert_cnt, "inserting:");
        size_t local_inserted = 0;
        size_t local_deleted = 0;
        
        auto insert_start = std::chrono::high_resolution_clock::now();
        for (size_t i=0; i<insert_vec.size(); i++) {
            // process a delete if necessary
            if (applied_deletes < deletes && delete_idx < delete_batch_size && gsl_rng_uniform(g_rng) < delete_prop) {

                if (deleted.find(delbuf[delete_idx]) == deleted.end()) {
                    tree->erase_one(delbuf[delete_idx]);
                    local_deleted++;
                    deleted.insert(delbuf[delete_idx]);
                }
                delete_idx++;
            }
            // insert the record;
            tree->insert(insert_vec[local_inserted].first, insert_vec[local_inserted].second);
            local_inserted++;
        }
        auto insert_stop = std::chrono::high_resolution_clock::now();

        applied_deletes += local_deleted;
        applied_inserts += local_inserted;
        total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(insert_stop - insert_start).count();
    } 

    progress_update(1.0, "inserting:");
    size_t throughput = (((double) (applied_inserts + applied_deletes) / (double) total_time) * 1e9);

    fprintf(stdout, "%ld\n", throughput);

    return continue_benchmark;
}


static void sample_benchmark(TreeMap *tree, size_t k, size_t trial_cnt)
{
    char progbuf[25];
    sprintf(progbuf, "sampling (%ld):", k);

    size_t batch_size = 100;
    size_t batches = trial_cnt / batch_size;
    size_t total_time = 0;

    std::vector<extension::skey_t> sample_set;
    sample_set.reserve(k);

    for (int i=0; i<batches; i++) {
        progress_update((double) (i * batch_size) / (double) trial_cnt, progbuf);
        auto start = std::chrono::high_resolution_clock::now();
        for (int j=0; j < batch_size; j++) {
            tree->range_sample(g_min_key, g_max_key, k, sample_set, g_rng);
        }
        auto stop = std::chrono::high_resolution_clock::now();

        total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    }

    progress_update(1.0, progbuf);

    size_t throughput = (((double)(trial_cnt * k) / (double) total_time) * 1e9);

    fprintf(stdout, "%zu %.0ld\n", k, throughput);
}


int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: btree_benchmark <filename> <record_count> <delete_proportion> [osm_data]\n");
        exit(EXIT_FAILURE);
    }

    std::string filename = std::string(argv[1]);
    size_t record_count = atol(argv[2]);
    double delete_prop = atof(argv[3]);
    bool use_osm = (argc == 5) ? atoi(argv[4]) : 0;

    double insert_batch = 0.1; 

    init_bench_env(record_count, true, use_osm);

    auto sampling_tree = TreeMap();

    std::fstream datafile;
    datafile.open(filename, std::ios::in);

    // warm up the tree with initial_insertions number of initially inserted
    // records
    size_t warmup_cnt = insert_batch * record_count;
    warmup(&datafile, &sampling_tree, warmup_cnt, delete_prop);

    size_t insert_cnt = record_count - warmup_cnt;

    insert_benchmark(&sampling_tree, &datafile, insert_cnt, delete_prop);

    size_t max_sample_size = 1000000;
    for (size_t sample_size =1; sample_size < max_sample_size; sample_size *= 10) {
        sample_benchmark(&sampling_tree, sample_size, 10000);
    }

    delete_bench_env();
    fflush(stdout);
    fflush(stderr);

    exit(EXIT_SUCCESS);
}
