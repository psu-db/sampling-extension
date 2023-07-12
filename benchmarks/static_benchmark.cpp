/*
 * benchmarks/static_benchmark.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#include "bench.h"

size_t g_insert_batch_size = 1000;

static void sample_benchmark(extension::ShardWSS *shard, size_t k, size_t trial_cnt)
{
    char progbuf[25];
    sprintf(progbuf, "sampling (%ld):", k);

    size_t batch_size = 100;
    size_t batches = trial_cnt / batch_size;
    size_t total_time = 0;

    extension::record_t sample_set[k];

    for (int i=0; i<batches; i++) {
        progress_update((double) (i * batch_size) / (double) trial_cnt, progbuf);
        auto start = std::chrono::high_resolution_clock::now();
        for (int j=0; j < batch_size; j++) {
            shard->get_samples(sample_set, k, nullptr, g_rng);
        }
        auto stop = std::chrono::high_resolution_clock::now();

        total_time += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    }

    progress_update(1.0, progbuf);

    size_t throughput = (((double)(trial_cnt * k) / (double) total_time) * 1e9);

    fprintf(stdout, "%.0ld\n", throughput);
}


int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: static_throughput <filename> <record_count> [osm_data]\n");
        exit(EXIT_FAILURE);
    }

    std::string filename = std::string(argv[1]);
    size_t record_count = atol(argv[2]);
    bool use_osm = (argc == 4) ? atoi(argv[3]) : 0;

    init_bench_env(record_count, true, use_osm);

    auto de_wss = extension::SamplingFramework(12000, 6, 1, 100, g_rng);

    std::fstream datafile;
    datafile.open(filename, std::ios::in);

    // warm up the tree with initial_insertions number of initially inserted records
    warmup(&datafile, &de_wss, record_count, 0);

    auto static_shard = de_wss.create_static_structure();

    sample_benchmark(static_shard, 1000, 10000);

    delete_bench_env();
    fflush(stdout);
    fflush(stderr);

    delete static_shard;

    exit(EXIT_SUCCESS);
}
