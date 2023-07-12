/*
 * tests/shard_tests.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#include <check.h>

#include "framework/WSS.h"
#include "framework/Framework.h"
#include "framework/Level.h"
#include "util/bf_config.h"

using namespace extension;
using namespace ds;

bool roughly_equal(int n1, int n2, size_t mag, double epsilon) {
    return ((double) std::abs(n1 - n2) / (double) mag) < epsilon;
}


gsl_rng *g_rng = gsl_rng_alloc(gsl_rng_mt19937);

static MutableBuffer *create_test_mbuffer(size_t cnt)
{
    auto mbuffer = new MutableBuffer(cnt, true, cnt, g_rng);

    for (size_t i = 0; i < cnt; i++) {
        extension::key_t key = rand();
        extension::value_t val = rand();

        mbuffer->append(key, val);
    }

    return mbuffer;
}

static MutableBuffer *create_weighted_mbuffer(size_t cnt)
{
    auto mbuffer = new MutableBuffer(cnt, true, cnt, g_rng);
    
    // Put in half of the count with weight two.
    extension::key_t key = 1;
    for (size_t i=0; i< cnt / 2; i++) {
        mbuffer->append(key, i, 2);
    }

    // put in a quarter of the count with weight four.
    key = 2;
    for (size_t i=0; i< cnt / 4; i++) {
        mbuffer->append(key, i, 4);
    }

    // the remaining quarter with weight eight.
    key = 3;
    for (size_t i=0; i< cnt / 4; i++) {
        mbuffer->append(key, i, 8);
    }

    return mbuffer;
}


static MutableBuffer *create_double_seq_mbuffer(size_t cnt, bool ts=false) 
{
    auto mbuffer = new MutableBuffer(cnt, true, cnt, g_rng);

    for (size_t i = 0; i < cnt / 2; i++) {
        extension::key_t key = i;
        extension::value_t val = i;

        mbuffer->append(key, val, 1.0, ts);
    }

    for (size_t i = 0; i < cnt / 2; i++) {
        extension::key_t key = i;
        extension::value_t val = i + 1;

        mbuffer->append(key, val, 1.0, ts);
    }

    return mbuffer;
}

START_TEST(t_mbuffer_init)
{
    auto mbuffer = new MutableBuffer(1024, true, 1024, g_rng);
    for (uint64_t i = 512; i > 0; i--) {
        uint32_t v = i;
        mbuffer->append(i, v);
    }
    
    for (uint64_t i = 1; i <= 256; ++i) {
        uint32_t v = i;
        mbuffer->append(i, v, 1.0, true);
    }

    for (uint64_t i = 257; i <= 512; ++i) {
        uint32_t v = i + 1;
        mbuffer->append(i, v);
    }

    BloomFilter* bf = new BloomFilter(BF_FPR, mbuffer->get_tombstone_count(), BF_HASH_FUNCS, g_rng);
    ShardWSS* shard = new ShardWSS(mbuffer, bf, false);
    ck_assert_uint_eq(shard->get_record_count(), 512);

    delete bf;
    delete mbuffer;
    delete shard;
}

START_TEST(t_shard_init)
{
    size_t n = 512;
    auto mbuffer1 = create_test_mbuffer(n);
    auto mbuffer2 = create_test_mbuffer(n);
    auto mbuffer3 = create_test_mbuffer(n);

    BloomFilter* bf1 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    BloomFilter* bf2 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    BloomFilter* bf3 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    auto shard1 = new ShardWSS(mbuffer1, bf1, false);
    auto shard2 = new ShardWSS(mbuffer2, bf2, false);
    auto shard3 = new ShardWSS(mbuffer3, bf3, false);

    BloomFilter* bf4 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    ShardWSS* shards[3] = {shard1, shard2, shard3};
    auto shard4 = new ShardWSS(shards, 3, bf4, false);

    ck_assert_int_eq(shard4->get_record_count(), n * 3);
    ck_assert_int_eq(shard4->get_tombstone_count(), 0);

    size_t total_cnt = 0;
    size_t shard1_idx = 0;
    size_t shard2_idx = 0;
    size_t shard3_idx = 0;

    for (size_t i = 0; i < shard4->get_record_count(); ++i) {
        auto rec1 = shard1->get_record_at(shard1_idx);
        auto rec2 = shard2->get_record_at(shard2_idx);
        auto rec3 = shard3->get_record_at(shard3_idx);

        auto cur_rec = shard4->get_record_at(i);

        if (shard1_idx < n && cur_rec->match(rec1)) {
            ++shard1_idx;
        } else if (shard2_idx < n && cur_rec->match(rec2)) {
            ++shard2_idx;
        } else if (shard3_idx < n && cur_rec->match(rec3)) {
            ++shard3_idx;
        } else {
           assert(false);
        }
    }

    delete mbuffer1;
    delete mbuffer2;
    delete mbuffer3;

    delete bf1;
    delete shard1;
    delete bf2;
    delete shard2;
    delete bf3;
    delete shard3;
    delete bf4;
    delete shard4;
}

START_TEST(t_get_lower_bound_index)
{
    size_t n = 10000;
    auto mbuffer = create_double_seq_mbuffer(n);

    ck_assert_ptr_nonnull(mbuffer);
    BloomFilter* bf = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    ShardWSS* shard = new ShardWSS(mbuffer, bf, false);

    ck_assert_int_eq(shard->get_record_count(), n);
    ck_assert_int_eq(shard->get_tombstone_count(), 0);

    auto tbl_records = mbuffer->sorted_output();
    for (size_t i=0; i<n; i++) {
        const record_t *tbl_rec = mbuffer->get_record_at(i);
        auto pos = shard->get_lower_bound(tbl_rec->key);
        ck_assert_int_eq(shard->get_record_at(pos)->key, tbl_rec->key);
        ck_assert_int_le(pos, i);
    }

    delete mbuffer;
    delete bf;
    delete shard;
}


START_TEST(t_full_cancelation)
{
    size_t n = 100;
    auto mbuffer = create_double_seq_mbuffer(n, false);
    auto mbuffer_ts = create_double_seq_mbuffer(n, true);
    BloomFilter* bf1 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    BloomFilter* bf2 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    BloomFilter* bf3 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);

    ShardWSS* shard = new ShardWSS(mbuffer, bf1, false);
    ShardWSS* shard_ts = new ShardWSS(mbuffer_ts, bf2, false);

    ck_assert_int_eq(shard->get_record_count(), n);
    ck_assert_int_eq(shard->get_tombstone_count(), 0);
    ck_assert_int_eq(shard_ts->get_record_count(), n);
    ck_assert_int_eq(shard_ts->get_tombstone_count(), n);

    ShardWSS* shards[] = {shard, shard_ts};

    ShardWSS* merged = new ShardWSS(shards, 2, bf3, false);

    ck_assert_int_eq(merged->get_tombstone_count(), 0);
    ck_assert_int_eq(merged->get_record_count(), 0);

    delete mbuffer;
    delete mbuffer_ts;
    delete bf1;
    delete bf2;
    delete bf3;
    delete shard;
    delete shard_ts;
    delete merged;
}
END_TEST


START_TEST(t_weighted_sampling)
{
    size_t n=1000;
    auto mbuffer = create_weighted_mbuffer(n);

    BloomFilter* bf = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    ShardWSS* shard = new ShardWSS(mbuffer, bf, false);

    extension::key_t lower_key = 0;
    extension::key_t upper_key = 5;

    size_t k = 1000;

    //char * buffer = new char[k*extension::record_size]();
    record_t* buffer = new record_t[k]();
    size_t cnt[3] = {0};
    for (size_t i=0; i<1000; i++) {
        shard->get_samples(buffer, k, nullptr, g_rng);

        for (size_t j=0; j<k; j++) {
            cnt[buffer[j].key - 1]++;
        }
    }

    ck_assert(roughly_equal(cnt[0] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[1] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[2] / 1000, (double) k/2.0, k, .05));

    delete[] buffer;
    delete shard;
    delete bf;
    delete mbuffer;
}
END_TEST


START_TEST(t_tombstone_check)
{
    size_t cnt = 1024;
    size_t ts_cnt = 256;
    auto mbuffer = new MutableBuffer(cnt + ts_cnt, true, ts_cnt, g_rng);

    std::vector<std::pair<extension::key_t, extension::value_t>> tombstones;

    extension::key_t key = 1000;
    extension::value_t val = 101;
    for (size_t i = 0; i < cnt; i++) {
        mbuffer->append(key, val);
        key++;
        val++;
    }

    // ensure that the key range doesn't overlap, so nothing
    // gets cancelled.
    for (size_t i=0; i<ts_cnt; i++) {
        tombstones.push_back({i, i});
    }

    for (size_t i=0; i<ts_cnt; i++) {
        mbuffer->append(tombstones[i].first, tombstones[i].second, 1.0, true);
    }

    BloomFilter* bf1 = new BloomFilter(100, BF_HASH_FUNCS, g_rng);
    auto shard = new ShardWSS(mbuffer, bf1, false);

    for (size_t i=0; i<tombstones.size(); i++) {
        ck_assert(shard->check_tombstone(tombstones[i].first, tombstones[i].second));
        ck_assert_int_eq(shard->get_rejection_count(), i+1);
    }

    delete shard;
    delete mbuffer;
    delete bf1;
}
END_TEST

Suite *unit_testing()
{
    Suite *unit = suite_create("ShardWSS Unit Testing");

    TCase *create = tcase_create("extension::ShardWSS constructor Testing");
    tcase_add_test(create, t_mbuffer_init);
    tcase_add_test(create, t_shard_init);
    tcase_set_timeout(create, 100);
    suite_add_tcase(unit, create);


    TCase *bounds = tcase_create("extension::ShardWSS::get_{lower,upper}_bound Testing");
    tcase_add_test(bounds, t_get_lower_bound_index);
    tcase_set_timeout(bounds, 100);   
    suite_add_tcase(unit, bounds);


    TCase *tombstone = tcase_create("extension::ShardWSS::tombstone cancellation Testing");
    tcase_add_test(tombstone, t_full_cancelation);
    suite_add_tcase(unit, tombstone);


    TCase *sampling = tcase_create("extension::ShardWSS::sampling Testing");
    tcase_add_test(sampling, t_weighted_sampling);

    suite_add_tcase(unit, sampling);

    TCase *check_ts = tcase_create("extension::ShardWSS::check_tombstone Testing");
    tcase_add_test(check_ts, t_tombstone_check);
    suite_add_tcase(unit, check_ts);

    return unit;
}

int shard_unit_tests()
{
    int failed = 0;
    Suite *unit = unit_testing();
    SRunner *unit_shardner = srunner_create(unit);

    srunner_run_all(unit_shardner, CK_NORMAL);
    failed = srunner_ntests_failed(unit_shardner);
    srunner_free(unit_shardner);

    return failed;
}


int main() 
{
    int unit_failed = shard_unit_tests();

    return (unit_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
