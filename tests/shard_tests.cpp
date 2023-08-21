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
#include "base/bf_config.h"

using namespace extension;
using namespace psudb;

bool roughly_equal(int n1, int n2, size_t mag, double epsilon) {
    return ((double) std::abs(n1 - n2) / (double) mag) < epsilon;
}

gsl_rng *g_rng = gsl_rng_alloc(gsl_rng_mt19937);

static MutableBuffer *create_test_buffer(size_t cnt)
{
    auto buffer = new MutableBuffer(cnt, true, cnt);

    for (size_t i = 0; i < cnt; i++) {
        skey_t key = rand();
        value_t val = rand();

        buffer->append(key, val);
    }

    return buffer;
}

static MutableBuffer *create_weighted_buffer(size_t cnt)
{
    auto buffer = new MutableBuffer(cnt, true, cnt);
    
    // Put in half of the count with weight two.
    skey_t key = 1;
    for (size_t i=0; i< cnt / 2; i++) {
        buffer->append(key, i, 2);
    }

    // put in a quarter of the count with weight four.
    key = 2;
    for (size_t i=0; i< cnt / 4; i++) {
        buffer->append(key, i, 4);
    }

    // the remaining quarter with weight eight.
    key = 3;
    for (size_t i=0; i< cnt / 4; i++) {
        buffer->append(key, i, 8);
    }

    return buffer;
}


static MutableBuffer *create_double_seq_buffer(size_t cnt, bool ts=false) 
{
    auto buffer = new MutableBuffer(cnt, true, cnt);

    for (size_t i = 0; i < cnt / 2; i++) {
        skey_t key = i;
        value_t val = i;

        buffer->append(key, val, 1.0, ts);
    }

    for (size_t i = 0; i < cnt / 2; i++) {
        skey_t key = i;
        value_t val = i + 1;

        buffer->append(key, val, 1.0, ts);
    }

    return buffer;
}

START_TEST(t_buffer_init)
{
    auto buffer = new MutableBuffer(1024, true, 1024);
    for (uint64_t i = 512; i > 0; i--) {
        uint32_t v = i;
        buffer->append(i, v);
    }
    
    for (uint64_t i = 1; i <= 256; ++i) {
        uint32_t v = i;
        buffer->append(i, v, 1.0, true);
    }

    for (uint64_t i = 257; i <= 512; ++i) {
        uint32_t v = i + 1;
        buffer->append(i, v);
    }

    BloomFilter<skey_t>* bf = new BloomFilter<skey_t>(BF_FPR, buffer->get_tombstone_count(), BF_HASH_FUNCS);
    ShardWSS* shard = new ShardWSS(buffer, bf, false);
    ck_assert_uint_eq(shard->get_record_count(), 512);

    delete bf;
    delete buffer;
    delete shard;
}

START_TEST(t_shard_init)
{
    size_t n = 512;
    auto buffer1 = create_test_buffer(n);
    auto buffer2 = create_test_buffer(n);
    auto buffer3 = create_test_buffer(n);

    BloomFilter<skey_t>* bf1 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    BloomFilter<skey_t>* bf2 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    BloomFilter<skey_t>* bf3 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    auto shard1 = new ShardWSS(buffer1, bf1, false);
    auto shard2 = new ShardWSS(buffer2, bf2, false);
    auto shard3 = new ShardWSS(buffer3, bf3, false);

    BloomFilter<skey_t>* bf4 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
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

        if (shard1_idx < n && *cur_rec == *rec1) {
            ++shard1_idx;
        } else if (shard2_idx < n && *cur_rec == *rec2) {
            ++shard2_idx;
        } else if (shard3_idx < n && *cur_rec == *rec3) {
            ++shard3_idx;
        } else {
           assert(false);
        }
    }

    delete buffer1;
    delete buffer2;
    delete buffer3;

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
    auto buffer = create_double_seq_buffer(n);

    ck_assert_ptr_nonnull(buffer);
    BloomFilter<skey_t>* bf = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    ShardWSS* shard = new ShardWSS(buffer, bf, false);

    ck_assert_int_eq(shard->get_record_count(), n);
    ck_assert_int_eq(shard->get_tombstone_count(), 0);

    auto tbl_records = buffer->sorted_output();
    for (size_t i=0; i<n; i++) {
        const record_t *tbl_rec = buffer->get_record_at(i);
        auto pos = shard->get_lower_bound(tbl_rec->key);
        ck_assert_int_eq(shard->get_record_at(pos)->key, tbl_rec->key);
        ck_assert_int_le(pos, i);
    }

    delete buffer;
    delete bf;
    delete shard;
}


START_TEST(t_full_cancelation)
{
    size_t n = 100;
    auto buffer = create_double_seq_buffer(n, false);
    auto buffer_ts = create_double_seq_buffer(n, true);
    BloomFilter<skey_t>* bf1 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    BloomFilter<skey_t>* bf2 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    BloomFilter<skey_t>* bf3 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);

    ShardWSS* shard = new ShardWSS(buffer, bf1, false);
    ShardWSS* shard_ts = new ShardWSS(buffer_ts, bf2, false);

    ck_assert_int_eq(shard->get_record_count(), n);
    ck_assert_int_eq(shard->get_tombstone_count(), 0);
    ck_assert_int_eq(shard_ts->get_record_count(), n);
    ck_assert_int_eq(shard_ts->get_tombstone_count(), n);

    ShardWSS* shards[] = {shard, shard_ts};

    ShardWSS* merged = new ShardWSS(shards, 2, bf3, false);

    ck_assert_int_eq(merged->get_tombstone_count(), 0);
    ck_assert_int_eq(merged->get_record_count(), 0);

    delete buffer;
    delete buffer_ts;
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
    auto buffer = create_weighted_buffer(n);

    BloomFilter<skey_t>* bf = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    ShardWSS* shard = new ShardWSS(buffer, bf, false);

    skey_t lower_key = 0;
    skey_t upper_key = 5;

    size_t k = 1000;

    record_t* buf = new record_t[k]();
    size_t cnt[3] = {0};
    for (size_t i=0; i<1000; i++) {
        shard->get_samples(buf, k, nullptr, g_rng);

        for (size_t j=0; j<k; j++) {
            cnt[buf[j].key - 1]++;
        }
    }

    ck_assert(roughly_equal(cnt[0] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[1] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[2] / 1000, (double) k/2.0, k, .05));

    delete[] buf;
    delete shard;
    delete bf;
    delete buffer;
}
END_TEST


START_TEST(t_tombstone_check)
{
    size_t cnt = 1024;
    size_t ts_cnt = 256;
    auto buffer = new MutableBuffer(cnt + ts_cnt, true, ts_cnt);

    std::vector<record_t> tombstones;

    skey_t key = 1000;
    value_t val = 101;
    for (size_t i = 0; i < cnt; i++) {
        buffer->append(key, val);
        key++;
        val++;
    }

    // ensure that the key range doesn't overlap, so nothing
    // gets cancelled.
    for (size_t i=0; i<ts_cnt; i++) {
        tombstones.push_back({i, (value_t) i});
    }

    for (size_t i=0; i<ts_cnt; i++) {
        buffer->append(tombstones[i].key, tombstones[i].value, 1.0, true);
    }

    BloomFilter<skey_t>* bf1 = new BloomFilter<skey_t>(100, BF_HASH_FUNCS);
    auto shard = new ShardWSS(buffer, bf1, false);

    for (size_t i=0; i<tombstones.size(); i++) {
        ck_assert(shard->check_tombstone(tombstones[i]));
        ck_assert_int_eq(shard->get_rejection_count(), i+1);
    }

    delete shard;
    delete buffer;
    delete bf1;
}
END_TEST

Suite *unit_testing()
{
    Suite *unit = suite_create("ShardWSS Unit Testing");

    TCase *create = tcase_create("extension::ShardWSS constructor Testing");
    tcase_add_test(create, t_buffer_init);
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
