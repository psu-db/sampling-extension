/*
 * tests/buffer_tests.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Revised BSD License.
 *
 */

#include <check.h>
#include <string>
#include <thread>
#include <gsl/gsl_rng.h>
#include <vector>
#include <algorithm>

#include "framework/MutableBuffer.h"

using namespace extension;
using namespace psudb;

START_TEST(t_create)
{
    auto buffer = new MutableBuffer(100, 50);

    ck_assert_ptr_nonnull(buffer);
    ck_assert_int_eq(buffer->get_capacity(), 100);
    ck_assert_int_eq(buffer->get_record_count(), 0);
    ck_assert_int_eq(buffer->is_full(), false);
    ck_assert_ptr_nonnull(buffer->sorted_output());
    ck_assert_int_eq(buffer->get_tombstone_count(), 0);
    ck_assert_int_eq(buffer->get_tombstone_capacity(), 50);

    delete buffer;
}
END_TEST


START_TEST(t_insert)
{
    auto buffer = new MutableBuffer(100, 50);

    skey_t key = 0;
    value_t val = 5;

    for (size_t i=0; i<99; i++) {
        ck_assert_int_eq(buffer->append(key, val, 1.0, false), 1);
        ck_assert_int_eq(buffer->check_tombstone({key, val}), 0);

        key++;
        val++;

        ck_assert_int_eq(buffer->get_record_count(), i+1);
        ck_assert_int_eq(buffer->get_tombstone_count(), 0);
        ck_assert_int_eq(buffer->is_full(), 0);
    }

    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 1);

    key++;
    val++;

    ck_assert_int_eq(buffer->is_full(), 1);
    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 0);

    delete buffer;

}
END_TEST


START_TEST(t_insert_tombstones)
{
    auto buffer = new MutableBuffer(100, 50);

    skey_t key = 0;
    value_t val = 5;
    size_t ts_cnt = 0;

    for (size_t i=0; i<99; i++) {
        bool ts = false;
        if (i % 2 == 0) {
            ts_cnt++;
            ts=true;
        }

        ck_assert_int_eq(buffer->append(key, val, 1.0, ts), 1);
        ck_assert_int_eq(buffer->check_tombstone({key, val}), ts);

        key++;
        val++;

        ck_assert_int_eq(buffer->get_record_count(), i+1);
        ck_assert_int_eq(buffer->get_tombstone_count(), ts_cnt);
        ck_assert_int_eq(buffer->is_full(), 0);
    }

    // inserting one more tombstone should not be possible
    ck_assert_int_eq(buffer->append(key, val, 1.0, true), 0);


    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 1);

    key++;
    val++;

    ck_assert_int_eq(buffer->is_full(), 1);
    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 0);

    delete buffer;
}
END_TEST


START_TEST(t_truncate)
{
    auto buffer = new MutableBuffer(100, 100);

    skey_t key = 0;
    value_t val = 5;
    size_t ts_cnt = 0;

    for (size_t i=0; i<100; i++) {
        bool ts = false;
        if (i % 2 == 0) {
            ts_cnt++;
            ts=true;
        }

        ck_assert_int_eq(buffer->append(key, val, 1.0, ts), 1);
        ck_assert_int_eq(buffer->check_tombstone({key, val}), ts);

        key++;
        val++;

        ck_assert_int_eq(buffer->get_record_count(), i+1);
        ck_assert_int_eq(buffer->get_tombstone_count(), ts_cnt);
    }

    ck_assert_int_eq(buffer->is_full(), 1);
    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 0);

    ck_assert_int_eq(buffer->truncate(), 1);

    ck_assert_int_eq(buffer->is_full(), 0);
    ck_assert_int_eq(buffer->get_record_count(), 0);
    ck_assert_int_eq(buffer->get_tombstone_count(), 0);
    ck_assert_int_eq(buffer->append(key, val, 1.0, false), 1);

    delete buffer;
}
END_TEST


START_TEST(t_sorted_output)
{
    size_t cnt = 100;

    auto buffer = new MutableBuffer(cnt, cnt/2);

    std::vector<skey_t> keys(cnt);
    for (size_t i=0; i<cnt-2; i++) {
        keys[i] = rand();
    }

    // duplicate final two records for tombstone testing
    // purposes
    keys[cnt-2] =  keys[cnt-3];
    keys[cnt-1] =  keys[cnt-2];

    value_t val = 12345;
    for (size_t i=0; i<cnt-2; i++) {
        buffer->append(keys[i], val, 1.0, false);
    }

    buffer->append(keys[cnt-2], val, 1.0, true);
    buffer->append(keys[cnt-1], val, 1.0, true);


    record_t *sorted_records = buffer->sorted_output();
    std::sort(keys.begin(), keys.end());

    for (size_t i=0; i<cnt; i++) {
        ck_assert_int_eq(sorted_records[i].key, keys[i]);
    }

    delete buffer;
}
END_TEST


void insert_records(std::vector<std::pair<skey_t, value_t>> *values, size_t start, size_t stop, MutableBuffer *buffer)
{
    for (size_t i=start; i<stop; i++) {
        buffer->append((*values)[i].first, (*values)[i].second, 1.0);
    }

}

START_TEST(t_multithreaded_insert)
{
    size_t cnt = 10000;
    auto buffer = new MutableBuffer(cnt, cnt/2);

    std::vector<std::pair<skey_t, value_t>> records(cnt);
    for (size_t i=0; i<cnt; i++) {
        records[i] = {rand(), rand()};
    }

    // perform a t_multithreaded insertion
    size_t thread_cnt = 8;
    size_t per_thread = cnt / thread_cnt;
    std::vector<std::thread> workers(thread_cnt);
    size_t start = 0;
    size_t stop = start + per_thread;
    for (size_t i=0; i<thread_cnt; i++) {
        workers[i] = std::thread(insert_records, &records, start, stop, buffer);
        start = stop;
        stop = std::min(start + per_thread, cnt);
    }

    for (size_t i=0; i<thread_cnt; i++) {
        if (workers[i].joinable()) {
            workers[i].join();
        }
    }

    ck_assert_int_eq(buffer->is_full(), 1);
    ck_assert_int_eq(buffer->get_record_count(), cnt);

    std::sort(records.begin(), records.end());
    record_t *sorted_records = buffer->sorted_output();
    for (size_t i=0; i<cnt; i++) {
        ck_assert_int_eq(sorted_records[i].key, records[i].first);
    }

    delete buffer;
}
END_TEST


Suite *unit_testing()
{
    Suite *unit = suite_create("MutableBuffer Unit Testing");
    TCase *initialize = tcase_create("framework::MutableBuffer Constructor Testing");
    tcase_add_test(initialize, t_create);

    suite_add_tcase(unit, initialize);


    TCase *append = tcase_create("extension::MutableBuffer::append Testing");
    tcase_add_test(append, t_insert);
    tcase_add_test(append, t_insert_tombstones);
    tcase_add_test(append, t_multithreaded_insert);

    suite_add_tcase(unit, append);


    TCase *truncate = tcase_create("extension::MutableBuffer::truncate Testing");
    tcase_add_test(truncate, t_truncate);

    suite_add_tcase(unit, truncate);


    TCase *sorted_out = tcase_create("extension::MutableBuffer::sorted_output");
    tcase_add_test(sorted_out, t_sorted_output);

    suite_add_tcase(unit, sorted_out);

    return unit;
}


int run_unit_tests()
{
    int failed = 0;
    Suite *unit = unit_testing();
    SRunner *unit_runner = srunner_create(unit);

    srunner_run_all(unit_runner, CK_NORMAL);
    failed = srunner_ntests_failed(unit_runner);
    srunner_free(unit_runner);

    return failed;
}


int main() 
{
    int unit_failed = run_unit_tests();

    return (unit_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

