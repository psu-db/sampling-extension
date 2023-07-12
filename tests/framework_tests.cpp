/*
 * tests/framework_tests.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#include <check.h>
#include <set>
#include <random>
#include <algorithm>

#include "framework/Framework.h"

using namespace extension;

gsl_rng *g_rng = gsl_rng_alloc(gsl_rng_mt19937);

bool roughly_equal(int n1, int n2, size_t mag, double epsilon) {
    double delta = ((double) std::abs(n1 - n2) / (double) mag);
    return delta < epsilon;
}

START_TEST(t_create)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);

    ck_assert_ptr_nonnull(dynamic_extension);
    ck_assert_int_eq(dynamic_extension->get_record_cnt(), 0);
    ck_assert_int_eq(dynamic_extension->get_height(), 0);

    delete dynamic_extension;
}
END_TEST


START_TEST(t_append)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);

    extension::key_t key = 0;
    extension::value_t val = 0;
    for (size_t i=0; i<100; i++) {
        ck_assert_int_eq(dynamic_extension->append(key, val, 1, false, g_rng), 1);
        key++;
        val++;
    }

    ck_assert_int_eq(dynamic_extension->get_height(), 0);
    ck_assert_int_eq(dynamic_extension->get_record_cnt(), 100);

    delete dynamic_extension;
}
END_TEST


START_TEST(t_append_with_mem_merges)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);

    extension::key_t key = 0;
    extension::value_t val = 0;
    for (size_t i=0; i<300; i++) {
        ck_assert_int_eq(dynamic_extension->append(key, val, 1, false, g_rng), 1);
        key++;
        val++;
    }

    ck_assert_int_eq(dynamic_extension->get_record_cnt(), 300);
    ck_assert_int_eq(dynamic_extension->get_height(), 1);

    delete dynamic_extension;
}
END_TEST


START_TEST(t_range_sample_mbuffer)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);

    extension::key_t key = 0;
    extension::value_t val = 0;
    for (size_t i=0; i<100; i++) {
        ck_assert_int_eq(dynamic_extension->append(key, val, 1, false, g_rng), 1);
        key++;
        val++;
    }

    extension::key_t lower_bound = 0;
    extension::key_t upper_bound = 100;

    //char sample_set[100*record_size];
    record_t sample_set[100];

    dynamic_extension->range_sample(sample_set, 100, g_rng);

    for(size_t i=0; i<100; i++) {
        ck_assert_int_le(sample_set[i].key, upper_bound);
        ck_assert_int_ge(sample_set[i].key, lower_bound);
    }

    delete dynamic_extension;
}
END_TEST


START_TEST(t_range_sample_levels)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);

    extension::key_t key = 0;
    extension::value_t val = 0;
    for (size_t i=0; i<300; i++) {
        ck_assert_int_eq(dynamic_extension->append(key, val, 1, false, g_rng), 1);
        key++;
        val++;
    }

    extension::key_t lower_bound = 0;
    extension::key_t upper_bound = 300;

    record_t sample_set[100];
    dynamic_extension->range_sample(sample_set, 100, g_rng);

    for(size_t i=0; i<100; i++) {
        ck_assert_int_le(sample_set[i].key, upper_bound);
        ck_assert_int_ge(sample_set[i].key, lower_bound);
    }

    delete dynamic_extension;
}
END_TEST

START_TEST(t_range_sample_weighted)
{
    auto dynamic_extension = new SamplingFramework(100, 2, 1, 100, g_rng);
    size_t n = 10000;

    std::vector<extension::key_t> keys;

    extension::key_t key = 1;
    for (size_t i=0; i< n / 2; i++) {
        keys.push_back(key);
    }

    // put in a quarter of the count with weight two.
    key = 2;
    for (size_t i=0; i< n / 4; i++) {
        keys.push_back(key);
    }

    // the remaining quarter with weight four.
    key = 3;
    for (size_t i=0; i< n / 4; i++) {
        keys.push_back(key);
    }

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::shuffle(keys.begin(), keys.end(), gen);

    for (size_t i=0; i<keys.size(); i++) {
        double weight;
        if (keys[i] == 1)  {
            weight = 2.0;
        } else if (keys[i] == 2) {
            weight = 4.0;
        } else {
            weight = 8.0;
        }

        dynamic_extension->append(keys[i], i, weight, false, g_rng);
    }
    size_t k = 1000;
    extension::key_t lower_key = 0;
    extension::key_t upper_key = 5;

    //char *buffer = new char[k*extension::record_size]();
    record_t* buffer = new record_t[k]();

    size_t cnt[3] = {0};
    for (size_t i=0; i<1000; i++) {
        dynamic_extension->range_sample(buffer, k, g_rng);

        for (size_t j=0; j<k; j++) {
            cnt[buffer[j].key - 1]++;
        }
    }

    ck_assert(roughly_equal(cnt[0] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[1] / 1000, (double) k/4.0, k, .05));
    ck_assert(roughly_equal(cnt[2] / 1000, (double) k/2.0, k, .05));

    delete dynamic_extension;
    delete[] buffer;
}
END_TEST


START_TEST(t_tombstone_merging_01)
{
    size_t reccnt = 100000;
    auto dynamic_extension = new SamplingFramework(100, 2, .01, 100, g_rng);

    std::set<std::pair<extension::key_t, extension::value_t>> records; 
    std::set<std::pair<extension::key_t, extension::value_t>> to_delete;
    std::set<std::pair<extension::key_t, extension::value_t>> deleted;

    while (records.size() < reccnt) {
        extension::key_t key = rand();
        extension::value_t val = rand();

        if (records.find({key, val}) != records.end()) continue;

        records.insert({key, val});
    }

    size_t deletes = 0;
    size_t cnt=0;
    for (auto rec : records) {
        //const char *key_ptr = (char *) &rec.first;
        //const char *val_ptr = (char *) &rec.second;
        ck_assert_int_eq(dynamic_extension->append(rec.first, rec.second, 1, false, g_rng), 1);

         if (gsl_rng_uniform(g_rng) < 0.05 && !to_delete.empty()) {
            std::vector<std::pair<extension::key_t, extension::value_t>> del_vec;
            std::sample(to_delete.begin(), to_delete.end(), std::back_inserter(del_vec), 3, std::mt19937{std::random_device{}()});

            for (size_t i=0; i<del_vec.size(); i++) {
                //const char *d_key_ptr = (char *) &del_vec[i].first;
                //const char *d_val_ptr = (char *) &del_vec[i].second;
                if (extension::DELETE_TAGGING) {
                    dynamic_extension->delete_record(del_vec[i].first, del_vec[i].second, g_rng);
                } else {
                    dynamic_extension->append(del_vec[i].first, del_vec[i].second, 1, true, g_rng);
                }
                deletes++;
                to_delete.erase(del_vec[i]);
                deleted.insert(del_vec[i]);
            }
        }

        if (gsl_rng_uniform(g_rng) < 0.25 && deleted.find(rec) == deleted.end()) {
            to_delete.insert(rec);
        }

        ck_assert(dynamic_extension->validate_tombstone_proportion());
    }

    ck_assert(dynamic_extension->validate_tombstone_proportion());

    delete dynamic_extension;
}
END_TEST


Suite *unit_testing()
{
    Suite *unit = suite_create("extension::SamplingFramework Unit Testing");

    TCase *create = tcase_create("extension::SamplingFramework::constructor Testing");
    tcase_add_test(create, t_create);
    suite_add_tcase(unit, create);

    TCase *append = tcase_create("extension::SamplingFramework::append Testing");
    tcase_add_test(append, t_append);
    tcase_add_test(append, t_append_with_mem_merges);
    suite_add_tcase(unit, append);

    TCase *sampling = tcase_create("extension::SamplingFramework::range_sample Testing");
    tcase_add_test(sampling, t_range_sample_mbuffer);
    tcase_add_test(sampling, t_range_sample_levels);
    tcase_add_test(sampling, t_range_sample_weighted);

    suite_add_tcase(unit, sampling);


    TCase *ts = tcase_create("extension::SamplingFramework::tombstone_compaction Testing");
    tcase_add_test(ts, t_tombstone_merging_01);
    tcase_set_timeout(ts, 500);
    suite_add_tcase(unit, ts);

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
