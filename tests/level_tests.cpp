/*
 * tests/level_tests.cpp
 *
 * Copyright (C) 2023 Douglas B. Rumbaugh <drumbaugh@psu.edu> 
 *                    Dong Xie <dongx@psu.edu>
 *
 * All rights reserved. Published under the Simplified BSD License.
 *
 */
#include <check.h>

#include "framework/WSS.h"
#include "framework/Level.h"
#include "util/bf_config.h"

using namespace extension;

gsl_rng *g_rng = gsl_rng_alloc(gsl_rng_mt19937);

static MutableBuffer *create_test_membuffer(size_t cnt)
{
    auto mbuffer = new MutableBuffer(cnt, true, 0, g_rng);

    for (size_t i = 0; i < cnt; i++) {
        extension::key_t key = rand();
        extension::value_t val = rand();

        mbuffer->append(key, val);
    }

    return mbuffer;
}


static MutableBuffer *create_double_seq_membuffer(size_t cnt) 
{
    auto mbuffer = new MutableBuffer(cnt, true, 0, g_rng);

    for (size_t i = 0; i < cnt / 2; i++) {
        extension::key_t key = i;
        extension::value_t val = i;

        mbuffer->append(key, val);
    }

    for (size_t i = 0; i < cnt / 2; i++) {
        extension::key_t key = i;
        extension::value_t val = i + 1;

        mbuffer->append(key, val);
    }

    return mbuffer;
}

START_TEST(t_memlevel_merge)
{
    auto tbl1 = create_test_membuffer(100);
    auto tbl2 = create_test_membuffer(100);

    auto base_level = new Level(1, 1, false);
    base_level->append_buffer(tbl1, g_rng);
    ck_assert_int_eq(base_level->get_record_cnt(), 100);

    auto merging_level = new Level(0, 1, false);
    merging_level->append_buffer(tbl2, g_rng);
    ck_assert_int_eq(merging_level->get_record_cnt(), 100);

    auto old_level = base_level;
    base_level = Level::merge_levels(old_level, merging_level, false, g_rng);

    delete old_level;
    delete merging_level;
    ck_assert_int_eq(base_level->get_record_cnt(), 200);

    delete base_level;
    delete tbl1;
    delete tbl2;
}


Suite *unit_testing()
{
    Suite *unit = suite_create("Level Unit Testing");

    TCase *merge = tcase_create("extension::Level::merge_level Testing");
    tcase_add_test(merge, t_memlevel_merge);
    suite_add_tcase(unit, merge);

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
