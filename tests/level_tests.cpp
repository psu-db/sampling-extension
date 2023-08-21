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

#include "framework/WIRS.h"
#include "framework/Level.h"
#include "base/bf_config.h"

using namespace extension;

gsl_rng *g_rng = gsl_rng_alloc(gsl_rng_mt19937);

static MutableBuffer *create_test_buffer(size_t cnt)
{
    auto buffer = new MutableBuffer(cnt, 0);

    for (size_t i = 0; i < cnt; i++) {
        skey_t key = rand();
        value_t val = rand();

        buffer->append(key, val);
    }

    return buffer;
}


static MutableBuffer *create_double_seq_buffer(size_t cnt) 
{
    auto buffer = new MutableBuffer(cnt, 0);

    for (size_t i = 0; i < cnt / 2; i++) {
        skey_t key = i;
        value_t val = i;

        buffer->append(key, val);
    }

    for (size_t i = 0; i < cnt / 2; i++) {
        skey_t key = i;
        value_t val = i + 1;

        buffer->append(key, val);
    }

    return buffer;
}

START_TEST(t_memlevel_merge)
{
    auto tbl1 = create_test_buffer(100);
    auto tbl2 = create_test_buffer(100);

    auto base_level = new Level(1, 1, false);
    base_level->append_buffer(tbl1);
    ck_assert_int_eq(base_level->get_record_cnt(), 100);

    auto merging_level = new Level(0, 1, false);
    merging_level->append_buffer(tbl2);
    ck_assert_int_eq(merging_level->get_record_cnt(), 100);

    auto old_level = base_level;
    base_level = Level::merge_levels(old_level, merging_level, false);

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
