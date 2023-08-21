cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(namespace "extension")
project("Practical Dynamic Extension for Sampling Indexes" VERSION 0.1.0)

set(debug false)
set(tests True)
set(bench true)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

if (debug) 
    add_compile_options(-g -O0)
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
    add_compile_options(-fsanitize=undefined)
    add_link_options(-fsanitize=undefined)
else()
    add_compile_options(-O3)
endif()


# Test build instructions
if (tests) 
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/tests")

    add_executable(buffer_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/buffer_tests.cpp)
    target_link_libraries(buffer_tests PUBLIC check subunit pthread gsl)
    target_include_directories(buffer_tests PRIVATE include external/psudb-common/include)

    add_executable(framework_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/framework_tests.cpp)
    target_link_libraries(framework_tests PUBLIC check subunit pthread gsl)
    target_include_directories(framework_tests PRIVATE include external/psudb-common/include)

    add_executable(level_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/level_tests.cpp)
    target_link_libraries(level_tests PUBLIC check subunit pthread gsl)
    target_include_directories(level_tests PRIVATE include external/psudb-common/include)

    add_executable(shard_tests ${CMAKE_CURRENT_SOURCE_DIR}/tests/shard_tests.cpp)
    target_link_libraries(shard_tests PUBLIC check subunit gsl)
    target_include_directories(shard_tests PRIVATE include external/psudb-common/include)
endif()

# Benchmark build instructions
if (bench)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_CURRENT_SOURCE_DIR}/bin/benchmarks")

    add_executable(btree_benchmark ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/btree_benchmark.cpp)
    target_link_libraries(btree_benchmark PUBLIC pthread gsl)
    target_include_directories(btree_benchmark PRIVATE include external/psudb-common/include)

    add_executable(framework_benchmark ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/framework_benchmark.cpp)
    target_link_libraries(framework_benchmark PUBLIC pthread gsl)
    target_include_directories(framework_benchmark PRIVATE include external/psudb-common/include)

    add_executable(static_benchmark ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/static_benchmark.cpp)
    target_link_libraries(static_benchmark PUBLIC pthread gsl)
    target_include_directories(static_benchmark PRIVATE include external/psudb-common/include)
endif()
