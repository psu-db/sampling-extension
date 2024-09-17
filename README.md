# Practical Dynamic Extension for Sampling Indexes

This repository contains the source code associated with the SIGMOD 2024 
conference paper: Practical Dynamic Extension for Sampling Indexes, by
Douglas B. Rumbaugh and Dong Xie. The contents of this repository are a
result of a clean-up/restructuring effort to make the code easier to work
with. The original code used for the paper can be found in the 
[sampling-extension-original](https://github.com/psu-db/sampling-extension-original)
repository. Be warned--it's a mess.

Note that this clean-up only includes code for WSS and WIRS use-cases. For IRS,
this specialized implementation has been superceeded by the implementation discussed
in the VLDB 2024 paper: Towards Systematic Index Dynamization, by
Douglas B. Rumbaugh, Dong Xie, and Zhuoyue Zhao. See 
[the repository for this paper](https://github.com/psu-db/dynamic-extension)
for a cleaner/more modern implementation of dynamized IRS based on ISAM trees,
or the aforementioned original repository for the old-style implementation.

## Repository Organization
The project itself is split into two branches with specific implementations:

1. wss  -- code associated with the DE-WSS extended structure, based on Walker's
           Alias Structure for weighted set sampling.
2. wirs -- code associated with the DE-WIRS extended structure, based on an
           Alias-augmented B+Tree for weighted independent range sampling.

## Project Organization
The code is provided in the form of a header-only library, with all
relevant files located in the `include` directory. The `include/framework`
directory includes the implementation of the dynamization framework
itself, with `Framework.h` containing the user-facing framework object
and interface.  The files `Level.h` and `MutableBuffer.h` contain
the internal implementations of the multiple levels, and the buffer,
respectively. The `WSS.h` and `WIRS.h`, depending on the branch, contain
the shard implementation of the static data structure to be extended.

## Building the Project

### Dependencies
This project depends upon two external libraries that are not included
in the repository: `check` for unit testing, and `gsl` for random
number generation. These can be installed on Ubuntu using,
```
$ sudo apt install check libgsl-dev
```

Additionally, the project uses resources from the psudb-common library,
which has been included as a git submodule that must be initialized with,
```
$ git submodule update --init        
```
before the project can be build.

### Build Instructions
The project can be built using cmake once the dependencies have been 
installed/intialized,
```
$ mkdir build
$ cd build
$ cmake ..
$ make -j           
```
The `cmakelists.txt` file includes three significant flags that can
be set.  The `debug` flag will build with `-O3` if false, or `-g -O0`
with ASAN and UBSAN if true. The `tests` and `bench` flag control
the building of unit tests and benchmarks respectively.

Unit tests will be placed in `bin/tests/` and benchmarks in
`bin/benchmarks`.

## Framework Configuration and Usage
The framework itself is contained in `Framework.h` as `SamplingFramework`. 
The implementation is branch-specific, with the `wss` branch containing
the implementation of the framework for weighted set sampling, and `wirs`
containing weighted independent range sampling. 

### Configuration Parameters

The framework is configurable using a mix of compile-time and run-time
parameters. The compile time parameters are contained in `Framework.h`
and are,
| name | purpose | values |
|------|---------|--------|
|REJECTION_SAMPLING | configures whether sampling against the buffer is done using rejection, or with an initial full scan | true, false |
|LAYOUT_POLICY | toggle between the leveling and tiering layout policy | tue (leveling), false (tiering) |
|DELETE_POLICY | toggle between tombstone and tagging deletes | true (tagging), false (tombstone) | 

The runtime configuration is handled using constructor parameters. These
options are,
| name | purpose | values |
|------|---------|--------|
|buffer_cap| determines the record capacity of the buffer | any non-negative integer |
|scale_factor| adjusts the rate at which levels are added to the structure | any non-negative integer |
| max_tombstone_prop | the maximum proportion of deleted records on a level before a proactive-reconstruction is triggered | \[0,1\], with 1 disabling the behavior |
| max_rejection_rate | the maximum proportion of sampling rejections on a level before a proactive-reconstruciton is triggered | \[0,1\] with 1 disabling the behavior |

More information on all of these configuration parameters can be found in
the paper itself, including benchmark results showing their performance
effects.

### Methods

The framework allows for the insertion and deletion of records,
as well as sampling queries, to be executed. The datatypes are
hard-coded into the library (see the [more modern dynamization
framework](https://github.com/psu-db/dynamic-extension) for a
version that uses C++ templates instead). The record format is
defined in `include/base/record.h` and the key/value/weight types in
`include/base/types.h` and can be adjusted here as needed, though testing
was only done using the format currently defined, which uses a uint64_t
key, uint32_t value, uint32_t header, and uint64_t weight.

The methods supported by the framework are as follows,

1. `append(key, value, weight, tombstone)`
    Insert a new weighted key-value pair into the structure. `tombstone`
    should be set to `false` to insert a record. If tombstone-based deletes
    are in use, a record can be deleted by inserting it using this method,
    with `tombstone` set to `true`. When deleting a record via tombstones,
    the key and value must match exactly the record to be deleted, but the
    weight value is ignored.

    Returns 1 if the insert succeeds, and 0 if it does not.

2. `delete_record(key, value)`
    When TAGGING is used as the delete policy, this method performs a lookup
    for a record matching the specified key and value and marks it as deleted
    in its header when found. 

    Returns 1 if the delete succeeds and 0 if it does not. The only way for
    the delete to fail is if a record matching the specified key and value
    is not found.

3. `range_sample(sample_set, sample_sz, rng)`
   Draws a weighted sample of `sample_sz` records iid from the inserted
   records and places them into the buffer specified by `sample_set`. If
   `rng` is null, or the provided buffer is not large enough, the behavior
   of this method is undefined.

4. `create_static_structure()`
   Creates and returns a single, combined instance of the static data
   structure being dynamized containing all of the records currently
   within the framework.

The framework also supports some accessor methods for querying information
about the structure, including record count (including deleted records
and tombstones), tombstone count, height (i.e., the number of levels),
auxilliary memory usage (memory used for structures other than the static
one being extended), and memory usage (memory used for data structures
being extended)
