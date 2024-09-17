# Practical Dynamic Extension for Sampling Indexes

This repository contains the source code associated with the SIGMOD 2024 
conference paper: Practical Dynamic Extension for Sampling Indexes, by
Douglas B. Rumbaugh and Dong Xie. The contents of this repository are a
result of a clean-up/restructuring effort to make the code easier to work
with. The original code can be found in the 
[sampling-extension-original](https://github.com/psu-db/sampling-extension-original)
repository.

Note that this clean-up only includes code for WSS and WIRS use-cases. For IRS,
this specialized implementation has been superceeded by the implementation discussed
in the VLDB 2024 paper: Towards Systematic Index Dynamization, but
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
The code is provided in the form of a header-only library, with all relevant
files located in the `include` directory. 


