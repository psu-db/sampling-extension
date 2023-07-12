# Practical Dynamic Extension for Sampling Indexes

This repository contains the source code associated with the SIGMOD 2024 
conference paper: Practical Dynamic Extension for Sampling Indexes, by
Douglas B. Rumbaugh and Dong Xie. 

## Repository Organization
The project itself is split into four branches with specific implementations:

1. wss -- code associated with the DE-WSS extended structure, based on Walker's
          Alias Structure for weighted set sampling.
2. wirs -- code associated with the DE-WIRS extended structure, based on an
           Alias-augmented B+Tree for weighted independent range sampling.
3. irs -- code associated with the DE-IRS extended structure, based on a
          static ISAM tree, with both in-memory and on-disk implementations, 
          for independent range sampling.
4. concurrent -- code associated with a concurrent version of the DE-IRS
                 extended structure, for independent range sampling.

## Project Organization
The code is provided in the form of a header-only library, with all relevant
files located in the `include` directory. 


