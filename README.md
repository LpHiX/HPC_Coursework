# HPC_Coursework

This repository contains the code I've written for the HPC coursework assignment.

# Notes to self:
Add `-fopt-info-vec-missed` flag when I want to debug vectorization.

# TODO:
1. make sure both --forcing and --test don't work together
    1. Addendum, this doesn't work, ask chris cantwell as test has a default value
2. Pull initializing out of the solver
    1. This means both u, f
3. Refactor testing
4. Have a consisten naming scheme
6. Swap i j k, x, y, z
7. file read and write MPI
8. write MPI tests
9. rewrite serial to use MPI code
10. add initializer to mpi