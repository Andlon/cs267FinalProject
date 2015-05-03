# cs267FinalProject
Parallelizing Gstat for CS267 Final Project

## Serial algorithm pseudo code

1. Read in measurements (value and coordinates)
2. Determine all the pairs of points
3. For each pair of points
    - Calculate squared difference of values
    - Calculate distance from coordinates
    - Add the squared differnce to the appropriate distance bin
4. Divide the sum of squared differences by the number of pairs in bin
5. Print output

## Compiler requirements
This software uses some small subset of recent features from C++14 which your compiler must support. In particular, it should work on gcc version >= 4.9.

## Compilation instructions (Hopper)
From root directory:

    module switch PrgEnv-pgi PrgEnv-gnu
    cd implementations
    make

The output executable will be placed in implementations/bin. For instance, you can run the mpi-naive implementation:

    cd implementations/bin
    aprun -n <numcores> ./mpi-naive <input-path>
    
## Specifying compiler to use
The Makefile is set up to use 'CC' as a compiler. If you want to use a different compiler, you can override the MPICXX (defines the MPI C++ compiler to use) or CXX (defines the C++ compiler to use, for non-MPI modules) variables as arguments to make. For example, if you want to use the *mpic++* compiler:

    make MPICXX=mpic++

