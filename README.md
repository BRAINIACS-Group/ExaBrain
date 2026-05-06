This branch contains the code used for benchmarking within the dealii-X project. 
Severel parameters files with various system sizes are available and ready to be 
run with the "nonlinear-poro-viscoelasticity.cc" code.

To compile the code run:

    cmake .
    make all
    make release
    make

Run the examples using MPI:

    mpirun -np <number of processors> <path to nonlinear-poro-viscoleasticity executable> <number of threads> <path to parameter file>
