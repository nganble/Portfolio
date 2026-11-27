------------------------------------------------------------
COMPILATION
------------------------------------------------------------
A Makefile is included to manage the build process.

To compile all versions:
    make

To compile specific versions:
    make serial
    make mpi
    make cuda
    make mpi_cuda

------------------------------------------------------------
EXECUTION
------------------------------------------------------------
All executables accept the following command-line arguments:
./executable [N] [omega] [tolerance] [max_iterations]

Where:
  N              : Grid points per dimension (Total nodes = N^2)
  omega          : Relaxation parameter for SOR between 0.0 and 2.0, use 2.0 for optimal 
  tolerance      : Stopping criterion for the maximum residual
  max_iterations : Maximum iteration limit

1. Serial:
    ./poisson2d_serial 1024 2.0 1e-10 20000

2. MPI (Example using 2 ranks):
    mpirun -np 2 ./poisson2d_mpi 1024 2.0 1e-10 20000

3. CUDA:
    ./poisson2d_cuda 1024 2.0 1e-10 20000

4. Hybrid MPI + CUDA (Example using 2 ranks):
    mpirun -np 2 ./poisson2d_mpi_cuda 1024 2.0 1e-10 20000