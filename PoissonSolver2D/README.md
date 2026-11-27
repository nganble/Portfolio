# Parallel Solvers for the Poisson Equation in Two Dimensions

## Overview
This project solves the Poisson equation $u_{xx}+u_{yy}=f$ on the unit square $(0,1)^2$ with Dirichlet boundary conditions. The work includes:
- A prototype in MATLAB
- A serial implementation in C
- Parallel solvers using MPI, CUDA, and MPI+CUDA

## Features
- **Numerical method**:
  - Second-order 5-point finite difference approximation of the Laplacian.
  - Red-black SOR($\omega$): points of the same color are independent, so each half-step can be updated in parallel.
  - Input $\omega = 2.0$ selects the optimal value $\omega = 2/(1+\sin(\pi h))$. The iteration stops when the maximum residual is below the tolerance.

- **Implementation**:
  - The C codes share `poisson2d_common.h` (types, grid macros, test cases) and `timing.h` (timer).
  - MPI: the grid is split in the $y$-direction; each rank exchanges one ghost line with its neighbours (`MPI_Isend`/`MPI_Irecv`) after the red and after the black update.
  - CUDA: kernels for the red/black update, the residual, and the max reduction, with unified memory.
  - MPI+CUDA: same decomposition as MPI, one GPU per rank.
  - The test case is chosen at compile time with `-DTEST=EXACT_SOLUTION -DSOLUTION=POLY|GAUSSIAN` or `-DTEST=COMPUTE_SOLUTION` (default: a source at $(0.5,0.5)$ and a sink at $(0.25,0.25)$).

- **Verification**:
  - Polynomial (cubic in $x$ and $y$) exact solution: the scheme is exact (errors at the level of the tolerance).
  - Gaussian exact solution: second-order convergence (error ratio 3.9).
  - The serial and MPI codes give the same iterations, residuals, and errors for every number of processes. `src/runBaseline.sh` checks all cases against saved references.

## Results
Runs on the AiMOS cluster at RPI (IBM Power9 nodes with 32 cores and 4 NVIDIA V100 GPUs), tolerance $10^{-10}$:

**MPI, strong scaling** (runtime in seconds, speedup relative to $n_p=2$):

| $N$  | $n_p=2$ | $n_p=128$ | Speedup (ideal 64) |
|-----:|--------:|----------:|-------------------:|
| 512  | 4.18    | 0.156     | 26.9               |
| 1024 | 33.10   | 1.198     | 27.6               |

**MPI, weak scaling** ($N^2/n_p$ constant): the time per iteration grows from 0.61 ms ($N=256$, $n_p=2$) to 0.99 ms ($N=2048$, $n_p=128$), about 61%, because of communication.

**MPI+CUDA**: for $N \le 2048$ a single GPU is fastest; at $N=4096$ the runtime drops from 23.3 s (1 GPU) to 17.8 s (16 GPUs).

These timings were measured with an earlier version of the code, whose forcing was centered at $(\pm 1, 0)$. Full results and discussion are in [docs/poissonReport.pdf](docs/poissonReport.pdf).

## Repository layout
```text
src/         poisson2d_serial.c, poisson2d_mpi.c, poisson2d_cuda.cu, poisson2d_mpi_cuda.{c,cu},
             shared headers, Makefile, runBaseline.sh, Slurm scripts run_mpi.sh, run_mpi_cuda.sh
matlab/      main.m (SOR prototype), getSolution.m, plot scripts for the report figures
reference/   saved baseline results (serial code)
docs/        poissonReport.tex/.pdf (SOR theory, implementation, results), figures
```

## To run
**Requirements:** a C compiler, an MPI library (e.g. Open MPI); `nvcc` for the CUDA codes; MATLAB for the prototype.

**Build** (in `src/`):
```bash
make                          # serial and MPI
make cuda mpi-cuda            # CUDA codes (needs nvcc)
make POWER9=1 all             # IBM Power9 cluster (mpixlc, Power9 timer)
make clean; make CFLAGS="-DTEST=EXACT_SOLUTION -DSOLUTION=GAUSSIAN"   # other test case
```

**Run:** all codes take `N omega tol max_iter` ($N$ interior points per direction, $h = 1/(N+1)$):
```bash
./poisson2d_serial 1024 2.0 1e-10 20000
mpirun -np 4 ./poisson2d_mpi 1024 2.0 1e-10 20000
./poisson2d_cuda 1024 2.0 1e-10 20000
mpirun -np 4 ./poisson2d_mpi_cuda 1024 2.0 1e-10 20000
```

**Baseline** (serial and MPI, 3 cases, $N = 16, 32, 64, 128$; prints PASS/FAIL):
```bash
src/runBaseline.sh            # "save" writes new references
```

**Cluster jobs** (Slurm, from `src/`):
```bash
mkdir -p results
sbatch run_mpi.sh             # MPI strong and weak scaling
sbatch run_mpi_cuda.sh        # MPI+CUDA and CUDA strong scaling
```

## Documentation
[docs/poissonReport.pdf](docs/poissonReport.pdf): Jacobi, Gauss-Seidel and SOR methods, their convergence for the model problem, the parallel implementation, and the numerical results.
