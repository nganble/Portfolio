# Parallel Solvers for the Anisotropic Wave Equation

## Overview
This project solves the anisotropic wave equation $\rho w_{tt}=\nabla\cdot(m\nabla w) + f$ on a rectangular domain $(a_x,b_x)\times(a_y,b_y)$ with Dirichlet boundary conditions. The work includes:
- A prototype in MATLAB
- A serial implementation in C++
- Parallel solvers using OpenMP and MPI for speedup

## Features
- **Numerical method**:
  - Explicit scheme $\rho\, D_{+t}D_{-t}W = L_h W + f$, second-order accurate in space and time.
  - $L_h$ averages $m_{11}$, $m_{22}$ at the half points and uses pointwise centered differences for the cross terms $m_{12}$, $m_{21}$. SeeSee the [note](docs/anisotropicWaveScheme.pdf).
  - Time step: freeze the coefficients, do a von Neumann analysis, and take the minimum over the grid; $\Delta t = 0.9\,\Delta t_{\max}$.

- **Implementation**:
  - The C++ codes share the headers in `common/` (test cases, grid and array layout, stencil, output).
  - OpenMP: the update loop is parallelized over grid lines; loop constants are `firstprivate`.
  - MPI: the grid is split in the $y$-direction; each rank exchanges one ghost line with its neighbours using `MPI_Sendrecv`.
  - The arrays use a padded row stride to avoid cache conflicts.
  - The test case is chosen at compile time with `-DCOEFF=CONS|VAR` and `-DSOLUTION=POLY|TRIG|PULSE|PULSETRIG`.

- **Verification**:
  - Polynomial (degree 2) manufactured solution with constant coefficients: the scheme is exact (errors at machine precision).
  - Trigonometric and pulse (plane wave) exact solutions with constant coefficients: second-order convergence.
  - Polynomial manufactured solution with variable coefficients: second-order convergence.
  - The serial, OpenMP and MPI codes give the same results and agree with the MATLAB prototype to round-off. Baseline scripts check all 5 cases against saved references.

## Results
Speedup $S = T_{\text{serial}}/T_p$ on a workstation (AMD EPYC 9455P, 48 cores), $N_x = N_y = 3000$:

| Cores | OpenMP, const | OpenMP, var | MPI, const | MPI, var |
|------:|------:|------:|------:|------:|
| 8     | 5.5   | 7.9   | 7.2   | 7.8   |
| 16    | 11.2  | 14.8  | 13.4  | 14.7  |
| 32    | 21.6  | 26.9  | 20.8  | 24.8  |

Full tables and performance notes are in [docs/anisotropicWaveScheme.pdf](docs/anisotropicWaveScheme.pdf).

## Repository layout
```text
common/            shared C++ headers: test cases, grid/array layout, stencil, output, timer, options
serial/            wave2dSerial.C, Makefile, runBaselineSerial.sh
parallel_openmp/   wave2dOpenMP.C, Makefile, runBaselineOpenMP.sh, runTimingOpenMP.sh
parallel_mpi/      wave2dMPI.C, Makefile, runBaselineMPI.sh, runTimingMPI.sh
matlab/            anisotropicWave_2D.m (prototype), runBaselineMATLAB.m
reference/         saved baseline results (MATLAB, serial) and timing tables
docs/              anisotropicWaveScheme.pdf (scheme, test cases, results, speedup)
```

## To run
**Requirements:** a C++ compiler with OpenMP, an MPI library (e.g. Open MPI), MATLAB (tested with R2024a). On macOS, Apple clang has no OpenMP; the OpenMP Makefile uses Homebrew GCC (`g++-16`).

**MATLAB:** in `matlab/`, run `anisotropicWave_2D.m` (set `ms`, `coeffs`, `makeMovie` at the top) or `runBaselineMATLAB.m`.

**C++:** build in each code folder with `make` (default case: constant coefficients, polynomial solution). Other cases:
```bash
make clean; make CXXFLAGS="-O3 -Wall -DCOEFF=VAR -DSOLUTION=POLY"   # add -fopenmp for OpenMP
```
Run:
```bash
./wave2dSerial -nx=80 -tFinal=0.5
./wave2dOpenMP -nx=80 -tFinal=0.5 -numThreads=4
mpirun -np 4 ./wave2dMPI -nx=80 -tFinal=0.5
```
Options: `-nx`, `-ny`, `-tFinal`, `-debug`; `-numThreads` (OpenMP); `-saveMatlab`, `-matlabFileName` (MPI).

**Baselines** (5 cases, $N_x = 20, 40, 80, 160$; compare with `reference/`, print PASS/FAIL):
```bash
serial/runBaselineSerial.sh            # "save" writes new references
parallel_openmp/runBaselineOpenMP.sh
parallel_mpi/runBaselineMPI.sh
```

**Timing** (results saved in `reference/`):
```bash
MACHINE=workstation THREADS="1 2 4 8 16 32" parallel_openmp/runTimingOpenMP.sh
MACHINE=workstation NP="1 2 4 8 16 32" parallel_mpi/runTimingMPI.sh
COEFF=VAR TFINAL=0.01 parallel_openmp/runTimingOpenMP.sh       # variable coefficients
```
Other options: `NX` (default 3000; avoid sizes near a power of two), `TFINAL`, `REPEAT`.

## Documentation
[docs/anisotropicWaveScheme.pdf](docs/anisotropicWaveScheme.pdf): numerical scheme, test cases, convergence results, and parallel speedup.
