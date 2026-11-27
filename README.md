# Portfolio
This repo holds research and course projects in numerical methods for PDEs and parallel computing.

| Project | Description | Languages |
|---|---|---|
| [EigenWave](EigenWave/) | EigenWave algorithm: computing eigenvalues and eigenfunctions of the discrete Laplacian near a target frequency. The algorithm is based on solving an associated wave equation and filtering the solution in time. The code computes eigenpairs on an interval, a square, and an annulus | MATLAB |
| [AnisotropicWave](AnisotropicWave/) | Solving the anisotropic wave equation $\rho w_{tt}=\nabla\cdot(m\nabla w)+f$ with variable coefficients: explicit second-order finite difference scheme; serial, OpenMP, and MPI solvers verified against a MATLAB prototype; speedup of 21–27 on 32 cores | C++, OpenMP, MPI, MATLAB |
| [PoissonSolver2D](PoissonSolver2D/) | Solving the 2D Poisson equation by SOR with red-black ordering: serial, MPI, CUDA, and MPI+CUDA solvers; strong and weak scaling up to 128 MPI ranks and 16 GPUs | C, MPI, CUDA, MATLAB |

See the README in each folder for build instructions, verification, and results.
