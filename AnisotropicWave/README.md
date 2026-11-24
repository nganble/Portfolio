# Parallel Solver for Anisotropic Wave Equation

## Overview
This project focuses on solving the anisotropic wave equation using finite difference. The work includes:
- A serial implementation in C++
- Parallel solvers developed using OpenMP and MPI for speedup

## Features
- **Numerical Method**: Utilizes second-order finite difference approximations for spatial and time discretization.
- **Code Development**:
  - Serial implementation written in C++
  - Parallel versions developed using OpenMP and MPI
- **Verification**:
  - Tested using polynomial manufactured solutions (expected to match exactly)
  - Verified second-order convergence with exact solutions

## Results
- **Manufactured Solution**: The computed solution matches the exact solution with no error.
- **Exact Solution**: Demonstrates second-order convergence, consistent with the numerical method used.
