#ifndef GET_CPU_H
#define GET_CPU_H

// ---------------------------------------------------------------------
// getCPU(): wall-clock time in seconds (the same kind of timer in all codes)
//   OpenMP : omp_get_wtime()
//   MPI    : MPI_Wtime()        (include <mpi.h> before this file)
//   serial : std::chrono::steady_clock
// ---------------------------------------------------------------------
#if defined(_OPENMP)
  #include <omp.h>
  inline double getCPU() { return omp_get_wtime(); }
#elif defined(MPI_VERSION)
  inline double getCPU() { return MPI_Wtime(); }
#else
  #include <chrono>
  inline double getCPU()
  {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
  }
#endif

#endif
