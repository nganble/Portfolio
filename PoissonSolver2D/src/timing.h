#ifndef TIMING_H
#define TIMING_H

/* ===================================================
	Configuration
====================================================*/

/* Default settings (can override with -D flags) */

#ifndef USE_MPI
#define USE_MPI 0
#endif

#ifndef USE_POWER9
#define USE_POWER9 0
#endif

/* ===================================================
	Include Headers
====================================================*/
#if USE_MPI
#include <mpi.h>
#else
#include <time.h>
#endif

/* ===================================================
	POWER9 System Timer
====================================================*/

# if USE_POWER9

static inline unsigned long long aimos_clock_read(void)
{
    unsigned int tbl, tbu0, tbu1;
    do {
        __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
        __asm__ __volatile__ ("mftb %0"  : "=r"(tbl));
        __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
    } while (tbu0 != tbu1);

    return (((unsigned long long)tbu0) << 32) | tbl;
}

#endif

/* ===================================================
	Unified Timing Options
====================================================*/

static inline double get_time(){

#if USE_POWER9

	/* Convert to seconds */
	const double freq = 512000000.0;
	return (double) aimos_clock_read() / freq;

#elif USE_MPI

	return MPI_Wtime();

#else

	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1e-9 * ts.tv_nsec;

#endif

}


#endif /* TIMING_H */