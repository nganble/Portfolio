#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>

// =====================================================================
// Output shared by the serial, OpenMP and MPI codes
// =====================================================================

// Header line: test case and run parameters.
// parName/parValue: optional extra parameter, e.g. "numThreads" (OpenMP) or "np" (MPI)
inline void printCase(FILE *file, const char *coeffName, const char *msName,
                      double xa, double xb, double ya, double yb, double cfl, double tFinal,
                      const char *parName = NULL, int parValue = 0)
{
	fprintf(file, " coeffs=%s, ms=%s, domain=[%g,%g]x[%g,%g], cfl=%4.2f, tFinal=%4.2f",
	        coeffName, msName, xa, xb, ya, yb, cfl, tFinal);
	if (parName != NULL)
		fprintf(file, ", %s=%d", parName, parValue);
	fprintf(file, "\n");
}

// Result lines at t = tFinal:
//   - the same format as MATLAB (anisotropicWave_2D.m), plus the cpu time
//   - the BASELINE line in full precision (read by the baseline scripts)
inline void printResults(const char *coeffName, const char *msName, double tFinal,
                         int nx, int ny, int numSteps, double dt,
                         double maxNorm, double maxErr, double maxRelErr, double cpuTime)
{
	printf("%s: t=%4.2f Nx=%3d Ny=%3d Nt=%3d dt=%8.2e maxNorm=%8.2e maxRelErr=%8.2e cpuTime=%9.2e(s)\n",
	       msName, tFinal, nx, ny, numSteps, dt, maxNorm, maxRelErr, cpuTime);
	printf("BASELINE: coeffs=%s ms=%s nx=%d numSteps=%d dt=%.16e maxErr=%.16e maxRelErr=%.16e\n",
	       coeffName, msName, nx, numSteps, dt, maxErr, maxRelErr);
}

#endif
