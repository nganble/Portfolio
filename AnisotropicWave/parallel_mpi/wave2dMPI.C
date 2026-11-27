// =====================================================================
// ANISOTROPIC WAVE EQUATION IN TWO DIMENSIONS: MPI code
//   rho w_tt = div(m grad w) + f  on (xa,xb) x (ya,yb), Dirichlet boundary conditions,
//   explicit second-order scheme (see docs/anisotropicWaveScheme.tex)
//   The grid is split in the y-direction: each rank owns a block of rows plus one
//   ghost line above and below, exchanged with MPI_Sendrecv after each step.
//
// Build:   make                         (default test case: COEFF=CONS, SOLUTION=POLY)
//          make CXXFLAGS="-O3 -Wall -DCOEFF=VAR -DSOLUTION=POLY"
// Run:     mpirun -np 4 ./wave2dMPI -nx=80 -tFinal=0.5
//          mpirun -np 4 ./wave2dMPI -nx=3000 -tFinal=0.1      (timing)
// Options: -nx=<i> -ny=<i> -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] -matlabFileName=<s>
//          saveMatlab=1: write a Matlab file with the parameters, =2: also x, y, u and the error
//          debug=1: write debug/debugFileNp<np>Proc<rank>.debug (the folder is created if needed)
//          debug=2: errors after the first steps, debug=3: error of U1, debug=4: grid points (nx <= 16)
// =====================================================================

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <mpi.h>
#include <sys/stat.h>   // mkdir
#include <errno.h>
// getLocalIndexBounds(): split the rows between the ranks
#include "getLocalIndexBounds.h"
// parseCommand(): read the command line options
#include "../common/parseCommand.h"
#define REAL_MIN DBL_MIN

// Define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;
using std::min;

// getCPU(): wall-clock time (MPI_Wtime)
#include "../common/getCPU.h"

// writeMatlabArray(): write an array to a Matlab readable file
#include "writeMatlabArray.h"

// Position of point (i1,i2) in the local arrays: rows are stored with the padded stride ld,
// counted from the first stored local row nd2a_l
#define IDX(i1, i2) ((i1) - nd1a + ld * ((i2) - nd2a_l))
// Array layout and macros x, y, up, uc, un, err, rho, m11, m12, m21, m22; allocateVectors, freeVectors
#include "../common/grid.h"
// L_h (the spatial operator) as the macro LH(u, i1, i2)
#include "../common/stencil.h"
// printCase, printResults (header and result lines)
#include "../common/output.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	// Determine the rank of the current process
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	// Determine the number of processes
	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// debug = 4: print the grid points (nx <= 16)
	// debug = 3: check the approximation for U1
	// debug = 2: check the error after the first time steps
	// debug = 1: write to the debug files
	int debug = 0;
	Real xa = -1.0, xb = 1.0; // domain is [xa,xb] X [ya,yb]
	Real ya = -1.0, yb = 1.0;

	Real tFinal = 0.5;
	Real cfl = 0.9;
	int nx = 100, ny = nx;
	int saveMatlab = 0; // 1 = save a matlab file, 2 = save solution too
	string matlabFileName = "anisotropicWave2DMPI.m";

	bool echo;
	if (myRank==0)
		echo = true;
	else
		echo = false;
	string line;
	// parse command line
	for( int i=1; i<argc; i++ )
	{
		line = argv[i];
		if( parseCommand( line, "-nx=", nx, echo) ){ny = nx;}
		else if ( parseCommand( line, "-ny=", ny, echo) ){}
		else if( parseCommand( line, "-debug=", debug, echo) ){}
		else if( parseCommand( line, "-tFinal=", tFinal, echo) ){}
		else if( parseCommand( line, "-saveMatlab=", saveMatlab, echo) ){}
		else if( parseCommand( line, "-matlabFileName=", matlabFileName, echo) ){}
	}

	// Debug files: one per rank, in the folder debug/
	FILE *debugFile = NULL;
	if (debug > 0)
	{
		// rank 0 creates the folder debug/ if it does not exist yet
		if (myRank == 0 && mkdir("debug", 0755) != 0 && errno != EEXIST)
			printf("Warning: could not create the folder debug/\n");
		MPI_Barrier(MPI_COMM_WORLD);   // the folder exists before any rank opens its file
		char debugFileName[80];
		snprintf(debugFileName, sizeof(debugFileName), "debug/debugFileNp%dProc%d.debug", np, myRank);
		debugFile = fopen(debugFileName, "w");
		if (debugFile == NULL)          // checked on every rank
		{
			printf("Error: rank %d cannot open %s\n", myRank, debugFileName);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		fprintf(debugFile, "------------------- DebugFile -------------------\n");
		fprintf(debugFile, "------------------- np=%d, myRank=%d-------------\n\n", np, myRank);
	}
	if (myRank == 0)
		printf("Usage: wave2dMPI -nx=<i> -ny=<i> -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] -matlabFileName=<s>\n");

	// Coefficients (COEFF) and exact solutions (SOLUTION); included here, inside main,
	// because omega is evaluated at (xa,ya)
	#include "../common/testCases.h"

	const int numGhost = 1;  // one ghost layer (holds the parallel ghost lines)
	const int n1a  = 0;
	const int n1b  = n1a + nx;
	const int nd1a = n1a - numGhost;
	const int nd1b = n1b + numGhost;
	const int nd1  = nd1b - nd1a+1;

	const int n2a  = 0;
	const int n2b  = n2a + ny;
	const int nd2a = n2a - numGhost;
	const int nd2b = n2b + numGhost;
	const int nd2  = nd2b - nd2a + 1;
	const int ld = paddedStride(nd1); // padded row stride (common/grid.h)
	const int nd   = ld * nd2; // number of stored values of a global array
	int i1, i2;
	Real dx[2];
	dx[0] = (xb - xa) / nx;
	dx[1] = (yb - ya) / ny;

	// Get local index bound in the y direction
	int ny_l, n2a_l, n2b_l;
	getLocalIndexBounds(myRank, np, ny, ny_l, n2a_l, n2b_l);
	const int nd2a_l = n2a_l - numGhost;
	const int nd2b_l = n2b_l + numGhost;
	const int nd2_l  = nd2b_l - nd2a_l+1;
	const int nd_l   = ld * nd2_l; // number of stored values of a local array
	// Interior rows owned by this rank (the global boundary rows are set by the boundary conditions)
	const int n2a_int = max(n2a_l, n2a + 1);
	const int n2b_int = min(n2b_l, n2b - 1);

	// Local arrays: grid points, solution at three time levels, error, coefficients
	Real *x_p[2], *u_p[3], *errCPU_p;
	Real *rho_p, *m11_p, *m12_p, *m21_p, *m22_p;
	// Call the allocate function
	allocateVectors(&x_p[0], &x_p[1], &u_p[0], &u_p[1], &u_p[2], &errCPU_p,
									&rho_p, &m11_p, &m12_p, &m21_p, &m22_p, nd_l);

	// Fill in local grid points
	for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
	{
		for (i1 = nd1a; i1 <= nd1b; i1++)
		{
			x(i1, i2) = xa + (i1 - n1a) * dx[0]; // Fill in x
			y(i1, i2) = ya + (i2 - n2a) * dx[1]; // Fill in y
		}
	}
	// Check if the local grid points is correct
	if (debug > 3 && nx <= 16)
	{
		for (int ifile = 0; ifile <= 1; ifile++)
		{
			FILE *file = ifile == 0 ? stdout : debugFile;
			if ((ifile == 0 && myRank == 0) ||  // write to terminal if myrank=0
					(ifile == 1))                   // write to debug file
			{
				// Print out vector x
				fprintf (file, "nx=%d [n1a,n1b]=[%d,%d] [nd1a,nd1b]=[%d,%d]\n",
									nx, n1a, n1b, nd1a, nd1b);
				for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
				{
					fprintf (file, "x(%d:%d,%d)=[", nd1a, nd1b, i2);
					for (i1 = nd1a; i1 <= nd1b; i1++)
					{
						fprintf (file, "%4.2f ", x(i1, i2));

					}
					fprintf (file, "]\n");
				}

				// Print out vector y
				fprintf (file, "ny_l=%d [n2a_l,n2b_l]=[%d,%d] [nd2a_l,nd2b_l]=[%d,%d]\n",
									ny_l, n2a_l, n2b_l, nd2a_l, nd2b_l);
				for ( i2 = nd2a_l; i2 <= nd2b_l; i2++ )
				{
					fprintf(file, "y(%d:%d,%d)=[", nd1a, nd1b, i2);
					for( i1 = nd1a; i1 <= nd1b; i1++ )
					{
						fprintf(file, "%4.2f ", y(i1, i2));

					}
					fprintf(file, "]\n");
				}

			}

		}

		// Finish the program (optional)
		freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
								rho_p, m11_p, m12_p, m21_p, m22_p);
		if (debugFile)
			fclose(debugFile);
		MPI_Finalize();
		return 0;

	} // End for checking local grid points

	// Fill in rho and matrix m
	for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
	{
		for (i1 = n1a; i1 <= n1b; i1++)
		{
			rho(i1,i2) = RHO(x(i1,i2),y(i1,i2));
			m11(i1,i2) = M11(x(i1,i2),y(i1,i2));
			m12(i1,i2) = M12(x(i1,i2),y(i1,i2));
			m21(i1,i2) = M21(x(i1,i2),y(i1,i2));
			m22(i1,i2) = M22(x(i1,i2),y(i1,i2));
		}
	}


	// Time step: dt = cfl * min_ij sqrt(rho) / sqrt( m11/dx^2 + |m12+m21|/(4 dx dy) + m22/dy^2 )
	Real dt_l = 10.;
	for (i2 = n2a_l; i2 <= n2b_l; i2++) // Don't need to compute at ghost points
	{
		for (i1 = n1a; i1 <= n1b; i1++)   // Don't need to compute at ghost points
		{
			Real ind = sqrt(rho(i1, i2)) / sqrt(m11(i1,i2) / (dx[0] * dx[0])
									+ fabs(m12(i1, i2) + m21(i1,i2))/(4 * dx[0] * dx[1])
									+ m22(i1, i2) / (dx[1] * dx[1]));
			dt_l = min(dt_l, ind);
		}
	}

	// Gather dt from all processors to find the smallest one
	Real dt = dt_l;
	MPI_Allreduce(&dt_l, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	dt = cfl * dt;
	int numSteps = ceil(tFinal / dt);
	dt = tFinal / numSteps; // adjust dt to reach the final time

	// Initial conditions: U0 (exact solution at t = 0) and U1 (first step)
	Real t = 0.;
	int prev = 0, cur = 1;
	for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
	{
		for (i1 = nd1a; i1 <= nd1b; i1++)
		{
			up(i1, i2) = UTRUE(x(i1, i2), y(i1, i2), t);
		}
	}
	// Approximation for U1 at the interior points: Taylor expansion and the PDE
	for (i2 = n2a_int; i2 <= n2b_int; i2++)
	{
		for (i1 = n1a + 1; i1 <= n1b - 1; i1++)
		{
			uc(i1,i2) = up(i1, i2) + dt * UTRUET(x(i1, i2), y(i1, i2), t) + (dt * dt / (2.0 * rho(i1, i2))) * (LH(up, i1, i2) + FORCE(x(i1, i2), y(i1, i2), t) );
		}
	}
	// Update boundary
	for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
	{
		uc(n1a, i2) = UTRUE(x(n1a, i2), y(n1a, i2), dt); // Left
		uc(n1b, i2) = UTRUE(x(n1b, i2), y(n1b, i2), dt); // Right
	}
	if (np == 1) // only one processor
	{
		for (i1 = nd1a; i1 <= nd1b; i1++)
		{
			uc(i1, n2a_l) = UTRUE(x(i1, n2a_l), y(i1, n2a_l), dt); // Bottom
			uc(i1, n2b_l) = UTRUE(x(i1, n2b_l), y(i1, n2b_l), dt); // Top
		}
	}
	else  // we have more than one processor
	{
		// Dirichlet bc's at the bottom (rank 0) and top (rank np-1) boundaries
		if (myRank == 0)
			for (i1 = nd1a; i1 <= nd1b; i1++)
				uc(i1, n2a_l) = UTRUE(x(i1, n2a_l), y(i1, n2a_l), dt); // Bottom
		if (myRank == np-1)
			for (i1 = nd1a; i1 <= nd1b; i1++)
				uc(i1, n2b_l) = UTRUE(x(i1, n2b_l), y(i1, n2b_l), dt); // Top
		// Exchange parallel ghost lines with the neighbours (MPI_PROC_NULL at the ends).
		// MPI_Sendrecv does not rely on MPI buffering, so it cannot deadlock.
		const int below = (myRank > 0)      ? myRank - 1 : MPI_PROC_NULL;
		const int above = (myRank < np - 1) ? myRank + 1 : MPI_PROC_NULL;
		// send the first line to rank-1, receive the top ghost line from rank+1
		MPI_Sendrecv(&uc(n1a, n2a_l), nx + 1, MPI_DOUBLE, below, 1,
		             &uc(n1a, nd2b_l), nx + 1, MPI_DOUBLE, above, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// send the last line to rank+1, receive the bottom ghost line from rank-1
		MPI_Sendrecv(&uc(n1a, n2b_l), nx + 1, MPI_DOUBLE, above, 1,
		             &uc(n1a, nd2a_l), nx + 1, MPI_DOUBLE, below, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Check error for approximation of U1
	if (debug > 2)
	{
		Real maxErr = 0.0;
		Real maxNorm = 0.0;
		Real *testErr_p = new Real [nd1*nd2_l];
		#define testErr(i1, i2) testErr_p[(i1 - nd1a) + nd1 * (i2 - nd2a_l)]
		for( i2 = n2a_l; i2 <= n2b_l; i2++ )
		{
			for( i1 = n1a; i1 <= n1b; i1++ )
			{
				testErr(i1, i2) = fabs(uc(i1,i2) - UTRUE(x(i1,i2), y(i1,i2), dt));
				maxErr = max(testErr(i1, i2), maxErr);
				maxNorm = max(uc(i1, i2), maxNorm);
			}
		}
		// Gather maxNorm in all ranks
		Real maxNorm_g = maxNorm;
		MPI_Allreduce(&maxNorm, &maxNorm_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		maxErr /= max(maxNorm_g,REAL_MIN); // relative error
		// Gather max error in all ranks
		Real maxErr_g = maxErr;
		MPI_Allreduce(&maxErr, &maxErr_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		for (int ifile = 0; ifile <= 1; ifile++)
		{
			FILE *file = ifile == 0 ? stdout : debugFile;
			if ((ifile == 0 && myRank == 0) ||  // write to terminal if myRank==0
					(ifile==1))											// //write to debug file
			{
				fprintf(file, "\nDebug: check error for the approximation of U1: t=%9.2e\n", t+dt);
				if (nx<=16){
					for( i2 = n2a_l; i2 <= n2b_l; i2++ )
					{
						fprintf(file, "error(%d:%d,%d)=[", n1a, n1b, i2);
						for( i1 = n1a; i1 <= n1b; i1++ )
						{
							fprintf(file, "%9.2e ", testErr(i1,i2));
						}
						fprintf(file, "]\n");
					}
				}

				fprintf(file, "On myRank=%d: maxNormL=%9.2e maxRelErrL=%9.2e,  Gathered solution: maxNormG=%9.2e maxRelErrG=%9.2e\n",
					myRank, maxNorm, maxErr, maxNorm_g, maxErr_g);

			}
		}
		delete [] testErr_p;

		// Finish the program (optional)
		freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
								rho_p, m11_p, m12_p, m21_p, m22_p);
		if (debugFile)
			fclose(debugFile);
		MPI_Finalize();
		return 0;

	} // End for debugging the approximation of U1

	for (int ifile = 0; ifile <= 1; ifile++)
	{
		FILE *file = ifile == 0 ? stdout : debugFile;
		if ((ifile == 0 && myRank == 0) || // write to terminal if myRank==0
			 (ifile == 1 && debug > 0))    	 // write to debugFile if debug>0
		{
			fprintf(file, "----- Solving Anisotropic Wave Equation in two dimensions ------\n");
			fprintf(file, " saveMatlab=%d, matlabFileName=%s \n",saveMatlab,matlabFileName.c_str());
			fprintf(file, " nx=%d ny=%d dt=%6.2e numSteps=%d tFinal=%6.2f\n", nx, ny,dt, numSteps, tFinal);
			printCase(file, coeffName, msName, xa, xb, ya, yb, cfl, tFinal, "np", np);
		}
	}

	// ---------------------------- TIME-STEPPING LOOP ------------------------------
	Real cpu1 = getCPU();
	for(int n=0; n<numSteps-1; n++)
	{
		t = (n + 1) * dt;  // current time
		const int prev = n % 3;         // previous time level
		const int cur  = (n + 1) % 3;     // current time level
		const int next = (n + 2) % 3;     // next time level

		// Update the interior points
		for(i2 = n2a_int; i2 <= n2b_int; i2++)
		{
			for(i1 = n1a + 1; i1 <= n1b - 1; i1++)
			{
				un(i1,i2) = 2 * uc(i1, i2) - up(i1, i2) + (dt * dt / rho(i1, i2)) * (LH(uc, i1, i2) + FORCE(x(i1, i2), y(i1, i2), t) );
			}
		}

		// --- boundary conditions ---
		for (i2 = nd2a_l; i2 <= nd2b_l; i2++)
		{
			un(n1a, i2) = UTRUE(x(n1a, i2), y(n1a, i2), t + dt); // Left
			un(n1b, i2) = UTRUE(x(n1b, i2), y(n1b, i2), t + dt); // Right
		}
		if (np == 1) // only one processor
		{
			for (i1 = nd1a; i1 <= nd1b; i1++)
			{
				un(i1, n2a_l) = UTRUE(x(i1, n2a_l), y(i1, n2a_l), t + dt); // Bottom
				un(i1, n2b_l) = UTRUE(x(i1, n2b_l), y(i1, n2b_l), t + dt); // Top
			}
		}
		else				 // more than one processor
		{
			// Dirichlet bc's at the bottom (rank 0) and top (rank np-1) boundaries
			if (myRank == 0)
				for (i1 = nd1a; i1 <= nd1b; i1++)
					un(i1, n2a_l) = UTRUE(x(i1, n2a_l), y(i1, n2a_l), t + dt); // Bottom
			if (myRank == np-1)
				for (i1 = nd1a; i1 <= nd1b; i1++)
					un(i1, n2b_l) = UTRUE(x(i1, n2b_l), y(i1, n2b_l), t + dt); // Top
			// Exchange parallel ghost lines with the neighbours (MPI_PROC_NULL at the ends).
			// MPI_Sendrecv does not rely on MPI buffering, so it cannot deadlock.
			const int below = (myRank > 0)      ? myRank - 1 : MPI_PROC_NULL;
			const int above = (myRank < np - 1) ? myRank + 1 : MPI_PROC_NULL;
			// send the first line to rank-1, receive the top ghost line from rank+1
			MPI_Sendrecv(&un(n1a, n2a_l), nx + 1, MPI_DOUBLE, below, n,
			             &un(n1a, nd2b_l), nx + 1, MPI_DOUBLE, above, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// send the last line to rank+1, receive the bottom ghost line from rank-1
			MPI_Sendrecv(&un(n1a, n2b_l), nx + 1, MPI_DOUBLE, above, n,
			             &un(n1a, nd2a_l), nx + 1, MPI_DOUBLE, below, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Check error after first time step
		if (debug > 1 && n<=2)
		{
			Real maxErr = 0.0;
			Real maxNorm = 0.0;
			Real *testErr_p = new Real [nd1*nd2_l];
			#define testErr(i1, i2) testErr_p[(i1 - nd1a) + nd1 * (i2 - nd2a_l)]
			for( i2 = n2a_l; i2 <= n2b_l; i2++ )
			{
				for( i1 = n1a; i1 <= n1b; i1++ )
				{
					testErr(i1, i2) = fabs(un(i1, i2) - UTRUE(x(i1, i2), y(i1, i2), t + dt));
					maxErr = max(testErr(i1, i2), maxErr);
					maxNorm = max(un(i1, i2), maxNorm);
				}
			}
			// Gather maxNorm in all ranks
			Real maxNorm_g = maxNorm;
			MPI_Allreduce(&maxNorm, &maxNorm_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			maxErr /= max(maxNorm_g, REAL_MIN); // relative error
			// Gather max error in all ranks
			Real maxErr_g = maxErr;
			MPI_Allreduce(&maxErr, &maxErr_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			for (int ifile = 0; ifile <= 1; ifile++)
			{
				FILE *file = ifile == 0 ? stdout : debugFile;
				if ((ifile == 0 && myRank == 0) ||  // write to terminal if myRank==0
						(ifile==1))											// //write to debug file
				{
					fprintf(file, "\nDebug: check error after one step: n=%d\n", n);
					if (nx<=16){
						for( i2 = n2a_l; i2 <= n2b_l; i2++ )
						{
							fprintf(file, "error(%d:%d,%d)=[", n1a, n1b, i2);
							for( i1 = n1a; i1 <= n1b; i1++ )
							{
								fprintf(file, "%9.2e ", testErr(i1,i2));
							}
							fprintf(file, "]\n");
						}
					}

					fprintf(file, "On myRank=%d: maxNormL=%9.2e maxRelErrL=%9.2e,  Gathered solution: maxNormG=%9.2e maxRelErrG=%9.2e\n",
													myRank, maxNorm, maxErr, maxNorm_g, maxErr_g);

				}
			}

			if (n == numSteps-2)
			{
				delete [] testErr_p;
				// Finish the program (optional)
				freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
										rho_p, m11_p, m12_p, m21_p, m22_p);
				if (debugFile)
					fclose(debugFile);
				MPI_Finalize();
				return 0;
			}


		} // End for debug after first time-step

	}
	// ---------------------- END TIME-STEPPING LOOP ------------------------------
	Real cpuTimeStep = getCPU() - cpu1;

	// --- compute errors at tFinal ---
	t += dt; // tFinal
	if ((fabs(t - tFinal) > 1e-3 * dt / tFinal) && myRank == 0){
		printf("ERROR: AFTER TIME_STEPPING: t = %16.8e IS NOT EQUAL TO tFinal = %16.8e\n", t, tFinal);
	}
	cur = numSteps % 3;
	Real maxErr = 0.0;
	Real maxNorm = 0.0;
	for (i2 = n2a_l; i2 <= n2b_l; i2++)
	{
		for (i1 = n1a; i1 <= n1b; i1++)
		{
			err(i1, i2) = fabs(uc(i1, i2) - UTRUE(x(i1, i2), y(i1, i2), tFinal));
			maxErr = max(err(i1, i2), maxErr);
			maxNorm = max(uc(i1, i2), maxNorm);
		}
	}

	// Gather maxNorm in all ranks
	Real maxNorm_g = maxNorm;
	MPI_Allreduce(&maxNorm, &maxNorm_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	// Gather the absolute max error in all ranks (for the baseline)
	Real maxErrAbs_g = maxErr;
	MPI_Allreduce(&maxErr, &maxErrAbs_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	maxErr /= max(maxNorm_g, REAL_MIN); // relative error
	// Gather max error in all ranks
	Real maxErr_g = maxErr;
	MPI_Allreduce(&maxErr, &maxErr_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	// Gather max cpuTime in all ranks
	Real cpuTimeStep_g = cpuTimeStep;
	MPI_Allreduce(&cpuTimeStep, &cpuTimeStep_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for (int ifile = 0; ifile <= 1; ifile++)
	{
		FILE *file = ifile == 0 ? stdout : debugFile;
		if ((ifile == 0 && myRank == 0) ||  // write to terminal if myRank==0
				(ifile == 1 && debug > 0))		  // write to debug file
		{
			fprintf(file, "\nAfter time-stepping loop: n=%d\n", numSteps);
			fprintf(file, "numSteps=%d nx=%d maxNormL=%8.2e maxRelErrL=%8.2e\n cpuTimeStepL=%9.2e(s)\n",
							numSteps, nx, maxNorm, maxErr, cpuTimeStep);
			fprintf(file, "np=%d myRank=%d Gathered solution: maxNormG=%8.2e maxRelErrG=%8.2e cpuTimeStepG=%9.2e(s)\n",
					 		np, myRank, maxNorm_g, maxErr_g, cpuTimeStep_g);
		}
	}
	if (myRank == 0)
	{
		// Result lines (global values; cpu time = max over ranks)
		printResults(coeffName, msName, tFinal, nx, ny, numSteps, dt, maxNorm_g, maxErrAbs_g, maxErr_g, cpuTimeStep_g);
	}

	// Gather all data from various processes to rank 0
	if (myRank == 0)
	{
		// Allocate space to hold grid points
		Real *xg_p = new Real[nd];
		Real *yg_p = new Real[nd];
		// Allocate space to hold the solution
		Real *ug_p = new Real[nd];
		// Allocate space to hold error
		Real *errg_p = new Real[nd];
		#define xg(i1, i2) xg_p[(i1 - nd1a) + ld * (i2 - nd2a)]
		#define yg(i1, i2) yg_p[(i1 - nd1a) + ld * (i2 - nd2a)]
		#define ug(i1, i2) ug_p[(i1 - nd1a) + ld * (i2 - nd2a)]
		#define errg(i1, i2) errg_p[(i1 - nd1a) + ld * (i2 - nd2a)]
		// Copy the local data on rank 0 to the global arrays
		for (i2=nd2a_l; i2<=nd2b_l; i2++)
		{
			for (i1=nd1a; i1<=nd1b; i1++)
			{
				xg(i1, i2) = x(i1, i2);
				yg(i1, i2) = y(i1, i2);
				ug(i1, i2) = uc(i1, i2);
				errg(i1, i2) = err(i1, i2);
			}
		}
		// Request the data from each of the other rank, appending as we go.
		for (int rank = 1; rank <= np-1; rank++)
		{
			// Need to get local index for each rank in the y-direction, since we are now on rank 0
			int ny_local, n2a_local, n2b_local;
			getLocalIndexBounds(rank, np, ny, ny_local, n2a_local, n2b_local);
			int dataLen;
			if (rank < np-1)
				dataLen = ld*(ny_local+1);
			else
				dataLen = ld*(ny_local+2); // the last rank also sends its top ghost line
			// Receive data from other rank
			MPI_Recv(&(xg(nd1a,n2a_local)), dataLen, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&(yg(nd1a,n2a_local)), dataLen, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&(ug(nd1a,n2a_local)), dataLen, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&(errg(nd1a,n2a_local)), dataLen, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		// Print the gathered error (debug, small grids only)
		if (debug > 0 && nx <= 10)
		{
			printf("After gathering:\n");
			for( i2 = n2a; i2 <= n2b; i2++)
			{
				printf("errGlobal(%d:%d,%d)=[", n1a, n1b, i2);
				for( i1 = n1a; i1 <= n1b; i1++)
				{
					printf("%9.2e ", errg(i1, i2));
				}
				printf("]\n");
			}
		}
		// --- OPTIONALLY write a matlab file for plotting in matlab ---
		if (saveMatlab)
		{
			FILE *matlabFile = fopen(matlabFileName.c_str(), "w");
			fprintf(matlabFile, "%% File written by wave2dMPI.C\n");
			fprintf(matlabFile,"xa=%g; xb=%g; ya=%g; yb=%g; t=%g; maxErr=%10.3e; cpuTimeStep=%10.3e;\n",
				xa, xb, ya, yb, tFinal, maxErr_g, cpuTimeStep_g);
			fprintf(matlabFile, "n1a=%d; n1b=%d; nd1a=%d; nd1b=%d;\n", n1a, n1b, nd1a, nd1b);
			fprintf(matlabFile, "n2a=%d; n2b=%d; nd2a=%d; nd2b=%d;\n", n2a, n2b, nd2a, nd2b);
			fprintf(matlabFile, "dx=%14.6e; dy=%14.6e; numGhost=%d;\n", dx[0], dx[1], numGhost);

			if (saveMatlab > 1)
			{
				writeMatlabArray( matlabFile, xg_p, "x", nd1a, nd1b, nd2a, nd2b, ld);
				writeMatlabArray( matlabFile, yg_p, "y", nd1a, nd1b, nd2a, nd2b, ld);
				writeMatlabArray( matlabFile, ug_p, "u", nd1a, nd1b, nd2a, nd2b, ld);
				writeMatlabArray( matlabFile, errg_p, "err", nd1a, nd1b, nd2a, nd2b, ld);
			}
			fclose(matlabFile);
			printf("Wrote file [%s]\n",matlabFileName.c_str());
		}

		delete [] xg_p;
		delete [] yg_p;
		delete [] ug_p;
		delete [] errg_p;
	} // End statement for rank 0
	else
	{
		// All ranks other than 0 must send their interior data to rank 0
		int dataLen;
		if (myRank < np-1)
			dataLen = ld*(ny_l+1);
		else
			dataLen = ld*(ny_l+2);
		MPI_Send(&(x(nd1a,n2a_l)), dataLen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(y(nd1a,n2a_l)), dataLen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(uc(nd1a,n2a_l)), dataLen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(err(nd1a,n2a_l)), dataLen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
							rho_p, m11_p, m12_p, m21_p, m22_p);

	if (debugFile)
		fclose(debugFile);
	MPI_Finalize();
	return 0;
}
