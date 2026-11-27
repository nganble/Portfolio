// =====================================================================
// ANISOTROPIC WAVE EQUATION IN TWO DIMENSIONS: serial code
//   rho w_tt = div(m grad w) + f  on (xa,xb) x (ya,yb), Dirichlet boundary conditions,
//   explicit second-order scheme (see docs/anisotropicWaveScheme.tex)
//
// Build:   make                         (default test case: COEFF=CONS, SOLUTION=POLY)
//          make CXXFLAGS="-O3 -Wall -DCOEFF=VAR -DSOLUTION=POLY"
// Run:     ./wave2dSerial -nx=80 -tFinal=0.5
//          ./wave2dSerial -nx=3000 -tFinal=0.1          (timing)
// Options: -nx=<i> -ny=<i> -tFinal=<f> -debug=<i>
//          debug=1: error after the first time step, debug=2: error of U1, debug=3: print m22
// =====================================================================

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "../common/parseCommand.h"

// Define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;
using std::min;

// getCPU(): wall-clock time
#include "../common/getCPU.h"

// Position of point (i1,i2) in an array: rows are stored with the padded stride ld
#define IDX(i1,i2) ((i1)-nd1a+ld*((i2)-nd2a))
// Array layout and macros x, y, up, uc, un, err, rho, m11, m12, m21, m22; allocateVectors, freeVectors
#include "../common/grid.h"
// L_h (the spatial operator) as the macro LH(u, i1, i2)
#include "../common/stencil.h"
// printCase, printResults (header and result lines)
#include "../common/output.h"

int main(int argc, char *argv[])
{
	printf("Usage: wave2dSerial -nx=<i> -ny=<i> -tFinal=<f> -debug=<i>\n");

	int debug = 0;

	Real xa = -1.0, xb = 1.0; // domain is [xa,xb] X [ya,yb]
	Real ya = -1.0, yb = 1.0;

	Real tFinal = 0.5;
	Real cfl = 0.9;
	int nx = 100, ny = nx;

	string line;
	// parse command line
	for( int i=1; i<argc; i++ )
	{
		line = argv[i];
		if( parseCommand( line,"-nx=",nx) ){ny = nx;}
		else if ( parseCommand( line,"-ny=",ny) ){}
		else if( parseCommand( line,"-debug=",debug) ){}
		else if( parseCommand( line, "-tFinal=",tFinal) ){}
	}

	// Coefficients (COEFF) and exact solutions (SOLUTION); included here, inside main,
	// because omega is evaluated at (xa,ya)
	#include "../common/testCases.h"

	const int numGhost = 1;  // one ghost layer (same layout as the MPI code; not used here)
	const int n1a  = 0;
	const int n1b  = n1a + nx;
	const int nd1a = n1a-numGhost;
	const int nd1b = n1b+numGhost;
	const int nd1  = nd1b-nd1a+1;

	const int n2a  = 0;
	const int n2b  = n2a + ny;
	const int nd2a = n2a-numGhost;
	const int nd2b = n2b+numGhost;
	const int nd2  = nd2b-nd2a+1;
	const int ld = paddedStride(nd1); // padded row stride (common/grid.h)
	const int nd   = ld*nd2; // number of stored values per array
	int i1, i2;
	Real dx[2];
	dx[0] = (xb-xa)/nx;
	dx[1] = (yb-ya)/ny;

	// Arrays: grid points, solution at three time levels, error, coefficients
	Real *x_p[2], *u_p[3], *errCPU_p;
	Real *rho_p, *m11_p, *m12_p, *m21_p, *m22_p;
	// Call the allocate function
	allocateVectors(&x_p[0], &x_p[1], &u_p[0], &u_p[1], &u_p[2], &errCPU_p,
									&rho_p, &m11_p, &m12_p, &m21_p, &m22_p, nd);

	// Fill in grid points
	for (i2=nd2a; i2<=nd2b; i2++)
	{
		for (i1=nd1a; i1<=nd1b; i1++)
		{
			x(i1,i2) = xa + (i1-n1a)*dx[0]; // Fill in x
			y(i1,i2) = ya + (i2-n2a)*dx[1]; // Fill in y
		}
	}

	// Fill in rho and matrix m
	for (i2=n2a; i2<=n2b; i2++)
	{
		for (i1=n1a; i1<=n1b; i1++)
		{
			rho(i1,i2) = RHO(x(i1,i2),y(i1,i2));
			m11(i1,i2) = M11(x(i1,i2),y(i1,i2));
			m12(i1,i2) = M12(x(i1,i2),y(i1,i2));
			m21(i1,i2) = M21(x(i1,i2),y(i1,i2));
			m22(i1,i2) = M22(x(i1,i2),y(i1,i2));
		}
	}
	// Check results for rho and m
	if (debug>2)
	{
		for( i2=n2a; i2<=n2b; i2++ )
		{
			printf("m22(%d:%d,%d)=[", n1a, n1b, i2);
			for( i1=n1a; i1<=n1b; i1++ )
			{
				printf("%4.2f ", m22(i1,i2));

			}
			printf("]\n");
		}
		return 0;
	}


	// Time step: dt = cfl * min_ij sqrt(rho) / sqrt( m11/dx^2 + |m12+m21|/(4 dx dy) + m22/dy^2 )
	Real dt = 10.;
	for (i2=n2a; i2<=n2b; i2++)
	{
		for (i1=n1a; i1<=n1b; i1++)
		{
			Real ind = sqrt(rho(i1,i2))/sqrt( m11(i1,i2)/(dx[0]*dx[0])
									+ fabs((m12(i1,i2) + m21(i1,i2)))/(4*dx[0]*dx[1])
									+ m22(i1,i2)/(dx[1]*dx[1]) );
			dt = min(dt, ind);
		}
	}
	dt = cfl*dt;
	int numSteps = ceil(tFinal/dt);
	dt = tFinal/numSteps; // adjust dt to reach the final time

	// Initial conditions: U0 (exact solution at t = 0) and U1 (first step)
	Real t = 0.;
	int prev = 0, cur = 1;
	for (i2=nd2a; i2<=nd2b; i2++)
	{
		for (i1=nd1a; i1<=nd1b; i1++)
		{
			up(i1,i2)= UTRUE(x(i1,i2),y(i1,i2),t);
		}
	}
	// Approximation for U1 at the interior points: Taylor expansion and the PDE
	for (i2=n2a+1; i2<=n2b-1; i2++)
	{
		for (i1=n1a+1; i1<=n1b-1; i1++)
		{
			uc(i1,i2) = up(i1, i2) + dt * UTRUET(x(i1, i2), y(i1, i2), t) + (dt * dt / (2.0 * rho(i1, i2))) * (LH(up, i1, i2) + FORCE(x(i1, i2), y(i1, i2), t) );
		}
	}
	// Update boundary
	for (i2 = n2a; i2 <= n2b; i2++)
	{
		uc(n1a, i2) = UTRUE(x(n1a, i2), y(n1a, i2), dt);
		uc(n1b, i2) = UTRUE(x(n1b, i2), y(n1b, i2), dt);
	}
	for (i1 = n1a; i1 <= n1b; i1++)
	{
		uc(i1, n2a) = UTRUE(x(i1, n2a), y(i1, n2a), dt);
		uc(i1, n2b) = UTRUE(x(i1, n2b), y(i1, n2b), dt);
	}

	// Check the approximation
	if(debug>1)
	{
		Real maxErr = 0;
		Real err;
		for( i2=n2a; i2<=n2b; i2++ )
		{
			for( i1=n1a; i1<=n1b; i1++ )
			{
				err = fabs(uc(i1,i2) - UTRUE(x(i1,i2), y(i1,i2), dt));
				maxErr = max(err, maxErr);
			}
		}
		printf("Debug: check error for the approximation of U1: maxErr=%6.2e\n", maxErr);
		return 0;
	}

	printf("----- Solving Anisotropic Wave Equation in two dimensions ------\n");
	printf(" nx=%d ny=%d dt=%6.2e numSteps=%d tFinal=%6.2f\n", nx, ny,dt, numSteps, tFinal);
	printCase(stdout, coeffName, msName, xa, xb, ya, yb, cfl, tFinal);

	// ---------------------------- TIME-STEPPING LOOP ------------------------------
	Real cpu1 = getCPU();
	for(int n=0; n<numSteps-1; n++)
	{
		t = (n+1)*dt;  // current time
		const int prev = n%3;         // previous time level
		const int cur  = (n+1)%3;     // current time level
		const int next = (n+2)%3;     // next time level
		// Update the interior points
		for(i2=n2a+1; i2<=n2b-1; i2++)
		{
			for(i1=n1a+1; i1<=n1b-1; i1++)
			{
				un(i1,i2) = 2 * uc(i1, i2) - up(i1, i2) + (dt * dt / rho(i1, i2)) * (LH(uc, i1, i2) + FORCE(x(i1, i2), y(i1, i2), t) );
			}
		}

		// --- boundary conditions ---
		for (i2=n2a; i2<=n2b; i2++)
		{
			un(n1a,i2) = UTRUE(x(n1a,i2),y(n1a,i2),t+dt);
			un(n1b,i2) = UTRUE(x(n1b,i2),y(n1b,i2),t+dt);
		}
		for (i1=n1a; i1<=n1b; i1++)
		{
			un(i1,n2a) = UTRUE(x(i1,n2a),y(i1,n2a),t+dt);
			un(i1,n2b) = UTRUE(x(i1,n2b),y(i1,n2b),t+dt);
		}

		// Check error after first time step
		if (debug>0 && n<1)
		{
			Real maxErr = 0;
			Real err;

			for( i2=n2a; i2<=n2b; i2++ )
			{
				for( i1=n1a; i1<=n1b; i1++ )
				{
					Real ue = UTRUE(x(i1,i2), y(i1,i2), t+dt);
					err = fabs(un(i1,i2) - ue);
					maxErr = max(err, maxErr);
					if (nx <= 10)
						printf("err[%d,%d]=%6.2e, uc=%4.7f, ue=%4.7f\n", i1, i2, err, un(i1,i2), ue);
				}
			}
			printf("Debug: after the first time-step:n=%d maxErr=%6.2e\n", n, maxErr);
			return 0;
		}

	}
	// ---------------------- END TIME-STEPPING LOOP ------------------------------
	Real cpuTimeStep = getCPU() - cpu1;

	// --- errors at t = tFinal ---
	t += dt; // tFinal
	if (fabs(t-tFinal) > 1e-3*dt/tFinal){
		printf("ERROR: AFTER TIME_STEPPING: t = %16.8e IS NOT EQUAL TO tFinal = %16.8e\n", t, tFinal);
	}
	cur = numSteps%3;
	Real maxErr = 0.;
	Real maxNorm = 0.;
	for( i2=n2a; i2<=n2b; i2++ )
	{
		for( i1=n1a; i1<=n1b; i1++ )
		{
			err(i1,i2) = fabs(uc(i1,i2) - UTRUE(x(i1,i2), y(i1,i2),tFinal));
			maxErr = max(err(i1,i2), maxErr);
			maxNorm = max(uc(i1,i2), maxNorm);
		}
	}
	Real maxRelErr = maxErr/maxNorm;
	printResults(coeffName, msName, tFinal, nx, ny, numSteps, dt, maxNorm, maxErr, maxRelErr, cpuTimeStep);

	// Free allocated vectors
	freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
							rho_p, m11_p, m12_p, m21_p, m22_p);

	return 0;
}
