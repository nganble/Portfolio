// g++ -O3 -o wave2dSerial main_serial.C
// ./main_serial -nx=80 -tFinal=0.5 -debug=0
// ./main_serial -nx=2048 -tFinal=0.1 -debug=0 (test for speedup)
// =====================================================================
//
// ANISOTROPIC WAVE EQUATION IN TWO DIMENSIONS SERIAL VERSION 
//
// =====================================================================
// Solver the anisotropic equation in two-dimension

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "parse_command.h"

// Define a new type "Real" which is equivalent to a "double"
typedef double Real;

#include <string>
using std::string;
using std::max;
using std::min;

#include <ctime>
// Function return the current wall-clock time in seconds
inline double getCPU()
{
	return( 1.0*std::clock() )/CLOCKS_PER_SEC;
}


// function to write an array to a matlab reabable file:
// #include "writeMatlabArray.h"


// Macros to use later
#define x(i1,i2) x_p[0][(i1-nd1a)+nd1*(i2-nd2a)]
#define y(i1,i2) x_p[1][(i1-nd1a)+nd1*(i2-nd2a)]
#define up(i1,i2) u_p[prev][(i1-nd1a)+nd1*(i2-nd2a)]
#define uc(i1,i2) u_p[cur ][(i1-nd1a)+nd1*(i2-nd2a)]
#define un(i1,i2) u_p[next][(i1-nd1a)+nd1*(i2-nd2a)]
#define err(i1,i2) errCPU_p[(i1-nd1a)+nd1*(i2-nd2a)]
#define rho(i1,i2) rho_p[(i1-nd1a)+nd1*(i2-nd2a)]
#define m11(i1,i2) m11_p[(i1-nd1a)+nd1*(i2-nd2a)]
#define m12(i1,i2) m12_p[(i1-nd1a)+nd1*(i2-nd2a)]
#define m21(i1,i2) m21_p[(i1-nd1a)+nd1*(i2-nd2a)]
#define m22(i1,i2) m22_p[(i1-nd1a)+nd1*(i2-nd2a)]

// Function to allocate vectors
void allocateVectors( 
			double** x_p, double** y_p,
			double** up_p, double** un_p,
			double** uc_p, double** errCPU,
			double** rho_p, double** m11_p,
			double** m12_p, double** m21_p,
			double** m22_p, int nd);

// Function to free all allocated vectors
void freeVectors( 
			double* x_p, double* y_p,
			double* up_p, double* un_p,
			double* uc_p, double* errCPU_p,
			double* rho_p, double* m11_p,
			double* m12_p, double* m21_p,
			double* m22_p);


int main(int argc, char *argv[])
{
	const Real pi = 4.*atan2(1.,1.);

	printf("Usage: wave2d -nx=<i> -tFinal=<f> -debug=<i> -saveMatlab=[0|1|2] -"
			"matlabFile=<s>\n");
	
	const int numberOfDimensions = 2;
	int debug = 0;
	// Real kappa = .1;
	Real xa = 0.0, xb = 1.0; // domain is [xa,xb] X [ya,yb]
	Real ya = 0.0, yb = 1.0;

	Real tFinal = 0.5;
	Real cfl = 0.9;
	int nx = 100, ny = nx;
	int saveMatlab = 0; // 1 = save a matlab file, 2 = save solution too
	string matlabFileName = "wave2d.m";

	string line;
	// parse command line
	for( int i=1; i<argc; i++ )
	{
		line = argv[i];
		if( parseCommand( line,"-nx=",nx) ){ny = nx;}
		else if ( parseCommand( line,"-ny=",ny) ){}
		else if( parseCommand( line,"-debug=",debug) ){}
		else if( parseCommand( line, "-tFinal=",tFinal) ){}
		else if( parseCommand( line,"-saveMatlab=",saveMatlab) ){}
		else if( parseCommand( line,"-matlabFileName=",matlabFileName) ){}
	}

	// Define the option for coefficients
	#define CONS  1
	#define VAR   2
	#ifndef COEFF
		#define COEFF CONS
		//#define COEFF VAR
	#endif
	
	#if COEFF == CONS
		#define RHO(x, y)   (1.0)
		#define M11(x, y)   (1.0)
		#define M11x(x, y)  (0.0)
		#define M12(x, y)   (1.0)
		#define M12x(x, y)  (0.0)
		#define M21(x, y)   (0.0)
		#define M21y(x, y)  (0.0)
		#define M22(x, y)   (2.0)
		#define M22y(x, y)  (0.0)

	#elif COEFF == VAR
		const Real a11 = 0.3;
		const Real a12 = 0.2;
		const Real a21 = 0.0;
		const Real a22 = 0.6;
		#define RHO(x, y)  (1.0 * sin((x) + (y)) + 1.1)
		#define M11(x, y)  (a11 * ((x) + (y)) + 1.0)
		#define M11x(x, y) (a11)
		#define M12(x, y)  (0.0)
		#define M12x(x, y) (0.0)
		#define M21(x, y)  (2.0)
		#define M21y(x, y) (0.0)
		#define M22(x, y)  (a22 * (2.0 * (x) + (y)) + 2.0)
		#define M22y(x, y) (a22)

	#endif

	// Define ms solution
	#define TRIG  1
	#define POLY  2
	#define PULSE 3
	#define PULSETRIG 4
	#ifndef SOLUTION
		//#define SOLUTION TRIG
		#define SOLUTION POLY
	#endif

	#if SOLUTION == TRIG
		const Real kx = 1;
		const Real ky = 0;
		const Real omega = sqrt(M11(xa,ya)*kx*kx + (M12(xa,ya) + M21(xa,ya))*kx*ky +
														M22(xa,ya)*ky*ky)/sqrt(RHO(xa,ya));
		#define UTRUE(x, y, t)   (3.0 * cos(kx * (x) + ky * (y) - omega * (t)))
		#define UTRUET(x, y, t)  (3.0 * omega * sin(kx * (x) + ky * (y) - omega * (t)))
		#define FORCE(x, y, t)   (0.0)
	#elif SOLUTION == POLY
		// Coefficients for manufactured solution
		static const Real a0 = .5, a1 = .7, a2 = .9;  // Time coefficients
		static const Real b0 = .5, b1 = .8, b2 = .6;  // Space coefficients in x
		static const Real c0 = .3, c1 = .25, c2 = .2; // Space coefficients in y
		#define UTRUE(x,y,t)   ( (b0 + (x)*( b1 + (x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUET(x,y,t)  ( ( b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*( a1 + 2.*(t)*a2) )
		#define UTRUEX(x,y,t)  ( ( b1 + 2.*(x)*b2)*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEY(x,y,t)  ( (b0 + (x)*( b1 + (x)*b2 ))*(c1 + 2.*(y)*c2)*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUETT(x,y,t) ( (b0 + (x)*( b1 +(x)*b2 ))*(c0 + (y)*( c1 + (y)*c2 ))*(2.*a2 + 0.*(t)) )
		#define UTRUEXY(x,y,t) ( (b1 + 2.*(x)*b2)*(c1 + 2.*(y)*c2)*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEXX(x,y,t) ( ( 2.*b2 )*(c0 + (y)*( c1 + (y)*c2 ))*( a0 + (t)*( a1 + (t)*a2 )) )
		#define UTRUEYY(x,y,t) ( (b0 + (x)*( b1 + (x)*b2 ))*( 2.*c2 )*( a0 + (t)*( a1 + (t)*a2 )) )
		#define FORCE(x,y,t)   (RHO(x, y) * UTRUETT(x, y, t) -\
													 (M11(x, y) * UTRUEXX(x, y, t) + M22(x, y) * UTRUEYY(x, y, t)) - \
												 	 (M12(x, y) + M21(x, y)) * UTRUEXY(x, y, t) - \
												 	 (M11x(x, y) + M21y(x, y)) * UTRUEX(x, y, t) - \
												 	 (M12x(x, y) + M22y(x, y)) * UTRUEY(x, y, t))
	#endif

	const int numGhost = 0;  // use Dirichlet bc's
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
	const int nd   = nd1*nd2; // total number of grid points
	int i1, i2;
	Real dx[2];
	dx[0] = (xb-xa)/nx;
	dx[1] = (yb-ya)/ny;

	/*
	int *gridIndex_p = new int[2*numberOfDimensions];
	#define gridIndexRange(side,axis) gridIndex_p[(side)+(axis)*numberOfDimensions]
	gridIndexRange(0,0) = n1a; gridIndexRange(1,0) = n1b;
	gridIndexRange(0,1) = n2a; gridIndexRange(1,1) = n2b;
	*/

	// int numIntPar = 12;
	// int numRealPar = 15;

	// Declare vectors
	Real *x_p[2], *u_p[3], *errCPU_p;
	// Real *rpar_p; int *ipar_p;
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
	

	// Time-step restriction: // Need to check again
	Real dt = 10.;
	for (i2=n2a; i2<=n2b; i2++)
	{
		for (i1=n1a; i1<=n1b; i1++)
		{
			/*
			Real ind = sqrt(RHO(x(i1,i2),y(i1,i2)))/sqrt( M11(x(i1,i2),y(i1,i2))/(dx[0]*dx[0])
									+ (M12(x(i1,i2),y(i1,i2))+M21(x(i1,i2),y(i1,i2)))/(4*dx[0]*dx[1])
									+ M22(x(i1,i2),y(i1,i2))/(dx[1]*dx[1]) );

			// printf("min[%d,%d]=%4.5f\n", i1, i2,ind);
			*/
			Real ind = sqrt(rho(i1,i2))/sqrt( m11(i1,i2)/(dx[0]*dx[0])
									+ (m12(i1,i2) + m21(i1,i2))/(4*dx[0]*dx[1])
									+ m22(i1,i2)/(dx[1]*dx[1]) );
			// printf("min[%d,%d]=%4.5f\n", i1, i2,ind);
			dt = min(dt, ind);
		}
	}
	dt = cfl*dt;
	int numSteps = ceil(tFinal/dt);
	dt = tFinal/numSteps; // adjust dt to reach the final time
	// printf("dx=%4.5f dy=%4.5f dt=%4.5f numSteps=%d\n", dx[0], dx[1], dt, numSteps);
	

	// Vectors to save integer parameters
	/*
	ipar_p[0]  = n1a; ipar_p[1] = n1b; ipar_p[2] = nd1a; ipar_p[3] = nd1b; ipar_p[4] = nd1;
	ipar_p[5]  = n2a; ipar_p[6] = n2b; ipar_p[7] = nd2a; ipar_p[8] = nd2b; ipar_p[9] = nd2;
	ipar_p[10] = nd; ipar_p[11] = numSteps;
	// Vector to save real parameters
	rpar_p[0] = rx; rpar_p[1] = ry; rpar_p[2] = kappa; rpar_p[3] = dx[0]; rpar_p[4] = dx[1];
	rpar_p[5] = dt; rpar_p[6] = c0; rpar_p[7] = c1; rpar_p[8] = c2; rpar_p[9] = b0; rpar_p[10] = b1;
	rpar_p[11] = b2; rpar_p[12] = a0; rpar_p[13] = a1; rpar_p[14] = a2;
	*/

	// update U0 and U1 to get started
	Real t = 0.;
	int prev = 0, cur = 1;
	for (i2=nd2a; i2<=nd2b; i2++)
	{
		for (i1=nd1a; i1<=nd1b; i1++)
		{
			up(i1,i2)= UTRUE(x(i1,i2),y(i1,i2),t);
		}
	}
	// Approximation for U1
	// Use Taytor expansion and the PDE for interior points
	for (i2=n2a+1; i2<=n2b-1; i2++)
	{
		for (i1=n1a+1; i1<=n1b-1; i1++)
		{
			uc(i1,i2) = up(i1, i2) + dt * UTRUET(x(i1, i2), y(i1, i2), t) + (dt * dt / (2.0 * rho(i1, i2))) *
			    ((0.5 * (m11(i1 + 1, i2) + m11(i1, i2)) * up(i1 + 1, i2) -
					  0.5 * (m11(i1 + 1, i2) + 2.0 * m11(i1, i2) + m11(i1 - 1, i2)) * up(i1, i2) +
					  0.5 * (m11(i1 - 1, i2) + m11(i1, i2)) * up(i1 - 1, i2) ) / (dx[0] * dx[0])
			+    (0.5 * (m12(i1 + 1, i2) + m12(i1, i2)) * (up(i1 + 1, i2 + 1) - up(i1 + 1, i2 - 1)) - 
					  0.5 * (m12(i1 - 1, i2) + m12(i1, i2)) * (up(i1 - 1, i2 + 1) - up(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1])
			+    (0.5 * (m21(i1, i2 + 1) + m21(i1, i2)) * (up(i1 + 1, i2 + 1) - up(i1 - 1, i2 + 1)) -
					  0.5 * (m21(i1, i2 - 1) + m21(i1, i2)) * (up(i1 + 1, i2 - 1) - up(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1])
			+    (0.5 * (m22(i1, i2 + 1) + m22(i1, i2)) * up(i1, i2 + 1) -
					  0.5 * (m22(i1, i2 + 1) + 2.0 * m22(i1, i2) + m22(i1, i2 - 1)) * up(i1, i2) +
					  0.5 * (m22(i1, i2 - 1) + m22(i1, i2)) * up(i1, i2 - 1) ) / (dx[1] * dx[1]) + FORCE(x(i1, i2), y(i1, i2), t) );
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
				Real ue = UTRUE(x(i1,i2), y(i1,i2), dt);
				err = fabs(uc(i1,i2) - UTRUE(x(i1,i2), y(i1,i2), dt));
				maxErr = max(err, maxErr);
				// printf("err[%d,%d]=%6.2e, uc=%4.7f, ue=%4.7f\n", i1, i2, err, uc(i1,i2), ue);
			}
		}
		printf("Debug: check error for the approximation of U1: maxErr=%6.2e", maxErr);
		return 0;
	}
	
	printf("----- Solving Anisotropic Wave Equation in two dimensions ------\n");
	printf(" saveMatlab=%d, matlabFileName=%s \n",saveMatlab,matlabFileName.c_str());
	printf(" nx=%d ny=%d dt=%6.2e numSteps=%d tFinal=%6.2f\n", nx, ny,dt, numSteps, tFinal);

	// ---------------------------- TIME-STEPPING LOOP ------------------------------
	Real cpu1 = getCPU();
	for(int n=0; n<numSteps-1; n++)
	{
		t = (n+1)*dt;  // current time
		const int prev = n%3;         // previous time level
		const int cur  = (n+1)%3;     // current time level
		const int next = (n+2)%3;     // next time level
		// printf("perv=%d cur=%d next=%d\n", prev, cur, next);
		// Update the interior points
		for(i2=n2a+1; i2<=n2b-1; i2++)
		{
			for(i1=n1a+1; i1<=n1b-1; i1++)
			{
				un(i1,i2) = 2 * uc(i1, i2) - up(i1, i2) + (dt * dt / rho(i1, i2))*
					((0.5 * (m11(i1 + 1, i2) + m11(i1, i2)) * uc(i1 + 1, i2) -
					  0.5 * (m11(i1 + 1, i2) + 2.0 * m11(i1, i2) + m11(i1 - 1, i2)) * uc(i1, i2) +
					  0.5 * (m11(i1 - 1, i2) + m11(i1, i2)) * uc(i1 - 1, i2) ) / (dx[0] * dx[0])
			+    (0.5 * (m12(i1 + 1, i2) + m12(i1, i2)) * (uc(i1 + 1, i2 + 1) - uc(i1 + 1, i2 - 1)) - 
					  0.5 * (m12(i1 - 1, i2) + m12(i1, i2)) * (uc(i1 - 1, i2 + 1) - uc(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1])
			+    (0.5 * (m21(i1, i2 + 1) + m21(i1, i2)) * (uc(i1 + 1, i2 + 1) - uc(i1 - 1, i2 + 1)) -
					  0.5 * (m21(i1, i2 - 1) + m21(i1, i2)) * (uc(i1 + 1, i2 - 1) - uc(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1])
			+    (0.5 * (m22(i1, i2 + 1) + m22(i1, i2)) * uc(i1, i2 + 1) -
					  0.5 * (m22(i1, i2 + 1) + 2.0 * m22(i1, i2) + m22(i1, i2 - 1)) * uc(i1, i2) +
					  0.5 * (m22(i1, i2 - 1) + m22(i1, i2)) * uc(i1, i2 - 1) ) / (dx[1] * dx[1]) + FORCE(x(i1, i2), y(i1, i2), t) );
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

	// --- compute errors ---
	// --------- Check the error ---------
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
	maxErr = maxErr/maxNorm;
	printf("CPU: numSteps=%d nx=%d maxNorm=%8.2e maxRelErr=%8.2e cpuTime=%9.2e(s)\n",
					numSteps, nx, maxNorm, maxErr, cpuTimeStep);

	

	// Free allocated vectors
	// delete [] gridIndex_p;
	freeVectors(x_p[0], x_p[1], u_p[0], u_p[1], u_p[2], errCPU_p,
							rho_p, m11_p, m12_p, m21_p, m22_p);
	
	return 0;
}


// Definition of functions declared above
void allocateVectors( 
			double** x_p, double** y_p,
			double** up_p, double** un_p,
			double** uc_p, double** errCPU,
			double** rho_p, double** m11_p,
			double** m12_p, double** m21_p,
			double** m22_p,
			int nd){
	
	*x_p    = new double [nd];
	*y_p    = new double [nd];
	*up_p   = new double [nd];
	*un_p   = new double [nd];
	*uc_p   = new double [nd];
	*errCPU = new double [nd];
	*rho_p  = new double [nd];
	*m11_p  = new double [nd];
	*m12_p  = new double [nd];
	*m21_p  = new double [nd];
	*m22_p  = new double [nd];
	
}

// Free all allocated vectors
void freeVectors( 
			double* x_p, double* y_p,
			double* up_p, double* un_p,
			double* uc_p, double* errCPU_p,
			double* rho_p, double* m11_p,
			double* m12_p, double* m21_p,
			double* m22_p){
	
	delete [] x_p;
	delete [] y_p;
	delete [] up_p;
	delete [] un_p;
	delete [] uc_p;
	delete [] errCPU_p;
	delete [] rho_p;
	delete [] m11_p;
	delete [] m12_p;
	delete [] m21_p;
	delete [] m22_p;
}
