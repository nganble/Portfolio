#ifndef POISSON2D_COMMON_H
#define POISSON2D_COMMON_H

#include<math.h>

/* ===================================================
	Basic Types
====================================================*/
typedef double real;
typedef unsigned long long ticks;


/* ===================================================
	Core Data Structures
====================================================*/
/* Struct for domain in 2d */
typedef struct {
	real xa, xb, ya, yb;
} domain_t;


/* Grid description */
typedef struct {
    int num_ghost;
    int n1a, n1b, nd1a, nd1b, nd1;
    int n2a, n2b, nd2a, nd2b, nd2;
    int nd;
    int ny_l, n2a_l, n2b_l;
    int nd2a_l, nd2b_l;
    int nd2_l, nd_l;
    real dx, dy, h, h_sq;
} grid_t;

/* Grid description used for kernel */
typedef struct {
    int nd1, nd_l;
    int nd1a, nd1b;
    int nd2a_l, nd2b_l;
    real h_sq;
} grid_kernel_t;

/* Indexing metadata */
typedef struct {
	int nd1a, nd2a, nd1;
} idx_t;

/* Solver parameters */
typedef struct {
    real omega;
    real tol;
    int max_iter;
    int iter_count;
} solver_par_t;



/* ===================================================
	Constants and Flags
====================================================*/
#define RED 0
#define BLACK 1
#define BLOCK_MAX 65535
#define PI 3.14159265358979323846
#define SQ_ROOT_PI 1.77245385090551602729

/* ===================================================
	Macros for grid points and grid functions
====================================================*/
#define x(i1) x_p[(i1) - grid.nd1a]
#define y(i2) y_p[(i2) - grid.nd2a_l]
#define u_k(i1, i2) u_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define res(i1, i2) r_p[( (i1) - grid.nd1a) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define err(i1, i2) err_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]


/* ===================================================
	Host / Device Compatibility
====================================================*/
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

/* ===================================================
	Indexing Utilities (not used)
====================================================*/
/* Convert (i1, i2) to linear index */
static inline HD int idx(int i1, int i2, const idx_t *f){
	return (i1 - f->nd1a) + f->nd1 * (i2 - f->nd2a);
}


/* ===================================================
	Problem Configuration
====================================================*/
/* Test modes */
#define EXACT_SOLUTION 1
#define COMPUTE_SOLUTION 2

#ifndef TEST
#define TEST COMPUTE_SOLUTION
#endif

/* Solution types */
#define POLY 1
#define GAUSSIAN 2

#ifndef SOLUTION
#define SOLUTION POLY
#endif

/* Analytical solution */
#if SOLUTION == POLY
	static const real a0 = 0.5, a1 = 0.7, a2 = -0.3, a3 = -0.2;
	static const real b0 = 0.75, b1 = -0.45, b2 = -0.3, b3 = 0.25;

	static inline HD real uexact(real x, real y){
		return 	( a0 + x * ( a1 + x * ( a2 + x * a3 ) ) ) *
				( b0 + y * ( b1 + y * ( b2 + y * b3 ) ) );
	}

	static inline HD real uexact_xx(real x, real y){
		return 	(2.0 * a2 + 6.0 * x * a3) *
				( b0 + y * ( b1 + y * ( b2 + y * b3 ) ) );
	}

	static inline HD real uexact_yy(real x, real y){
		return 	( a0 + x * ( a1 + x * ( a2 + x * a3 ) ) ) *
				(2.0 * b2 + 6.0 * y * b3);
	}

#elif SOLUTION == GAUSSIAN
	static const real kx = 2.0, ky = 2.0;
	static const real x_0 = 0.3, y_0 = 0.7;

	static inline HD real uexact(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

	static inline HD real uexact_xx(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * kx * kx * dx * dx - 2 * kx) * 
				exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

	static inline HD real uexact_yy(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * ky * ky * dy * dy - 2 * ky) * 
				exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

#endif

/* Forcing and Boundary Conditions */
#if TEST == EXACT_SOLUTION
	static inline HD real forcing(real x, real y){
		return uexact_xx(x, y) + uexact_yy(x, y);
	}

	#define HAS_EXACT_SOLUTION 1

#endif

#if TEST == COMPUTE_SOLUTION
	static const real lambda = 100.0;
	static inline HD real forcing(real x, real y){
		real dxp = x + 1.0;
		real dxm = x - 1.0;
		return ( 10.0 * lambda / (SQ_ROOT_PI) ) * 
		( exp(-lambda * lambda * (dxp * dxp + y * y)) 
		- exp(-lambda * lambda * (dxm * dxm + y * y)) );
	}
	static inline real HD gxa(real y){
		return 0.0;
	}
	static inline real HD gxb(real y){
		return 0.0;
	}
	static inline real HD gya(real x){
		return 0.0;
	}
	static inline real HD gyb(real x){
		return 0.0;
	}

	#define HAS_EXACT_SOLUTION 0

#endif


#endif /*.POISSON2D_COMMON_H */