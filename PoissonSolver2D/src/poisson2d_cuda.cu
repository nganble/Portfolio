// nvcc -O3 poisson2d_cuda.cu -o poisson2d_cuda
// ./poisson2d_cuda 8 1.0 1e-10 1000
// ./poisson2d_cuda 1024 2.0 1e-10 20000
// ----------------------------------------------------------
//
//			POISSON EQUATION IN 2D (CUDA VERSION)
// 
// ----------------------------------------------------------

#include <stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include<cuda.h>
#include<cuda_runtime.h>


/* New type 'real' which is equivalent to a 'double' */
typedef double real;

/* Struct for domain in 2d */
typedef struct {
	real xa, xb, ya, yb;
} domain_t;

/* Struct for 2d grid */
typedef struct {
	int num_ghost;
	int n1a, n1b, nd1a, nd1b, nd1;
	int n2a, n2b, nd2a, nd2b, nd2;
	int nd;
	real dx, dy, h, h_sq;
} grid_t;

/* Struct for SOR solver parameters */
typedef struct {
	real omega;
	real tol;
	int max_iter;
	int iter_count;
} solver_par_t;

#define BLOCK_MAX 65535
#define RED 0
#define BLACK 1
#define PI 3.14159265358979323846
#define SQ_ROOT_PI 1.77245385090551602729

/* Macros for grid points and grid functions */
#define x(i1) x_p[(i1) - grid.nd1a]
#define y(i2) y_p[(i2) - grid.nd2a]
#define u_k(i1, i2) u_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a )]
#define res(i1, i2) r_p[( (i1) - grid.nd1a) + grid.nd1 * ( (i2) - grid.nd2a )]
#define err(i1, i2) err_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a )]


#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

static inline HD void compute_launch_config(size_t nd, int *blocks, int *threads){
	*threads = 256;

	*blocks = (nd + (*threads) - 1) / (*threads);

	if (*blocks < 1) {
        *blocks = 1;
    }

    if (*blocks > BLOCK_MAX){
        *blocks = BLOCK_MAX;
    }	
}

#ifndef USE_POWER9
#define USE_POWER9 1
#endif

#ifndef USE_MPI
#define USE_MPI 0
#endif

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

/* Test solutions and options */
#define EXACT_SOLUTION 1
#define COMPUTE_SOLUTION 2

#ifndef TEST
#define TEST COMPUTE_SOLUTION
#endif

#define POLY 1
#define GAUSSIAN 2

#ifndef SOLUTION
#define SOLUTION POLY
#endif

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

	static inline real uexact(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

	static inline HD real uexact_xx(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * kx * kx * dx * dx - 2 * kx) * exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

	static inline HD real uexact_yy(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * ky * ky * dy * dy - 2 * ky) * exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

#endif

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
	static inline HD real gxa(real y){
		return 0.0;
	}
	static inline HD real gxb(real y){
		return 0.0;
	}
	static inline HD real gya(real x){
		return 0.0;
	}
	static inline HD real gyb(real x){
		return 0.0;
	}

	#define HAS_EXACT_SOLUTION 0

#endif

// #if TEST == COMPUTE_SOLUTION && defined(SOLUTION)
// 	#warning "SOLUTION is defined but not used in COMPUTE mode"
// #endif


/* ****** Functions declaration ******** */
/* ************************************* */

/* ------- Allocate vector ------------- */
void allocate_vectors(
	real **x_p			/* out */,
	real **y_p			/* out */,
	real **u_p			/* out */,
	real **r_p			/* out */,
	grid_t g 			/* in */);			

/* ------- Free vector ----------------- */
void free_vectors(
	real *x_p			/* in */,
	real *y_p			/* in */,
	real *u_p			/* in */,
	real *r_p			/* in */);

/* ------- Display the grid -------------*/
void display_grid(
	real *v 		/* in */,
	grid_t g 		/* in */);


/* ------- Get command-line args ------- */
void get_arguments(
	const int  	argc		/* in */,
	char 		*argv[]		/* in */,
	int  		*n			/* out */,
	real 		*omega		/* out */,
	real 		*tol		/* out */,
	int  		*k 			/* out */);

/* ------- Make 2d grid ---------------- */
grid_t make_grid(
	domain_t	dom 	/* in */,
	int 		nx		/* in */,
	int 		ny		/* in */);

/* Function compute maximum */
void get_max(
	real *u_p 		/* in */, 
	real *max_val 	/* out */,
	int nd 			/* in */);

/* ****** CUDA kernel declaration ****** */
/* ************************************* */

/* CUDA kernel for updating red/black points */
__global__
void red_black_step(
	real	*u_p 	/* out */,
	real	*x_p 	/* in */,
	real	*y_p 	/* in */, 
	grid_t  grid 	/* in */,
	real 	omega 	/* in */,
	int 	color 	/* in */);

/* CUDA kernel for computing the residual */
__global__
void compute_residual(
	real	*u_p 	/* in */,
	real	*r_p 	/* out */,
	real	*x_p 	/* in */,
	real	*y_p 	/* in */, 
	grid_t  grid 	/* in */);

/* CUDA kernel for max reduction */
__global__
void reduce_max(
	real *u_p 			/* in */, 
	real *block_max 	/* out */,
	int nd 				/* in */);

// -------------------- main -------------------- //
int main(int argc, char *argv[]){

	/* -------- Options to debug ------ */
	int debug = 1;
	bool display_mesh = false;
	bool display_ini = false;
	bool display_res = false;
	bool display_err = false;

	// Problem domains
	const domain_t dom ={
		.xa = 0.0, .xb = 1.0,
		.ya = 0.0, .yb = 1.0
	};

	// Default values but can get different values from the command line
	int nx = 9;

	solver_par_t solver_par;
	solver_par.iter_count = 0;

	real max_res = 0.0;

	
	// Get input arguments
	get_arguments(argc, argv, &nx, &solver_par.omega, &solver_par.tol, &solver_par.max_iter);

	// Set ny=nx for this particular problem
	nx += 1;
	int ny = nx;
	
	// Make grid
	grid_t grid = make_grid(dom, nx, ny);

	printf("**************************************************************************\n");

	// Use optimal value for omega if solver_par.omega = 2.0
	if (fabs(solver_par.omega - 2.0) < 1e-15 ){
		solver_par.omega = 2 / (1 + sin(PI * grid.h));
		printf("Use optimal value for omega.\n");
	}

	printf("-----------------------------------------------------------------------\n");
	printf("2D Poisson Solver (CUDA): nx=ny=%d\nomega=%g, tol=%.2e, max_iter=%d\n",
		nx - 1, solver_par.omega, solver_par.tol, solver_par.max_iter);
	printf("----------------------------------------------------------------------\n");

	real *x_p, *y_p, *u_p, *r_p;

	allocate_vectors(&x_p, &y_p, &u_p, &r_p, grid);


	// Fill in grid points
	for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
		x(i1) = dom.xa + (i1 - grid.n1a) * grid.dx;
	}

	for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2++){
		y(i2) = dom.ya + (i2 -grid.n2a) * grid.dy;
 	}

	if ( (debug == 1) && (display_mesh == true) ){
		printf("nd1a=%d, nd1b=%d, nd1=%d, dx=%f, xa=%.1f xb=%.1f\n",
				grid.nd1a, grid.nd1b, grid.nd1, grid.dx, dom.xa, dom.xb);
		printf("nd2a=%d, nd2b=%d, nd2=%d, dy=%f, ya=%.1f yb=%.1f\n",
				grid.nd2a, grid.nd2b, grid.nd2, grid.dy, dom.ya, dom.yb);
		printf("-------------------------- grid vector for x --------------------------\n");
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
			printf("x[%d]=%.3f\n", i1, x(i1));
		}
		printf("-------------------------- grid vector for y -------------------------- \n");
		for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2++){
			printf("y[%d]=%.3f\n", i2, y(i2));
		}
		
	}

	/* ------------ Fill in boundary points ---------------- */
	// For computing the actuall solution --> boundary conditions are zeros
	#if HAS_EXACT_SOLUTION
		// Left boundary and right boundaries
		for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2++){
			u_k(grid.n1a, i2) = uexact(dom.xa, y(i2));
			u_k(grid.n1b, i2) = uexact(dom.xb, y(i2));
		}

		// Bottom and top boundaries
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
			u_k(i1, grid.n2a) = uexact(x(i1), dom.ya);
			u_k(i1, grid.n2b) = uexact(x(i1), dom.yb);
		}

		if ( (debug == 1) && (display_ini == true) ){
			printf("-------------- initial guess u_k with boundaries fill in -------------- \n");
			display_grid(u_p, grid);
		}
	#endif

	// Prefetch the arrays to GPU
	cudaMemPrefetchAsync(x_p, grid.nd1 * sizeof(real), 0);
	cudaMemPrefetchAsync(y_p, grid.nd2 * sizeof(real), 0);
	cudaMemPrefetchAsync(u_p, grid.nd * sizeof(real), 0);
	cudaMemPrefetchAsync(r_p, grid.nd * sizeof(real), 0);
	cudaDeviceSynchronize(); // Ensure prefetch finished

	int threads, blocks;
	compute_launch_config(grid.nd, &blocks, &threads);
	
	// =========================== Launch kernel for iterative solver ====================== //
	double t_start = get_time();
	for (int k = 1; k <= solver_par.max_iter; k++){

		solver_par.iter_count += 1;

		// --------------- Update red points ----------------- //
		red_black_step<<<blocks, threads>>>(u_p, x_p, y_p, grid, solver_par.omega, RED);

		// if (debug == 1){
		// 	printf("-------------- one step: filling in red point -------------- \n");
		// 	display_grid(u_p, grid);
		// 	break;
		// }

		// ---------- Update black points ----------------- //
		red_black_step<<<blocks, threads>>>(u_p, x_p, y_p, grid, solver_par.omega, BLACK);

		// if (debug == 1){
		// 	printf("-------------- one step: filling in black point -------------- \n");
		// 	display_grid(u_p, grid);
		// 	break;
		// }

		// Compute the residual after few steps (e.g. 5 or 10 iterations)
		if ( (k % 1) == 0 ){

			/* Residual before update */
			// display_grid(r_p, grid);

			compute_residual<<<blocks, threads>>>(u_p, r_p, x_p, y_p, grid);

			/* Check whether the kernal launched successfully */
			// cudaError_t err = cudaGetLastError();
			// if (err != cudaSuccess){
			// 	printf("Launch error: %s\n", cudaGetErrorString(err));
			// 	break;
			// }

			/* Compute maximum residual (UNCOMMENT AFTER DEBUGGING) */
			max_res = 0;
			get_max(r_p, &max_res, grid.nd);

			
			/* Recompute the residual for checking (serial code) COMMMENT LATER */
			// max_res = 0;
			// for (int i2 = nd2a + 1; i2 <= nd2b - 1; i2++){
			// 	for (int i1 = nd1a + 1; i1 <= nd1b - 1; i1++){
			// 		res(i1, i2) = u_k(i1 + 1, i2)  + u_k(i1 - 1, i2) - 4 * u_k(i1, i2) +
			// 		u_k(i1, i2 + 1) + u_k(i1, i2 - 1) - h_sq * forcing(x(i1, i2), y(i1, i2));

			// 		max_res = fmax( max_res, fabs(res(i1, i2)) );
			// 	}
			// }


			if (max_res < solver_par.tol){
				printf("SOR converged after %d iterations!\nMax residual: %.2e, ", k, max_res);
				/* Residual after update */
				// display_grid(r_p, grid);
				break;
			}
		}
		
	} /* End SOR loop */

	double t_end = get_time();

	if ( (debug == 1) && (display_res == true) ){
		printf("-------------- SOR with w=%.2f: Residual -------------- \n", solver_par.omega);
		display_grid(r_p, grid);
		printf("Max residual: %.e\n", max_res);
		
	}

	if (solver_par.iter_count == solver_par.max_iter){
		printf("SOR failed to converge after %d iterations! Max residual: %.2e, ", solver_par.iter_count, max_res);
	}

	double run_time = t_end - t_start;
	printf("solver time: %.3f(s)\n", run_time);
	printf("Time per iter: %.4f(ms)\n", 1000.0 * run_time / solver_par.iter_count);
	
	/* ------ Compare with exact solution (if exists) ------- */
	#if HAS_EXACT_SOLUTION
		real max_err = 0;

		real *err_p;
		err_p = (real*) calloc(grid.nd, sizeof(*err_p)); // jsut call on CPU

		for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2++){
			for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
				err(i1, i2) = fabs(uexact(x(i1), y(i2)) - u_k(i1, i2));
				max_err = fmax( max_err, err(i1, i2) );
			}
		}
		if ( (debug == 1) && (display_err == true) ){
			printf("-------------------- Maximum Error -------------------- \n");
			display_grid(err_p, grid);
		}
		printf("Max error: %.2e\n", max_err);

		free(err_p);
	#endif

	printf("\n");

	// Free all allocated vectors
	free_vectors(x_p, y_p, u_p, r_p);

} /* main */

/* ******************** CPU FUNCTIONS ******************** */
/* ******************************************************* */

void allocate_vectors(
	real **x_p			/* out */,
	real **y_p			/* out */,
	real **u_p			/* out */,
	real **r_p			/* out */,
	grid_t grid			/* in */){

	/* Allocate Unified Memory -- accessible from CPU or GPU */
	cudaMallocManaged(x_p, grid.nd1 * sizeof(real));
	cudaMallocManaged(y_p, grid.nd2 * sizeof(real));
	cudaMallocManaged(u_p, grid.nd * sizeof(real));
	cudaMallocManaged(r_p, grid.nd * sizeof(real));

	cudaMemset(*u_p, 0.0, grid.nd * sizeof(real));		// set all entires of u_p to be zeros
	cudaMemset(*r_p, 0.0, grid.nd * sizeof(real));		// set all entires of r_p to be zeros

}	/* allocate_vectors */

// ---------------------------------------------------------------------------- //
void free_vectors(
	real *x_p			/* in */,
	real *y_p			/* in */,
	real *u_p			/* in */,
	real *r_p			/* in */){

	cudaFree(x_p);
	cudaFree(y_p);
	cudaFree(u_p);
	cudaFree(r_p);

}	/* free_vectors */

// ---------------------------------------------------------------------------- //
void display_grid(
	real *v 		/* in */,
	grid_t grid 	/* in */){

	for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2 ++){
		printf("v^T[:,%d]= [", i2);
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1 ++){
			printf("%.g ", v[i1 - grid.nd1a + grid.nd1 * (i2 - grid.nd2a)]);
		}
		printf("]\n");

	}

}	/* display_grid_function */

// ---------------------------------------------------------------------------- //
void get_arguments(
	const int  	argc		/* in */,
	char 		*argv[]		/* in */,
	int  		*n			/* out */,
	real 		*omega		/* out */,
	real 		*tol		/* out */,
	int  		*k 			/* out */){

	if (argc != 5){
		printf("Usage: ./poisson2d n omega tol k!\n");
		exit(EXIT_FAILURE);
	}

	*n 		= strtol(argv[1], NULL, 10);
	*omega 	= strtod(argv[2], NULL);
	*tol 	= strtod(argv[3], NULL);
	*k		= strtol(argv[4], NULL, 10);

}	/* get_arguments */

/* ------- Make 2d grid ------- */
grid_t make_grid(
	domain_t	dom 	/* in */,
	int 		nx		/* in */,
	int 		ny		/* in */){

	grid_t grid;

	grid.num_ghost = 0; // Assume Dichlet boundary conditions

	grid.n1a 	= 0;
	grid.n1b 	= grid.n1a + nx;
	grid.nd1a 	= grid.n1a - grid.num_ghost;
	grid.nd1b   = grid.n1b + grid.num_ghost;
	grid.nd1 	= grid.nd1b - grid.nd1a + 1;

	grid.n2a 	= 0;
	grid.n2b 	= grid.n2a + ny;
	grid.nd2a 	= grid.n2a - grid.num_ghost;
	grid.nd2b 	= grid.n2b + grid.num_ghost;
	grid.nd2 	= grid.nd2b - grid.nd2a + 1;

	grid.nd = grid.nd1 * grid.nd2;
	grid.dx = (dom.xb - dom.xa) / nx;
	grid.dy = (dom.yb - dom.ya) / ny;
	grid.h  = grid.dx;
	grid.h_sq  = grid.h * grid.h;

	return grid;
} /* make_grid */


/* ******************** CUDA UTILITY ********************* */
/* ******************************************************* */
// ---------------------------------------------------------------------------- //
/* Compute the max by launching the kernel several times */
void get_max(real *u_p, real *max_val, int nd){

	real *in = u_p;
	int n = nd;
	real *out;
	while (n > 1){
		int threads, blocks;
		compute_launch_config(n, &blocks, &threads);

		size_t shared_mem = threads * sizeof(real);

		cudaMallocManaged(&out, blocks * sizeof(real));

		reduce_max<<<blocks, threads, shared_mem>>>(in, out, n);

		cudaDeviceSynchronize();

		if( in != u_p){
			cudaFree(in);
		}

		in = out;
		n = blocks;

	}

	*max_val = in[0];

	cudaFree(in);
	
} /* get_max */


/* ******************** CUDA KERNELS ********************* */
/* ******************************************************* */

// ---------------------------------------------------------------------------- //
__global__
void red_black_step(real *u_p, real	*x_p, real *y_p, grid_t grid, real omega, int color){

	// Unpack necsserary values
	int nd1   = grid.nd1;
    int nd 	  = grid.nd;
    int nd1a  = grid.nd1a;
    int nd1b  = grid.nd1b;
    int nd2a  = grid.nd2a;
    int nd2b  = grid.nd2b;
    real h_sq = grid.h_sq;

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	// Total number of threads in a grid
	int stride = blockDim.x * gridDim.x;
	// Using grid striding
	for (int idx = index; idx < nd; idx += stride){

		int global_idx = idx;
		int i2 = global_idx / nd1 + nd2a;
		int i1 = global_idx - nd1 * (i2 - nd2a) + nd1a;

		if ( (nd1a + 1 <= i1 &&  i1 <= nd1b - 1) &&
			 (nd2a + 1 <= i2 &&  i2 <= nd2b - 1)){

			// Updating the point matching the color
			if ( (i1 + i2) % 2 == color){

				real fval = forcing(x(i1), y(i2));

				real uval = u_k(i1 + 1, i2)  + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1);

				// This way may be more efficient
				u_k(i1, i2) = (1 - omega) * u_k(i1, i2) +
					(omega / 4) * (uval - h_sq * fval);			

			}		

		}

	}

} /* red_black_step */

// ---------------------------------------------------------------------------- //
/* CUDA kernel for computing the residual */
__global__
void compute_residual(real	*u_p, real *r_p, real *x_p, real *y_p, grid_t grid){

	// Unpack necsserary values
	int nd1   = grid.nd1;
    int nd 	  = grid.nd;
    int nd1a  = grid.nd1a;
    int nd1b  = grid.nd1b;
    int nd2a  = grid.nd2a;
    int nd2b  = grid.nd2b;
    real h_sq = grid.h_sq;

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	// Total number of threads in a grid
	int stride = blockDim.x * gridDim.x;
	// Using grid striding
	for (int idx = index; idx < nd; idx += stride){

		int global_idx = idx;
		int i2 = global_idx / nd1 + nd2a;
		int i1 = global_idx - nd1 * (i2 - nd2a) + nd1a;

		if ( (nd1a + 1 <= i1 &&  i1 <= nd1b - 1) &&
			 (nd2a + 1 <= i2 &&  i2 <= nd2b - 1)){

			real fval = forcing(x(i1), y(i2));

			real uval = u_k(i1 + 1, i2)  + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1);

			res(i1, i2) = (uval - 4 * u_k(i1, i2))  - h_sq * fval;
			res(i1, i2) = res(i1, i2) / grid.h;
			// printf("res(%d,%d)=%.3f\n", i1, i2, res(i1, i2));	

		}
	}
}

// ---------------------------------------------------------------------------- //
/* CUDA kernel for max reduction */
__global__
void reduce_max(real *u, real *block_max, int nd){
	extern __shared__ real sdata[];

	int thread_idx = threadIdx.x;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;

	// Compute the thread local max reduction (grid striding)
	real local_max = 0.0;
	for (int i = index; i < nd; i += stride){
		local_max = fmax(local_max, fabs(u[i]));
	}

	// Load the local result into shared memory
	sdata[thread_idx] = local_max;
	__syncthreads();

	// Reduction in each block using shared memory (always assume blockDim.x is a power of 2)
	for (int s = blockDim.x / 2; s > 0; s >>= 1){
		if (thread_idx < s){
			sdata[thread_idx] = fmax(sdata[thread_idx], sdata[thread_idx + s]);
		}
		__syncthreads();
	}

	// Write result for this block
	if (thread_idx == 0){
		block_max[blockIdx.x] = sdata[0];
	}
}


