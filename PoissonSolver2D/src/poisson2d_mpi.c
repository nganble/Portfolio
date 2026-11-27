// mpicc -O3 -o -O3 poisson2d_mpi.c -o poisson2d_mpi
// mpirun -np 2 ./poisson2d_mpi 8 1.0 1e-10 20000
// mpirun -np 2 ./poisson2d_mpi 1024 2.0 1e-10 20000
// ----------------------------------------------------------
//
// 			POISSON EQUATION IN 2D (MPI VERSION)
// 
// ----------------------------------------------------------

#include <stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include<mpi.h>


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
	int ny_l, n2a_l, n2b_l; 	// added here
	int nd2a_l, nd2b_l;			// added here
	int nd2_l, nd_l;			// added here
	real dx, dy, h, h_sq;
} grid_t;

/* Struct for metadata */
typedef struct {
	int nd1a, nd2a, nd1;
} idx_t;

/* Struct for SOR solver parameters */
typedef struct {
	real omega;
	real tol;
	int max_iter;
	int iter_count;
} solver_par_t;


/* Macros for grid points and grid functions */
#define x(i1) x_p[(i1) - grid.nd1a]
#define y(i2) y_p[(i2) - grid.nd2a_l]
#define u_k(i1, i2) u_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define res(i1, i2) r_p[( (i1) - grid.nd1a) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define err(i1, i2) err_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]

#ifndef USE_POWER9
#define USE_POWER9 1
#endif

#ifndef USE_MPI
#define USE_MPI 1
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


/* Macro for printing debug info */
#define DEBUG_PRINT(flag, mgs) 	\
	if (debug && (flag)){ 	\
		debug_print_ranked(my_rank, np, mgs, u_p, grid, ny, solver_par); \
		goto cleanup; \
}

// May be better to use inline function
static inline int idx(int i1, int i2, const idx_t *f){
	return (i1 - f->nd1a) + f->nd1 * (i2 - f->nd2a);
}

#define PI 3.14159265358979323846
#define SQ_ROOT_PI 1.77245385090551602729


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

	static inline real uexact(real x, real y){
		return 	( a0 + x * ( a1 + x * ( a2 + x * a3 ) ) ) *
				( b0 + y * ( b1 + y * ( b2 + y * b3 ) ) );
	}

	static inline real uexact_xx(real x, real y){
		return 	(2.0 * a2 + 6.0 * x * a3) *
				( b0 + y * ( b1 + y * ( b2 + y * b3 ) ) );
	}

	static inline real uexact_yy(real x, real y){
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

	static inline real uexact_xx(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * kx * kx * dx * dx - 2 * kx) * exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

	static inline real uexact_yy(real x, real y){
		real dx = x - x_0;
		real dy = y - y_0;
		return (4 * ky * ky * dy * dy - 2 * ky) * exp(-kx * dx * dx) * exp(-ky * dy * dy);
	}

#endif

#if TEST == EXACT_SOLUTION
	static inline real forcing(real x, real y){
		return uexact_xx(x, y) + uexact_yy(x, y);
	}

	#define HAS_EXACT_SOLUTION 1

#endif

#if TEST == COMPUTE_SOLUTION
	static const real lambda = 100.0;
	static inline real forcing(real x, real y){
		real dxp = x + 1.0;
		real dxm = x - 1.0;
		return ( 10.0 * lambda / (SQ_ROOT_PI) ) * 
		( exp(-lambda * lambda * (dxp * dxp + y * y)) 
		- exp(-lambda * lambda * (dxm * dxm + y * y)) );
	}
	static inline real gxa(real y){
		return 0.0;
	}
	static inline real gxb(real y){
		return 0.0;
	}
	static inline real gya(real x){
		return 0.0;
	}
	static inline real gyb(real x){
		return 0.0;
	}

	#define HAS_EXACT_SOLUTION 0

#endif

// #if TEST == COMPUTE_SOLUTION && defined(SOLUTION)
// 	#warning "SOLUTION is defined but not used in COMPUTE mode"
// #endif

/* ****** Functions declaration ******** */
/* ************************************* */

/* ------------- Allocate vector ------------- */
void allocate_vectors(
	real **x_p			/* out */,
	real **y_p			/* out */,
	real **u_p			/* out */,
	real **r_p			/* out */,
	grid_t g 			/* in */);			

/* ------------- Free vector ----------------- */
void free_vectors(
	real *x_p			/* in */,
	real *y_p			/* in */,
	real *u_p			/* in */,
	real *r_p			/* in */);

/* --------- Display the grid ---------------- */
void display_grid(
	real *v 		/* in */,
	grid_t g 		/* in */);

/* ---------- Get command-line args ---------- */
void get_arguments(
	const int  	argc		/* in */,
	char 		*argv[]		/* in */,
	int 		my_rank		/* in */,
	int  		*n			/* out */,
	real 		*omega		/* out */,
	real 		*tol		/* out */,
	int  		*k 			/* out */);

/* ----------- Make 2d grid ------------------ */
grid_t make_grid(
	domain_t	dom 	/* in */,
	int 		nx		/* in */,
	int 		ny		/* in */);

/* ------- Get local bound index  ------------ */
void get_local_index_bounds(
	const int 	np			/* in */,
	const int 	my_rank		/* in */,
	const int 	ni			/* in */,
	int 		*ni_l		/* out */,
	int 		*nia_l		/* out */,
	int 		*nib_l		/* out */);

/* ----------- Print debug info  ------------- */
void debug_print_ranked(
	int 			my_rank 	/* in */,
	int 			np 			/* in */,
	const char 		*msg 		/* in */,
	real 			*u_p 		/* in */,
	grid_t 			grid 		/* in */,
	int 			ny 			/* in */,
	solver_par_t 	solver_par 	/* in */);

/* ------------- Main function --------------- */ 
int main(int argc, char *argv[]){

	int np;
	int my_rank;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int debug = 1;

	// Problem domains
	const domain_t dom ={
		.xa = 0.0, .xb = 1.0,
		.ya = 0.0, .yb = 1.0
	};
	
	// Default values but can get different values from the command line
	int nx;

	solver_par_t solver_par;
	solver_par.iter_count = 0;

	real max_res = 0.0;
	real max_res_g = 0.0;
	
	// Get input arguments
	get_arguments(argc, argv, my_rank, &nx, &solver_par.omega, &solver_par.tol, &solver_par.max_iter);

	// Set ny=nx for this particular problem
	nx += 1; 
	int ny = nx;

	// Make grid
	grid_t grid = make_grid(dom, nx, ny);

	if (my_rank == 0){
		printf("**************************************************************************\n");
	}

	// Use optimal value for omega if solver_par.omega = 2.0
	if (fabs(solver_par.omega - 2.0) < 1e-15 ){
		solver_par.omega = 2 / (1 + sin(PI * grid.h));
		if(my_rank == 0){
			printf("Use optimal value for omega.\n");
		}
	}

	if (my_rank == 0){
		printf("-----------------------------------------------------------------------\n");
		printf("2D Poisson Solver (MPI): nx=ny=%d, np=%d\nomega=%g, tol=%.2e, max_iter=%d\n",
			nx - 1, np, solver_par.omega, solver_par.tol, solver_par.max_iter);
		printf("----------------------------------------------------------------------\n");
	}
	
	// We will devide in the y direction 
	// int ny_l, n2a_l, n2b_l;
	get_local_index_bounds(np, my_rank, ny, &grid.ny_l, &grid.n2a_l, &grid.n2b_l);
	// Recompute local values (NEED TO MAKE IT CLEAN LATTER)
	int num_ghost_top =0, num_ghost_bottom = 0;

	if (my_rank > 0){
		num_ghost_bottom += 1;
	}
	if (my_rank < np - 1){
		num_ghost_top += 1;
	}
	grid.nd2a_l = grid.n2a_l - num_ghost_bottom;
	grid.nd2b_l = grid.n2b_l + num_ghost_top;
	grid.nd2_l 	= grid.nd2b_l - grid.nd2a_l + 1;
	grid.nd_l 	= grid.nd1 * grid.nd2_l;

	// -------- Options to debug ----------
	bool display_ginfo 		= false;
	bool display_mesh 		= false;
	bool display_ini 		= false;
	bool display_one_step_r	= false;
	bool display_one_step_b = false;
	bool display_comm 		= false;
	bool display_res 		= false;
	bool display_err 		= false;

	if (debug == 1 && display_ginfo == true){
		printf("Rank %d out of %d ranks: ny=%d, ny_l=%d, n2a_l=%d, n2b_l=%d, nd2a_l=%d, nd2b_l=%d, omega=%.2f\n",
			my_rank, np, ny, grid.ny_l, grid.n2a_l, grid.n2b_l, grid.nd2a_l, grid.nd2b_l, solver_par.omega);
		MPI_Finalize();
		return 0;
	}

	// Declare grid vectors
	real *x_p, *y_p, *u_p, *r_p;

	allocate_vectors(&x_p, &y_p, &u_p, &r_p, grid);

	/* ------------ Fill in grid points ---------------- */
	for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
		x(i1) = dom.xa + (i1 - grid.n1a) * grid.dx;
	}

	for (int i2 = grid.nd2a_l; i2 <= grid.nd2b_l; i2++){
		y(i2) = dom.ya + (i2 - grid.n2a) * grid.dy;
 	}

	if ( (my_rank==0) && (debug == 1) && (display_mesh == true) ){
		printf("Rank %d out of %d ranks:\n", my_rank, np);
		printf("nd_l=%d, nd2_l=%d, n2a_l=%d, n2b_l=%d, dx=%f\n",
				grid.nd_l, grid.nd2_l, grid.n2a_l, grid.n2b_l, grid.dx);
		// printf("nd2a=%d, nd2b=%d, nd2=%d, dy=%f, ya=%.1f yb=%.1f\n",
		// 		grid.nd2a, grid.nd2b, grid.nd2, grid.dy, dom.ya, dom.yb);

		printf("-------------------------- grid vector for x --------------------------\n");
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
			printf("x[%d]=%.3f\n", i1, x(i1));
		}

		printf("-------------------------- grid vector for y -------------------------- \n");
		for (int i2 = grid.nd2a_l; i2 <= grid.nd2b_l; i2++){
			printf("y[%d]=%.3f\n", i2, y(i2));
		}
		
	}

	/* ------------ Fill in boundary points ---------------- */
	// For computing the actual solution --> boundary conditions are zeros
	#if HAS_EXACT_SOLUTION
		// Left boundary and right boundaries
		for (int i2 = grid.nd2a_l; i2 <= grid.nd2b_l; i2++){
			u_k(grid.n1a, i2) = uexact(dom.xa, y(i2));
			u_k(grid.n1b, i2) = uexact(dom.xb, y(i2));
		}

		// Bottom and top boundaries	
		if (my_rank == 0){
			for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
				u_k(i1, grid.n2a) = uexact(x(i1), dom.ya);
			}
		}

		if (my_rank == np - 1){
			for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
				u_k(i1, grid.n2b) = uexact(x(i1), dom.yb);
			}
		}
		

		if ( (debug == 1) && (display_ini == true) ){
			printf("-------------- initial guess u_k with boundaries fill in -------------- \n");
			display_grid(u_p, grid);
			goto cleanup;
		}
	#endif

	double t_start = get_time();
	/* ------------------ SOR iterative solver --------------------- */
	for (int k = 1; k <= solver_par.max_iter; k++){

		solver_par.iter_count += 1;

		
		// Update red points first (UNCOMMENT AFTER DEBUGGING)
		for (int i2 = grid.nd2a_l + 1; i2 <= grid.nd2b_l - 1; i2++){
			for (int i1 = (i2 + 1) % 2 + 1; i1 <= grid.nd1b - 1; i1 += 2){
				u_k(i1, i2) = (1 - solver_par.omega) * u_k(i1, i2) +
				(solver_par.omega / 4) * ( u_k(i1 + 1, i2) + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1) -
				grid.h_sq * forcing(x(i1), y(i2)) );		
			}
		}

		// Check if red points are updated correctly
		DEBUG_PRINT(display_one_step_r, "-------------- one step: filling in red points --------------\n");

		// Communication between processes
		MPI_Request requests[4];
		int req_count = 0;
		if (my_rank != 0 ){
			// From my current rank post a receive to get data from my_rank - 1
			MPI_Irecv(&u_k(grid.nd1a, grid.nd2a_l), grid.nd1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
		}

		if (my_rank != np - 1){
			// From my current rank post a receive to get data from my_rank + 1
			MPI_Irecv(&u_k(grid.nd1a, grid.nd2b_l), grid.nd1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &requests[req_count++]);

		}

		if (my_rank != np - 1){
			// From my_rank send to my_rank + 1
			MPI_Isend(&u_k(grid.nd1a, grid.n2b_l), grid.nd1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
		}

		if (my_rank != 0){
			// From my_rank send to my_rank - 1
			MPI_Isend(&u_k(grid.nd1a, grid.n2a_l), grid.nd1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
		}

		// Wait for the communication to finish
		MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

		// Check the communication
		DEBUG_PRINT(display_comm, "-------------- after communication --------------\n");


		// Update black points (UNCOMMENT AFTER DEBUGGING)
		for (int i2 = grid.nd2a_l + 1; i2 <= grid.nd2b_l - 1; i2++){
			for (int i1 = i2 % 2 + 1; i1 <= grid.nd1b - 1; i1 += 2){				
				u_k(i1, i2) = (1 - solver_par.omega) * u_k(i1, i2) +
				(solver_par.omega / 4) * ( u_k(i1 + 1, i2) + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1) -
				grid.h_sq * forcing(x(i1), y(i2)) );				
			}
		}

		// Check if blacks points are updated correctly
		DEBUG_PRINT(display_one_step_b, "-------------- one step: filling in black points --------------\n");
		

		
		// Communication after black points got updated
		// MPI_Request requests[4];
		// int req_count = 0;
		req_count = 0;
		if (my_rank != 0 ){
			// From my current rank post a receive to get data from my_rank - 1
			MPI_Irecv(&u_k(grid.nd1a, grid.nd2a_l), grid.nd1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
		}

		if (my_rank != np - 1){
			// From my current rank post a receive to get data from my_rank + 1
			MPI_Irecv(&u_k(grid.nd1a, grid.nd2b_l), grid.nd1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &requests[req_count++]);

		}

		if (my_rank != np - 1){
			// From my_rank send to my_rank + 1
			MPI_Isend(&u_k(grid.nd1a, grid.n2b_l), grid.nd1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
		}

		if (my_rank != 0){
			// From my_rank send to my_rank - 1
			MPI_Isend(&u_k(grid.nd1a, grid.n2a_l), grid.nd1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
		}

		// Wait for the communication to finish
		MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

		
		// Compute the residual
		if (k % 1 == 0){
			max_res = 0;
			for (int i2 = grid.nd2a_l + 1; i2 <= grid.nd2b_l - 1; i2++){
				for (int i1 = grid.nd1a + 1; i1 <= grid.nd1b - 1; i1++){
					res(i1, i2) = ( u_k(i1 + 1, i2)  + u_k(i1 - 1, i2) - 4 * u_k(i1, i2) +
					u_k(i1, i2 + 1) + u_k(i1, i2 - 1) ) - grid.h_sq * forcing(x(i1), y(i2));
					res(i1, i2) = res(i1, i2) / grid.h;
					max_res = fmax( max_res, fabs(res(i1, i2)) );
				}
			}

			// Gather maximum residual in all ranks
			max_res_g = max_res;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allreduce(&max_res, &max_res_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			if (max_res_g < solver_par.tol){
				if (my_rank==0){
					printf("SOR converged after %d iterations!\nMax residual: %.2e, ", k, max_res_g);
				}
				break;
			}
		}


	} /* End SOR loop */
	double t_end = get_time();
	

	if ( (my_rank==0) && (debug == 1) && (display_res == true) ){
		printf("-------------- SOR with w=%.2f: Residual -------------- \n", solver_par.omega);
		display_grid(r_p, grid);
		printf("Max residual: %.e\n", max_res);
	}

	if (solver_par.iter_count == solver_par.max_iter && my_rank == 0){
		printf("SOR failed to converge after %d iterations! Max residual: %.2e, ", solver_par.iter_count, max_res_g);
	}

	double run_time = t_end - t_start;
	if (my_rank == 0){
		printf("solver time: %.4f(s)\n", run_time);
		printf("Time per iter: %.4f(ms)\n", 1000.0 * run_time / solver_par.iter_count);
	}

	#if HAS_EXACT_SOLUTION
		real max_err = 0.0;
		real max_err_g = 0.0;

		real *err_p;
		err_p = calloc(grid.nd_l, sizeof(*err_p));

		for (int i2 = grid.nd2a_l; i2 <= grid.nd2b_l; i2++){
			for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1++){
				err(i1, i2) = fabs(uexact(x(i1), y(i2)) - u_k(i1, i2));
				max_err = fmax( max_err, err(i1, i2) );
			}
		}
		MPI_Allreduce(&max_err, &max_err_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (my_rank == 0){
			if ( (debug == 1) && (display_err == true) ){
				printf("-------------------- Maximum Error -------------------- \n");
				display_grid(err_p, grid);
			}
			printf("Max error: %.2e\n", max_err_g);
		}

		free(err_p);

	#endif
	
	if (my_rank == 0){
		printf("\n");
	}

	cleanup:
		// Free all allocated vectors
		free_vectors(x_p, y_p, u_p, r_p);
		MPI_Finalize();
		return 0;

} /* main */


/* ****** Functions definition ********* */
/* ************************************* */

void allocate_vectors(
	real 	**x_p			/* out */,
	real 	**y_p			/* out */,
	real 	**u_p			/* out */,
	real 	**r_p			/* out */,
	grid_t 	grid 			/*in*/){

	*x_p = malloc(grid.nd1 * sizeof(real));
	*y_p = malloc(grid.nd2_l * sizeof(real));
	*u_p = calloc(grid.nd_l, sizeof(real));		// set all entires of u_p to be zeros
	*r_p = calloc(grid.nd_l, sizeof(real));		// set all entires of r_p to be zeros

}	/* allocate_vectors */


void free_vectors(
	real *x_p			/* in */,
	real *y_p			/* in */,
	real *u_p			/* in */,
	real *r_p			/* in */){

	free(x_p);
	free(y_p);
	free(u_p);
	free(r_p);

}	/* free_vectors */

void display_grid(
	real *v 		/* in */,
	grid_t grid 	/* in */){

	for (int i2 = grid.nd2a_l; i2 <= grid.nd2b_l; i2 ++){
		printf("v^T[:,%d]= [", i2);
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1 ++){
			printf("%.3f ", v[i1 - grid.nd1a + grid.nd1 * (i2 - grid.nd2a_l)]);
		}
		printf("]\n");

	}

}	/* display_grid_function */

void get_arguments(
	const int  	argc		/* in */,
	char 		*argv[]		/* in */,
	int 		my_rank		/* in */,
	int  		*n			/* out */,
	real 		*omega		/* out */,
	real 		*tol		/* out */,
	int  		*k 			/* out */){

	if (argc != 5){
		if (my_rank == 0){
			printf("Usage: ./poisson2d n omega tol k!\n");
			
		}
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		
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

void get_local_index_bounds(
	const int 	np			/* in */,
	const int 	my_rank		/* in */,
	const int 	ni			/* in */,
	int 		*ni_l		/* out */,
	int 		*nia_l		/* out */,
	int 		*nib_l		/* out */){

	*ni_l = (ni + 1) / np - 1;
	int r = (ni + 1) % np;

	// (ni+1) may not be a muliple of np, add extra points to the processors on the right side
	if (my_rank >= np - r){
		*ni_l += 1;
		*nia_l = (*ni_l) * (np - r)  + (my_rank - np + r) * ( (*ni_l) + 1);
	}
	else{
		*nia_l = my_rank * ( (*ni_l) + 1);
	}
	

	*nib_l = (*nia_l) + (*ni_l);

} /* get_local_index_bounds */

/* ----------- Print debug info  ------------- */
void debug_print_ranked(
	int 			my_rank 	/* in */,
	int 			np 			/* in */,
	const char 		*msg 		/* in */,
	real 			*u_p 		/* in */,
	grid_t 			grid 		/* in */,
	int 			ny 			/* in */,
	solver_par_t 	solver_par 	/* in */){

	for (int r = 0; r < np; r++){
		MPI_Barrier(MPI_COMM_WORLD);
		if (my_rank == r){
			printf("%s", msg);
			printf("Rank %d out of %d ranks: ny=%d, ny_l=%d,"
				" n2a_l=%d, n2b_l=%d, nd2a_l=%d, nd2b_l=%d, omega=%.2f\n",
				my_rank, np, ny, grid.ny_l, grid.n2a_l, 
				grid.n2b_l, grid.nd2a_l, grid.nd2b_l, solver_par.omega);
			display_grid(u_p, grid);
			fflush(stdout);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
} /* debug_print_ranked */

