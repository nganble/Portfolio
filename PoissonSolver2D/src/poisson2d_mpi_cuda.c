// =====================================================================
// POISSON EQUATION IN 2D: MPI+CUDA code: host part
//   Solves u_xx + u_yy = f on (0,1)^2 with Dirichlet boundary conditions,
//     5-point Laplacian, red-black SOR (omega = 2.0 -> optimal omega = 2/(1+sin(pi h))).
//   Same decomposition as poisson2d_mpi.c; each rank uses one GPU, the kernels
//   are in poisson2d_mpi_cuda.cu.
//
// Build:   make POWER9=1 mpi-cuda
// Run:     mpirun -np 4 ./poisson2d_mpi_cuda 1024 2.0 1e-10 20000
// Args:    N omega tol max_iter   (N interior points per direction, h = 1/(N+1); stop when max residual < tol)
// Cases:   -DTEST=EXACT_SOLUTION -DSOLUTION=POLY|GAUSSIAN, or -DTEST=COMPUTE_SOLUTION (default)
// =====================================================================

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mpi.h>



#include "poisson2d_common.h"

#include "poisson2d_cuda.h"

#include "poisson2d_mpi.h"

#ifndef USE_MPI
#define USE_MPI 1
#endif

#ifndef USE_POWER9
#define USE_POWER9 1
#endif

#include "timing.h"

/* Macros for grid points and grid functions */
#define x(i1) x_p[(i1) - grid.nd1a]
#define y(i2) y_p[(i2) - grid.nd2a_l]
#define u_k(i1, i2) u_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define res(i1, i2) r_p[( (i1) - grid.nd1a) + grid.nd1 * ( (i2) - grid.nd2a_l )]
#define err(i1, i2) err_p[( (i1) - grid.nd1a ) + grid.nd1 * ( (i2) - grid.nd2a_l )]


/* Macro for printing debug info */
#define DEBUG_PRINT(flag, mgs) 	\
	if (debug && (flag)){ 	\
		debug_print_ranked(my_rank, np, mgs, u_p, grid, ny, solver_par); \
		goto cleanup; \
}

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
	int nx = 9;

	solver_par_t solver_par;
	solver_par.iter_count = 0;

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
		printf("2D Poisson Solver (MPI-CUDA): nx=ny=%d np=%d\nomega=%g, tol=%.2e, max_iter=%d\n",
			nx - 1, np, solver_par.omega, solver_par.tol, solver_par.max_iter);
		printf("-----------------------------------------------------------------------\n");
	}

	
	// Split the grid in the y-direction 
	get_local_index_bounds(np, my_rank, ny, &grid.ny_l, &grid.n2a_l, &grid.n2b_l);
	// Recompute local values
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
#if HAS_EXACT_SOLUTION
	bool display_ini 		= false;
	bool display_err 		= false;
#endif
	bool display_one_step_r	= false;
	bool display_one_step_b = false;
	bool display_comm 		= false;
	bool display_res 		= false;


	if (debug == 1 && display_ginfo == true){
		printf("Rank %d out of %d ranks: ny=%d, ny_l=%d, n2a_l=%d, n2b_l=%d, nd2a_l=%d, nd2b_l=%d, omega=%.2f\n",
			my_rank, np, ny, grid.ny_l, grid.n2a_l, grid.n2b_l, grid.nd2a_l, grid.nd2b_l, solver_par.omega);
		MPI_Finalize();
		return 0;
	}

	// Declare grid vectors
	real *x_p, *y_p, *u_p, *r_p, *err_p;

	allocate_vectors(&x_p, &y_p, &u_p, &r_p, &err_p, grid);

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
			free_vectors(x_p, y_p, u_p, r_p, err_p);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return 0;
		}
	#endif

	// Get the number of GPUs
	int num_devices;
	cuda_get_device_count(&num_devices);

	// Rank r uses GPU (r mod number of GPUs on its node)
	int device = my_rank % num_devices;

	cuda_set_device(device);

	// Prefetch the arrays to GPU
	cuda_mem_prefetch(x_p, grid.nd_l, device);
	cuda_mem_prefetch(y_p, grid.nd_l, device);
	cuda_mem_prefetch(u_p, grid.nd_l, device);
	cuda_mem_prefetch(r_p, grid.nd_l, device);
	cuda_mem_prefetch(err_p, grid.nd_l, device);
	cuda_synchronize(); // Ensure prefetch finished

	real max_res;
	real max_res_g = 0.0;

	/* ------------------ SOR iterative solver --------------------- */
	double t_start = get_time();

	for (int k = 1; k <= solver_par.max_iter; k++){

		solver_par.iter_count += 1;

		
		// Launch the kernel to update RED POINTS
		cuda_red_black_step(u_p, x_p, y_p, grid, solver_par.omega, RED);

		cuda_synchronize(); // Need to wait for the job on GPU to done first

		MPI_Barrier(MPI_COMM_WORLD); // Do we need this? (Need to wait so that all ranks have done the update before we communicate)

		// Check if red points are updated correctly
		DEBUG_PRINT(display_one_step_r, "-------------- one step: filling in red points --------------\n");

		// Exchange data after red points got updated
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

		// Launch the kernel to update BLACK POINTS
		cuda_red_black_step(u_p, x_p, y_p, grid, solver_par.omega, BLACK);

		cuda_synchronize();

		MPI_Barrier(MPI_COMM_WORLD); // Do we need this?

		// Check if black points are updated correctly
		DEBUG_PRINT(display_one_step_b, "-------------- one step: filling in black points --------------\n");

		
		// Exchange data after black points got updated
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

		
		// Compute the residual every few steps (e.g. after every 5 or 10 iterations)
		if (k % 1 == 0){

			max_res = 0;
			cuda_compute_residual(u_p, x_p, y_p, r_p, grid);

			cuda_reduce_max(r_p, &max_res, grid.nd_l);

			cuda_synchronize();

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
		printf("run time: %.3f(s)\n", run_time);
		printf("Time per iter: %.4f(ms)\n", 1000.0 * run_time / solver_par.iter_count);
	}

	#if HAS_EXACT_SOLUTION
		real max_err = 0.0;
		real max_err_g = 0.0;
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

	#endif
	if (my_rank == 0){
		printf("\n");
	}

	cleanup:
		// Free all allocated vectors
		free_vectors(x_p, y_p, u_p, r_p, err_p);
		MPI_Finalize();
		return 0;

} /* main */


/* ===================================================
	Functions Implemetation
====================================================*/

/* ---------- Allocate memory ------- */
void allocate_vectors(
	real 	**x_p			/* out */,
	real 	**y_p			/* out */,
	real 	**u_p			/* out */,
	real 	**r_p			/* out */,
	real 	**err_p			/* out */,
	grid_t 	grid 			/* in  */){

	/* Allocate Unified Memory -- accessible from CPU or GPU */
	cuda_malloc_managed(x_p, grid.nd1);
	cuda_malloc_managed(y_p, grid.nd2_l);
	cuda_malloc_managed(y_p, grid.nd_l);
	cuda_malloc_managed(u_p, grid.nd_l);
	cuda_malloc_managed(r_p, grid.nd_l);
	cuda_malloc_managed(err_p, grid.nd_l);

	cuda_memset_zero(u_p, grid.nd_l);		// set all entires of u_p to be zeros
	cuda_memset_zero(r_p, grid.nd_l);		// set all entires of r_p to be zeros
	cuda_memset_zero(err_p, grid.nd_l);		// set all entries of err_p to be zeros
}	/* allocate_vectors */

/* ---------- Free memory ------- */
void free_vectors(
	real *x_p			/* in */,
	real *y_p			/* in */,
	real *u_p			/* in */,
	real *r_p			/* in */,
	real *err_p			/* in */){

	cuda_free(x_p);
	cuda_free(y_p);
	cuda_free(u_p);
	cuda_free(r_p);
	cuda_free(err_p);

}	/* free_vectors */

/* ------- Print grid function ------- */
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

/* ------- Get inputs arguments ------- */
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

/* ------- Get local index bounds ------- */
void get_local_index_bounds(
	const int 	np			/* in */,
	const int 	my_rank		/* in */,
	const int 	ni			/* in */,
	int 		*ni_l		/* out */,
	int 		*nia_l		/* out */,
	int 		*nib_l		/* out */){

	*ni_l = (ni + 1) / np - 1;
	int r = (ni + 1) % np;

	// (ni+1) may not be a multiple of np, add extra points to the processors on the right side
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

