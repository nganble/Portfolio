// gcc -O3 poisson2d_serial.c -o poisson2d_serial
// ./poisson2d_serial 8 1.0 1e-10 20000
// ./poisson2d_serial 1024 2.0 1e-10 20000
// ----------------------------------------------------------
//
// 			POISSON EQUATION IN 2D (SERIAL VERSION)
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


/* Shared types, grid macros, test cases (TEST, SOLUTION) and timer */
#include "poisson2d_common.h"
#include "timing.h"


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

/* ------- Main function -------------- */ 
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
	printf("2D Poisson Solver (SERIAL): nx=ny=%d\nomega=%g, tol=%.2e, max_iter=%d\n",
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

	double t_start = get_time();
	/* ------------------ SOR iterative solver --------------------- */
	for (int k = 1; k <= solver_par.max_iter; k++){

		solver_par.iter_count += 1;

		// Update red points first
		for (int i2 = grid.nd2a + 1; i2 <= grid.nd2b - 1; i2++){
			for (int i1 = (i2 + 1) % 2 + 1; i1 <= grid.nd1b - 1; i1 += 2){
				u_k(i1, i2) = (1 - solver_par.omega) * u_k(i1, i2) +
				(solver_par.omega / 4) * ( u_k(i1 + 1, i2) + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1) -
				grid.h_sq * forcing(x(i1), y(i2)) );		
			}
		}

		// if (debug == 1){
		// 	printf("-------------- one step: filling in red point -------------- \n");
		// 	display_grid(u_p, grid);
		// 	break;
		// }

		// Update black points
		for (int i2 = grid.nd2a + 1; i2 <= grid.nd2b - 1; i2++){
			for (int i1 = i2 % 2 + 1; i1 <= grid.nd1b - 1; i1 += 2){				
				u_k(i1, i2) = (1 - solver_par.omega) * u_k(i1, i2) +
				(solver_par.omega / 4) * ( u_k(i1 + 1, i2) + u_k(i1 - 1, i2) + u_k(i1, i2 + 1) + u_k(i1, i2 - 1) -
				grid.h_sq * forcing(x(i1), y(i2)) );				
			}
		}

		// if (debug == 1){
		// 	printf("-------------- one step: filling in black point -------------- \n");
		// 	display_grid(u_p, grid);
		// 	break;
		// }

		// Compute the residual after few steps (e.g. 5 or 10 iterations)
		if (k % 1 == 0){
			max_res = 0;
			for (int i2 = grid.nd2a + 1; i2 <= grid.nd2b - 1; i2++){
				for (int i1 = grid.nd1a + 1; i1 <= grid.nd1b - 1; i1++){
					res(i1, i2) = ( u_k(i1 + 1, i2)  + u_k(i1 - 1, i2) - 4 * u_k(i1, i2) +
					u_k(i1, i2 + 1) + u_k(i1, i2 - 1) ) - grid.h_sq * forcing(x(i1), y(i2));
					res(i1, i2) = res(i1, i2) / grid.h;
					max_res = fmax( max_res, fabs(res(i1, i2)) );
				}
			}

			if (max_res < solver_par.tol){
				printf("SOR converged after %d iterations!\nMax residual: %.2e, ", k, max_res);
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
		printf("SOR failed to converge after %d iterations!\nMax residual: %.2e, ", solver_par.iter_count, max_res);
	}

	double run_time = t_end - t_start;
	printf("solver time: %.4f(s)\n", run_time);
	printf("Time per iter: %.4f(ms)\n", 1000.0 * run_time / solver_par.iter_count);

	/* ------ Compare with exact solution (if exists) ------- */
	#if HAS_EXACT_SOLUTION
		real max_err = 0;

		real *err_p;
		err_p = calloc(grid.nd, sizeof(*err_p));

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

	cleanup:
		// Free all allocated vectors
		free_vectors(x_p, y_p, u_p, r_p);

} /* main */

void allocate_vectors(
	real 	**x_p			/* out */,
	real 	**y_p			/* out */,
	real 	**u_p			/* out */,
	real 	**r_p			/* out */,
	grid_t 	grid 			/* in  */){

	*x_p = malloc(grid.nd1 * sizeof(real));
	*y_p = malloc(grid.nd2 * sizeof(real));
	*u_p = calloc(grid.nd, sizeof(real));	// set all entires of u_p to be zeros
	*r_p = calloc(grid.nd, sizeof(real));	// set all entires of r_p to be zeros

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

	for (int i2 = grid.nd2a; i2 <= grid.nd2b; i2 ++){
		printf("v^T[:,%d]= [", i2);
		for (int i1 = grid.nd1a; i1 <= grid.nd1b; i1 ++){
			printf("%.g ", v[i1 - grid.nd1a + grid.nd1 * (i2 - grid.nd2a)]);
		}
		printf("]\n");

	}

}	/* display_grid_function */

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

	// The serial grid is one block: the local (_l) bounds used by the macros in
	// poisson2d_common.h are the global ones
	grid.nd2a_l = grid.nd2a;
	grid.nd2b_l = grid.nd2b;

	grid.nd = grid.nd1 * grid.nd2;
	grid.dx = (dom.xb - dom.xa) / nx;
	grid.dy = (dom.yb - dom.ya) / ny;
	grid.h  = grid.dx;
	grid.h_sq  = grid.h * grid.h;

	return grid;
} /* make_grid */
