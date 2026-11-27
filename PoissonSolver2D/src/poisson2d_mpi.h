#ifndef POISSON2D_MPI_H
#define POISSON2D_MPI_H

#include "poisson2d_common.h"

/* ===================================================
	Memory Management
====================================================*/
void allocate_vectors(real **x_p, real **y_p, real **u_p,
                      real **r_p, real **err_p, grid_t g);

void free_vectors(real *x_p, real *y_p, real *u_p,
                  real *r_p, real *err_p);


/* ===================================================
	Grid Utilities
====================================================*/

void display_grid(real *v, grid_t g);

grid_t make_grid(domain_t dom, int nx, int ny);

void get_local_index_bounds(int np, int my_rank, int ni,
                           int *ni_l, int *nia_l, int *nib_l);

/* ===================================================
	Input / Debug
====================================================*/

void get_arguments(int argc, char *argv[], int my_rank,
                   int *n, real *omega, real *tol, int *k);

void debug_print_ranked(int my_rank, int np, const char *msg,
                        real *u_p, grid_t grid, int ny,
                        solver_par_t solver_par);

#endif /* POISSON2D_MPI_H */