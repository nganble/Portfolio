#ifndef POISSON2D_CUDA_H
#define POISSON2D_CUDA_H

#include "poisson2d_common.h"

#ifdef __cplusplus

extern "C" {
#endif

/* ===================================================
	Memory management
====================================================*/
void cuda_malloc_managed(real **ptr, size_t n);
void cuda_memset_zero(real **ptr, size_t n);
void cuda_free(real *ptr);

/* ===================================================
	Device control
====================================================*/
void cuda_get_error();
void cuda_mem_prefetch(real *ptr, size_t n, int device);
void cuda_synchronize();
void cuda_get_device_count(int *num_device);
void cuda_get_current_device(int *device);
void cuda_set_device(int device);

/* ===================================================
	Kernels Wrappers
====================================================*/
void cuda_red_black_step(real *u_p, real *x_p, real *y_p, 
					grid_t grid, real omega, int color);

void cuda_compute_residual(real *u_p, real *x_p, real *y_p, 
					real *r_p, grid_t grid);

void cuda_reduce_max(real *u_p, real *block_max, int nd);


#ifdef __cplusplus
}

#endif


#endif /* POISSON2D_CUDA_H */
