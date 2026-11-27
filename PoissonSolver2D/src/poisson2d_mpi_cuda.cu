// =====================================================================
// POISSON EQUATION IN 2D: MPI+CUDA code: CUDA kernels and wrappers
//   Solves u_xx + u_yy = f on (0,1)^2 with Dirichlet boundary conditions,
//     5-point Laplacian, red-black SOR (omega = 2.0 -> optimal omega = 2/(1+sin(pi h))).
//   Called from poisson2d_mpi_cuda.c (declarations in poisson2d_cuda.h).
// =====================================================================

#include <stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include<cuda.h>
#include<cuda_runtime.h>


#include "poisson2d_common.h"
#include "poisson2d_cuda.h"



/* Configure the number of threads and blocks */
static inline void compute_launch_config(size_t nd, int *blocks, int *threads){
	*threads = 256;
	
	*blocks = (nd + (*threads) - 1) / (*threads);

	if (*blocks < 1) {
        *blocks = 1;
    }

    if (*blocks > BLOCK_MAX){
        *blocks = BLOCK_MAX;
    }	
}


/* ===================================================
	Kernels implementation
====================================================*/

/* ---------------------------------------------------------------------------- */
__global__
void red_black_step(real *u_p, real	*x_p, real *y_p, grid_t grid, real omega, int color){

	// Unpack necessary values
	int nd1   = grid.nd1;
    int nd 	  = grid.nd_l;  	// use local value
    int nd1a  = grid.nd1a;
    int nd1b  = grid.nd1b;
    int nd2a  = grid.nd2a_l; 	// use local value
    int nd2b  = grid.nd2b_l; 	// use local value
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

	// Unpack necessary values
	int nd1   = grid.nd1;
    int nd 	  = grid.nd_l;  	// use local value
    int nd1a  = grid.nd1a;
    int nd1b  = grid.nd1b;
    int nd2a  = grid.nd2a_l; 	// use local value
    int nd2b  = grid.nd2b_l; 	// use local value
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

			res(i1, i2) = (uval - 4 * u_k(i1, i2)) - h_sq * fval;
			res(i1, i2) = res(i1, i2) / grid.h;

		}
	}
} /* compute_residual */


/* ---------------------------------------------------------------------------- */
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
} /* reduce_max */

/* ===================================================
	Wrappers to call from CPU
====================================================*/

/*----------- Memory management ------------------- */
void cuda_malloc_managed(real **ptr, size_t n){

	cudaMallocManaged(ptr, n * sizeof(real));

} /* cuda_land_allocate */

void cuda_memset_zero(real **ptr, size_t n){

	cudaMemset(*ptr, 0.0, n * sizeof(real));

} /* cuda_land_memset */

void cuda_free(real *ptr){

	cudaFree(ptr);

} /* cuda_free */


/*----------- Device control ------------------- */
void cuda_get_error(){
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess){
		printf("Launch error: %s\n", cudaGetErrorString(err));

	} else {
		printf("Kernel launched succeed!\n");

	}
} /* cuda_get_error */

void cuda_mem_prefetch(real *ptr, size_t n, int device){

	cudaMemPrefetchAsync(ptr, n * sizeof(real), device);

} /* cuda_prefetch */

void cuda_synchronize(){

	cudaDeviceSynchronize();

} /* cuda_land_synchronize */


void cuda_get_device_count(int *num_devices){
	cudaGetDeviceCount(num_devices);
}

void cuda_get_current_device(int *device){
    cudaGetDevice(device);
}

void cuda_set_device(int device){
    cudaSetDevice(device);
}


/*----------- Kernels Wrappers ------------------- */
void cuda_red_black_step(real *u_p, real *x_p, real *y_p,
					grid_t grid, real omega, int color){

	int threads, blocks;
	compute_launch_config(grid.nd_l, &blocks, &threads);

	red_black_step<<<blocks, threads>>>(u_p, x_p, y_p, grid, omega, color);

} /* cuda_red_black_step */


void cuda_compute_residual(real *u_p, real *x_p, real *y_p,
						real *r_p, grid_t grid){

	int threads, blocks;
	compute_launch_config(grid.nd_l, &blocks, &threads);
	compute_residual<<<blocks, threads>>>(u_p,r_p, x_p, y_p, grid);

} /* cuda_land_comp_residual */


void cuda_reduce_max(real *u_p, real *max_val, int nd){

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
	
} /* cuda_reduce_max */

