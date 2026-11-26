#ifndef GRID_H
#define GRID_H

#include <cstddef>

// =====================================================================
// Array layout shared by the serial, OpenMP and MPI codes
//   - one ghost layer (numGhost = 1, set in the codes)
//   - rows of length nd1 are stored with a padded stride ld = paddedStride(nd1)
//   - all 11 arrays are stored in one block with staggered starts
// The including code defines IDX(i1,i2), the position of point (i1,i2) in an
// array (the MPI code counts the rows from its first local row).
// =====================================================================

// Row stride: a multiple of 8 values (64 bytes), but not a multiple of 512 values (4 KB),
// to avoid cache conflicts between neighbouring rows
inline int paddedStride(const int nd1)
{
	int ld = ((nd1 + 7)/8)*8;
	if (ld % 512 == 0) ld += 8;
	return ld;
}

// All 11 arrays in one block. Array k starts at k*slot, where slot is nd rounded up to a
// multiple of 512 values (4 KB) plus 520 values: the arrays then start at different
// offsets within a page, so the same point (i1,i2) of different arrays does not map to
// the same cache set (this made the code about 2 times slower).
inline void allocateVectors(
			double** x_p,   double** y_p,
			double** up_p,  double** un_p,
			double** uc_p,  double** errCPU,
			double** rho_p, double** m11_p,
			double** m12_p, double** m21_p,
			double** m22_p, int nd)
{
	const int slot = ((nd + 511)/512)*512 + 520;
	double *block = new double [11*(size_t)slot];
	double **arrays[11] = {x_p, y_p, up_p, un_p, uc_p, errCPU, rho_p, m11_p, m12_p, m21_p, m22_p};
	for (int k = 0; k < 11; k++)
		*arrays[k] = block + k*(size_t)slot;
}

// Free all arrays (x_p is the start of the block that holds all arrays)
inline void freeVectors(
			double* x_p,   double* y_p,
			double* up_p,  double* un_p,
			double* uc_p,  double* errCPU_p,
			double* rho_p, double* m11_p,
			double* m12_p, double* m21_p,
			double* m22_p)
{
	delete [] x_p;
}

// Array access at point (i1,i2) (IDX is defined by the including code)
#define x(i1, i2)   x_p[0][IDX(i1, i2)]
#define y(i1, i2)   x_p[1][IDX(i1, i2)]
#define up(i1, i2)  u_p[prev][IDX(i1, i2)]
#define uc(i1, i2)  u_p[cur ][IDX(i1, i2)]
#define un(i1, i2)  u_p[next][IDX(i1, i2)]
#define err(i1, i2) errCPU_p[IDX(i1, i2)]
#define rho(i1, i2) rho_p[IDX(i1, i2)]
#define m11(i1, i2) m11_p[IDX(i1, i2)]
#define m12(i1, i2) m12_p[IDX(i1, i2)]
#define m21(i1, i2) m21_p[IDX(i1, i2)]
#define m22(i1, i2) m22_p[IDX(i1, i2)]

#endif
