#ifndef STENCIL_H
#define STENCIL_H

// =====================================================================
// LH(u, i1, i2): L_h u at the grid point (i1,i2), the second-order approximation of
//   div(m grad u) = (m11 u_x)_x + (m12 u_y)_x + (m21 u_x)_y + (m22 u_y)_y
// (eq. 2DAnisotropicLh in docs/anisotropicWaveScheme.tex):
//   D+x[ M11_{i-1/2,j} D-x u ] + D+y[ M22_{i,j-1/2} D-y u ] + D0x[ M12 D0y u ] + D0y[ M21 D0x u ]
// with M11, M22 averaged to the half points and M12, M21 at the grid points.
// u is an array macro (up, uc); uses the macros m11..m22 (common/grid.h) and dx[].
// =====================================================================
#define LH(u, i1, i2) ( \
	(0.5 * (m11(i1 + 1, i2) + m11(i1, i2)) * u(i1 + 1, i2) - \
	0.5 * (m11(i1 + 1, i2) + 2.0 * m11(i1, i2) + m11(i1 - 1, i2)) * u(i1, i2) + \
	0.5 * (m11(i1 - 1, i2) + m11(i1, i2)) * u(i1 - 1, i2) ) / (dx[0] * dx[0]) \
	+    (m12(i1 + 1, i2) * (u(i1 + 1, i2 + 1) - u(i1 + 1, i2 - 1)) - \
	m12(i1 - 1, i2) * (u(i1 - 1, i2 + 1) - u(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1]) \
	+    (m21(i1, i2 + 1) * (u(i1 + 1, i2 + 1) - u(i1 - 1, i2 + 1)) - \
	m21(i1, i2 - 1) * (u(i1 + 1, i2 - 1) - u(i1 - 1, i2 - 1))) / (4 * dx[0] * dx[1]) \
	+    (0.5 * (m22(i1, i2 + 1) + m22(i1, i2)) * u(i1, i2 + 1) - \
	0.5 * (m22(i1, i2 + 1) + 2.0 * m22(i1, i2) + m22(i1, i2 - 1)) * u(i1, i2) + \
	0.5 * (m22(i1, i2 - 1) + m22(i1, i2)) * u(i1, i2 - 1) ) / (dx[1] * dx[1]) )

#endif
