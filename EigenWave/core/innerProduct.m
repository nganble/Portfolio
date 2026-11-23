function p = innerProduct(U, V, par)
% INNERPRODUCT: discrete L2 inner product (U,V)_h = sum_j w_j*U_j*V_j
%   p = innerProduct(U, V, par)
%   U, V : grid functions (vectors in 1D, arrays in 2D/Annulus), same size as par.weights
%   par  : struct from setup1D / setupRectangle / setupAnnulus
%   p    : scalar
%
%   Trapezoidal rule. All geometry dependence (half weights on the boundary,
%   the Jacobian r on the annulus, zero weight on ghost points) is in par.weights.


p = sum(par.weights.*U.*V, 'all');

end
