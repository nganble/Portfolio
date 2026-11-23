function a = normW(U, par)
% NORMW: discrete L2 norm ||U||_h = sqrt((U,U)_h)
%   a = normW(U, par)
%   Uses the weighted inner product in innerProduct.m.
%   Replaces dis_norm.m.

a = sqrt(innerProduct(U, U, par));

end
