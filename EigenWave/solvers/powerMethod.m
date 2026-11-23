function [beta, x, info] = powerMethod(S, x0, ip, maxit, tol, verbose)
% POWERMETHOD: eigenvalue of S with the largest magnitude, by power iteration
%   [beta, x, info] = powerMethod(S, x0, ip, maxit, tol, verbose)
%   S       : function handle, y = S(x) (one wave solve for EigenWave)
%   x0      : initial guess (vector of active unknowns)
%   ip      : function handle, weighted inner product ip(a,b)
%   maxit   : maximum number of iterations
%   tol     : stop when ||x_{k+1} - x_k|| < tol
%   verbose : print convergence ratios every 10 iterations (default true)
%   beta    : eigenvalue of S (filter value)
%   x       : eigenvector, normalized in the ip-norm
%   info    : numIter, numSolves, converged, res(k), err(k) (for convergence plots)
%
%   Replaces the power-iteration loop in main.m.

if nargin < 6, verbose = true; end
nrm = @(a) sqrt(ip(a,a));

res = zeros(maxit,1);                   % residual ||S x - beta x||
err = zeros(maxit,1);                   % change in eigenvector
converged = false;

x = x0/nrm(x0);
%-----------------------Start power-iteration loop---------------------%
for h = 1:maxit
    z    = S(x);                        % apply the filter: z = S x
    beta = ip(x,z);                     % Rayleigh quotient (x is normalized)
    res(h) = nrm(z-beta*x);             % residual
    % flip sign if x and z point in opposite directions
    if (beta < 0)
        x = -x;
    end
    z = z/nrm(z);                       % normalize z
    err(h) = nrm(z-x);                  % change in eigenvector

    % Print the convergence ratios every 10 iterations
    if (verbose && h > 1 && mod(h,10) == 0)
        fprintf('Iter %d: ratio_Res=%6.5f, ratio_Vec=%6.5f\n', ...
                h, res(h)/res(h-1), err(h)/err(h-1));
    end

    x = z;                              % update (z is the newest iterate)
    if (err(h) < tol)
        converged = true;
        break;
    end
end
%-----------------------End power-iteration loop-----------------------%

if (~converged)
    warning('powerMethod: not converged after %d iterations (err = %.2e)', maxit, err(h));
end

info.numIter   = h;
info.numSolves = h;                     % one application of S per iteration
info.converged = converged;
info.res       = res(1:h);
info.err       = err(1:h);
end
