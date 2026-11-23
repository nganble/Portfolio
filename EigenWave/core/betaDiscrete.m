function beta = betaDiscrete(lambda, par, opts)
% BETADISCRETE: discrete EigenWave filter function (trapezoidal rule in time)
%   beta = betaDiscrete(lambda, par, opts)
%   beta(lambda) = (2/T) dt sum_{n=0}^{Nt}' (cos(omega t^n) - 1/4) cos(lambda t^n)
%   (half weights at n = 0 and n = Nt), the same sum used in waveSolve.m.
%   par  : dt, Nt
%   opts : omega, T
%   Replaces dis_Beta.m.

dt = par.dt;
Nt = par.Nt;
t  = (0:Nt)*dt;                         % time levels (row)

% Trapezoidal weights times filter weights
w = cos(opts.omega*t) - 1/4;
w([1 end]) = w([1 end])/2;

% Sum over time levels for each lambda
beta = (2/opts.T)*dt*(cos(lambda(:)*t)*w(:));
beta = reshape(beta, size(lambda));

end
