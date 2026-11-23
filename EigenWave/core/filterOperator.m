function y = filterOperator(x, par, fun, opts, counter)
% FILTEROPERATOR: apply the EigenWave filter operator S_h=beta(L_h) to x
%   y = filterOperator(x, par, fun, opts)
%   y = filterOperator(x, par, fun, opts, counter)  % also counts wave solves
%   x, y : vectors of active unknowns (length par.Ngi), as eigs function expects
%   counter : optional SolveCounter, incremented once per wave solve
%   Replaces myFunction.m, expand_matrix.m, shrink_matrix.m

V0 = par.fromActive(x);     % active vector -> full grid function
V1 = zeros(size(V0));       % for EigenWave: zero initial velocity
U  = waveSolve(V0, V1, par, fun, opts);     % wave solve + time filter
y  = par.toActive(U);       % full grid function -> active vector

if (nargin > 4)
    counter.increment();    % one more wave solve
end

end