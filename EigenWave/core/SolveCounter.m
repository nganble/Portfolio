classdef SolveCounter < handle
% SOLVECOUNTER: counts the applications of the filter operator (wave solves)
%   counter = SolveCounter();
%   S = @(x) filterOperator(x, par, fun, opts, counter);
%   ... eigs(S,...) / powerMethod(S,...) / arnoldi(S,...)
%   counter.n   -> number of wave solves
%
%   A handle object is used because a function handle cannot change a
%   variable of the calling script (and no globals are used).

    properties
        n = 0;      % number of wave solves
    end

    methods
        function increment(obj)
            obj.n = obj.n + 1;
        end
        function reset(obj)
            obj.n = 0;
        end
    end
end
