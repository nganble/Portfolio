function beta = betaContinuous(lambda, opts)
% BETACONTINUOUS: continuous EigenWave filter function
%   beta = betaContinuous(lambda, opts)
%   beta(lambda) = sinc((omega-lambda)T) + sinc((omega+lambda)T) - (1/2) sinc(lambda T)
%   with sinc(u) = sin(u)/u and sinc(0) = 1.
%   opts : omega, T
%   Replaces con_Beta.m.

omega = opts.omega;
T     = opts.T;

beta = mysinc((omega-lambda)*T) + mysinc((omega+lambda)*T) - (1/2)*mysinc(lambda*T);

end

% ===== local functions =====
function s = mysinc(u)
s = ones(size(u));          % sinc(0) = 1
k = (u ~= 0);
s(k) = sin(u(k))./u(k);
end
