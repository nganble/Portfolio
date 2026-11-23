function [beta, X, info] = arnoldi(S, x0, ip, k, tol, reorth)
% ARNOLDI: k-step Arnoldi for S; returns the Ritz pairs that have converged
%   [beta, X, info] = arnoldi(S, x0, ip, k, tol, reorth)
%   S      : function handle, y = S(x) (one wave solve for EigenWave)
%   x0     : initial vector (vector of active unknowns)
%   ip     : function handle, weighted inner product ip(a,b)
%   k      : number of Arnoldi steps (= number of applications of S)
%   tol    : accept a Ritz pair if ||S x - beta x|| / |beta| < tol
%   reorth : true -> Gram-Schmidt twice (full reorthogonalization), default false
%   beta   : accepted Ritz values (eigenvalues of S)
%   X      : corresponding Ritz vectors (columns), normalized in the ip-norm
%   info   : numSteps, numSolves, H, Q, ritzValues, res (all Ritz pairs)
%
%   The residual uses the Arnoldi relation S Q_k = Q_{k+1} H_{k+1,k}:
%       ||S x_j - beta_j x_j|| = |h_{k+1,k}| |w_{k,j}|,
%   so no extra applications of S are needed (main.m used one per Ritz pair).
%
%   Replaces the Arnoldi loop in main.m.

if nargin < 6, reorth = false; end
nrm = @(a) sqrt(ip(a,a));

n = numel(x0);
Q = zeros(n, k+1);
H = zeros(k+1, k);
Q(:,1) = x0/nrm(x0);

%-------------------Start loop for Arnoldi algorithm------------------%
for j = 1:k
    v = S(Q(:,j));                      % apply the filter: v = S q_j
    % Modified Gram-Schmidt in the weighted inner product
    for i = 1:j
        H(i,j) = ip(Q(:,i), v);
        v = v - H(i,j)*Q(:,i);
    end
    % Optional second pass (keeps Q orthonormal, removes "ghost" copies)
    if (reorth)
        for i = 1:j
            hij = ip(Q(:,i), v);
            H(i,j) = H(i,j) + hij;
            v = v - hij*Q(:,i);
        end
    end
    H(j+1,j) = nrm(v);
    % Breakdown: the Krylov subspace is invariant, the Ritz pairs are exact
    if (H(j+1,j) < 1e-14*abs(H(1,1)))
        k = j;
        break;
    end
    Q(:,j+1) = v/H(j+1,j);
end
%-------------------End loop for Arnoldi algorithm--------------------%

% Ritz pairs from the k-by-k Hessenberg matrix
[W, D] = eig(H(1:k,1:k), 'vector');
if (max(abs(imag(D))) < 1e-10*max(abs(D)))
    D = real(D);  W = real(W);          % S is self-adjoint in ip: values are real
end
X = Q(:,1:k)*W;                         % Ritz vectors (ip-norm = ||W(:,j)||_2 = 1)

% Residual estimate for each Ritz pair (relative to |beta|)
res = abs(H(k+1,k)*W(k,:)).'./abs(D);

% Keep the converged Ritz pairs
keep = res < tol;
beta = D(keep);
X    = X(:,keep);
if (isempty(beta))
    warning('arnoldi: no Ritz pairs converged; increase k.');
end

info.numSteps   = k;
info.numSolves  = k;                    % one application of S per step
info.H          = H(1:k+1,1:k);
info.Q          = Q(:,1:k);
info.ritzValues = D;
info.res        = res;
end
