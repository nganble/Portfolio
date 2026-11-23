function par = setup1D(dm, m, opts)
% SETUP1D: grid, operators, and boundary handling for the 1D wave equation
%   par = setup1D(dm, m, opts)
%   dm   : domain [dm.LB, dm.RB]
%   m    : refinement level, Nx = 20*2^m
%   opts : c, cfl, T, scheme ('Explicit'|'Implicit'), mode ('EigenWave'|'testWaveSolver')

c = opts.c;

% ---- Mesh (same names as setupMesh.m) ----
par.Nx  = 20*2^m;                       % number of cells
par.Ngx = par.Nx+1;                     % total number of grid pts
par.Ngi = par.Nx-1;                     % number of active (interior) pts
par.dx  = (dm.RB-dm.LB)/par.Nx;         % grid spacing
par.x   = linspace(dm.LB,dm.RB,par.Nx+1)';
par.iax = 1;                            % index of boundary pt at x=a
par.ibx = par.iax+par.Nx;               % index of boundary pt at x=b
par.i1x = par.iax+1;                    % index of first interior pt
par.i2x = par.ibx-1;                    % index of last interior pt
par.I   = par.i1x:par.i2x;              % indices of interior pts
par.Ib  = par.iax:par.ibx;              % indices of interior + boundary pts

par.interior = false(par.Ngx,1);        % mask used by waveSolve
par.interior(par.I) = true;

% Local copies so the handles below don't carry the whole par struct
I = par.I;  x = par.x;  dx = par.dx;  Ngx = par.Ngx;
iax = par.iax;  ibx = par.ibx;

% ---- Discrete operators ----
DpxDmx     = @(u,I) (u(I+1)-2*u(I)+u(I-1))/(dx^2);  % interior values only
par.DpxDmx = DpxDmx;
par.Lh     = @(u) fillInterior(Ngx, I, DpxDmx(u,I)); % full grid fn, 0 on boundary

% ---- Forcing (full grid fn, 0 on boundary) ----
par.forcing = @(fun,t) fillInterior(Ngx, I, fun.f(x(I),t));

% ---- Dirichlet BCs: W_0 = ga(t), W_Nx = gb(t) ----
par.applyBC         = @(u,fun,t) setDirichlet(u, fun, t, iax, ibx);
par.applyImplicitBC = par.applyBC;      % boundary rows of A are identity rows

% ---- Time-stepping (same as end of setupMesh.m) ----
par.dt = opts.cfl*dx/c;
if strcmp(opts.scheme,'Implicit') && strcmp(opts.mode,'EigenWave')
    par.Nt = 10;
    par.dt = opts.T/par.Nt;
else
    par.Nt = ceil(opts.T/par.dt);
    par.dt = opts.T/par.Nt;             % adjust dt to reach T exactly
end

% ---- Implicit matrix, eq. (11) in waveSchemes.pdf ----
if strcmp(opts.scheme,'Implicit')
    r = (1/2)*(c*par.dt/dx)^2;
    e = ones(Ngx,1);
    A = spdiags([-r*e, (1+2*r)*e, -r*e], -1:1, Ngx, Ngx);
    A(iax,:) = 0;  A(iax,iax) = 1;      % boundary row at x=a
    A(ibx,:) = 0;  A(ibx,ibx) = 1;      % boundary row at x=b
    par.dA = decomposition(A,'banded');
end

% ---- Trapezoidal weights for the discrete inner product ----
par.weights = dx*ones(Ngx,1);
par.weights([iax ibx]) = dx/2;

% ---- Full grid fn <-> active unknowns (for eigs) ----
par.toActive   = @(u) u(I);
par.fromActive = @(v) fillInterior(Ngx, I, v);

% ---- Exact discrete eigenvalues, for checking ----
n = (1:par.Ngi)';
par.lambdaExact = (2/dx)*sin(n*pi*dx/(2*(dm.RB-dm.LB)));
end

% ===== local functions =====
function u = fillInterior(N, I, vals)
u = zeros(N,1);
u(I) = vals;
end

function u = setDirichlet(u, fun, t, iax, ibx)
u(iax) = fun.ga(t);
u(ibx) = fun.gb(t);
end
