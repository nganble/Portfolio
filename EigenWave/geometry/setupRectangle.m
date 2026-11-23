function par = setupRectangle(dm, m, opts)
% SETUPRECTANGLE: grid, operators, and boundary handling for the 2D wave equation
%   on the rectangle [dm.LB, dm.RB] x [dm.BB, dm.TB]
%   par = setupRectangle(dm, m, opts)
%   dm   : domain bounds LB, RB (x) and BB, TB (y)
%   m    : refinement level, Nx = Ny = 20*2^(m-1)
%   opts : c, cfl, T, scheme ('Explicit'|'Implicit'), mode ('EigenWave'|'testWaveSolver')
%
%   Grid functions are Ngx-by-Ngy arrays, u(ix,iy) ~ u(x_ix, y_iy).
%   Vector ordering (for the implicit matrix and eigs) is column-major, u(:),
%   i.e. ie = ix + Ngx*(iy-1), the same as par.eqn in setupMesh.m.

c = opts.c;

% ---- Mesh (same names as setupMesh.m) ----
par.Nx  = 20*2^(m-1);                   % number of cells in x
par.Ny  = 20*2^(m-1);                   % number of cells in y
par.dx  = (dm.RB-dm.LB)/par.Nx;         % grid spacing in x
par.dy  = (dm.TB-dm.BB)/par.Ny;         % grid spacing in y
par.iax = 1;                            % index of boundary pt at x=a
par.ibx = par.iax+par.Nx;               % index of boundary pt at x=b
par.iay = 1;                            % index of boundary pt at y=c
par.iby = par.iay+par.Ny;               % index of boundary pt at y=d
par.i1x = par.iax+1;                    % index of first interior pt in x
par.i2x = par.ibx-1;                    % index of last interior pt in x
par.i1y = par.iay+1;                    % index of first interior pt in y
par.i2y = par.iby-1;                    % index of last interior pt in y
par.I   = par.i1x:par.i2x;              % indices of interior pts in x
par.J   = par.i1y:par.i2y;              % indices of interior pts in y
par.Ib  = par.iax:par.ibx;              % indices of interior + boundary pts in x
par.Jb  = par.iay:par.iby;              % indices of interior + boundary pts in y
par.Ngx = par.Nx+1;                     % total number of grid pts in x
par.Ngy = par.Ny+1;                     % total number of grid pts in y
par.Ng  = par.Ngx*par.Ngy;              % total number of grid pts
par.Ngi = (par.Nx-1)*(par.Ny-1);        % number of active (interior) pts

% Grid coordinates (same formula as setupMesh.m)
xv = dm.LB+((par.Ib-par.iax)')*par.dx;
yv = dm.BB+((par.Jb-par.iay)')*par.dy;
[par.x, par.y] = ndgrid(xv, yv);        % Ngx-by-Ngy arrays

par.interior = false(par.Ngx,par.Ngy);  % mask used by waveSolve
par.interior(par.I,par.J) = true;

% Local copies so the handles below don't carry the whole par struct
I = par.I;  J = par.J;  Ib = par.Ib;  Jb = par.Jb;
dx = par.dx;  dy = par.dy;  Ngx = par.Ngx;  Ngy = par.Ngy;
iax = par.iax;  ibx = par.ibx;  iay = par.iay;  iby = par.iby;
x = par.x;  y = par.y;
xI = x(I,J);  yI = y(I,J);              % interior coordinates (for forcing)

% ---- Discrete operators ----
DpxDmx     = @(w,I,J) (w(I+1,J)-2*w(I,J)+w(I-1,J))/(dx^2);
DpyDmy     = @(w,I,J) (w(I,J+1)-2*w(I,J)+w(I,J-1))/(dy^2);
par.DpxDmx = DpxDmx;
par.DpyDmy = DpyDmy;
par.Lh     = @(w) fillInterior(Ngx, Ngy, I, J, DpxDmx(w,I,J)+DpyDmy(w,I,J));

% ---- Forcing (full grid fn, 0 on boundary) ----
par.forcing = @(fun,t) fillInterior(Ngx, Ngy, I, J, fun.f(xI,yI,t));

% ---- Dirichlet BCs on the four sides ----
bnd.iax = iax;  bnd.ibx = ibx;  bnd.iay = iay;  bnd.iby = iby;
bnd.Ib  = Ib;   bnd.Jb  = Jb;
bnd.yL  = y(iax,Jb);  bnd.yR = y(ibx,Jb);    % y along x=a and x=b
bnd.xB  = x(Ib,iay);  bnd.xT = x(Ib,iby);    % x along y=c and y=d
par.applyBC         = @(u,fun,t) setDirichlet(u, fun, t, bnd);
par.applyImplicitBC = par.applyBC;      % boundary rows of A are identity rows

% ---- Time-stepping (same as end of setupMesh.m, with ceil) ----
par.dt = (opts.cfl/c)/sqrt(1/(dx^2)+1/(dy^2));
if strcmp(opts.scheme,'Implicit') && strcmp(opts.mode,'EigenWave')
    par.Nt = 10;
    par.dt = opts.T/par.Nt;
else
    par.Nt = ceil(opts.T/par.dt);
    par.dt = opts.T/par.Nt;             % adjust dt to reach T exactly
end

% ---- Implicit matrix (same as formImplicitMatrix.m) ----
%   interior rows: (1+2r+2s) W_ij - r(W_i-1,j + W_i+1,j) - s(W_i,j-1 + W_i,j+1)
%   with r = (c dt/dx)^2/2, s = (c dt/dy)^2/2, i.e. A = I - (c dt)^2/2 * L_h;
%   boundary rows: identity
par.timeFactor = 0;                    % no matrix to factorize (explicit)
if strcmp(opts.scheme,'Implicit')
    Ng  = par.Ng;
    ex  = ones(Ngx,1);  ey = ones(Ngy,1);
    Dxx = spdiags([ex -2*ex ex], -1:1, Ngx, Ngx)/dx^2;
    Dyy = spdiags([ey -2*ey ey], -1:1, Ngy, Ngy)/dy^2;
    L   = kron(speye(Ngy),Dxx) + kron(Dyy,speye(Ngx));   % 5-point Laplacian, u(:) ordering
    P   = spdiags(double(par.interior(:)), 0, Ng, Ng);   % 1 on interior rows, 0 on boundary rows
    A   = P*(speye(Ng) - ((c*par.dt)^2/2)*L) + (speye(Ng) - P);
    tFactor = tic;
    par.dA = decomposition(A);
    par.timeFactor = toc(tFactor);     % time to factorize A
    par.sizeA = size(A,1);              % number of unknowns in A
end

% ---- Trapezoidal weights for the discrete inner product ----
wx = dx*ones(Ngx,1);  wx([iax ibx]) = dx/2;
wy = dy*ones(Ngy,1);  wy([iay iby]) = dy/2;
par.weights = wx*wy.';                  % Ngx-by-Ngy

% ---- Full grid fn <-> active unknowns (for eigs) ----
%   active ordering: column-major over (I,J), same as par.eqn1 in setupMesh.m
nI = numel(I);  nJ = numel(J);
par.toActive   = @(u) reshape(u(I,J), [], 1);
par.fromActive = @(v) fillInterior(Ngx, Ngy, I, J, reshape(v, nI, nJ));

% ---- Exact discrete eigenvalues, for checking ----
Lx = dm.RB-dm.LB;  Ly = dm.TB-dm.BB;
[kx, ky] = ndgrid(1:par.Nx-1, 1:par.Ny-1);
lam = sqrt((4/dx^2)*sin(kx*pi*dx/(2*Lx)).^2 + (4/dy^2)*sin(ky*pi*dy/(2*Ly)).^2);
par.lambdaExact = lam(:);
end

% ===== local functions =====
function u = fillInterior(Ngx, Ngy, I, J, vals)
u = zeros(Ngx,Ngy);
u(I,J) = vals;
end

function u = setDirichlet(u, fun, t, bnd)
% Same order as WaveHoltz.m: x-sides first, then y-sides (y-sides set the corners)
u(bnd.iax,bnd.Jb) = fun.gax(bnd.yL,t);   % x = a
u(bnd.ibx,bnd.Jb) = fun.gbx(bnd.yR,t);   % x = b
u(bnd.Ib,bnd.iay) = fun.gay(bnd.xB,t);   % y = c
u(bnd.Ib,bnd.iby) = fun.gby(bnd.xT,t);   % y = d
end
