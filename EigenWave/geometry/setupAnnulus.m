function par = setupAnnulus(dm, m, opts)
% SETUPANNULUS: grid, operators, and boundary handling for the 2D wave equation
%   on the annulus dm.LB <= r <= dm.RB, 0 <= theta <= 2*pi (polar grid)
%   par = setupAnnulus(dm, m, opts)
%   dm   : LB, RB (radii) and BB, TB (theta range, 0 and 2*pi)
%   m    : refinement level, Nr = 5*2^m, Nth = 10*2^m
%   opts : c, cfl, T, scheme ('Explicit'|'Implicit'), mode ('EigenWave'|'testWaveSolver')
%
%   Grid functions are Ngr-by-Ngth arrays, u(ir,ith) ~ u(r_ir, theta_ith).
%   Columns in theta (thNumGhost = 1):
%     ith = iaGhost (=1)       : ghost point theta = -dth   (copy of ith = i2th)
%     ith = iath..i2th         : theta = 0,...,2*pi-dth     (unknowns, periodic)
%     ith = ibth               : theta = 2*pi               (copy of ith = iath)
%   Dirichlet BCs at r = LB (ir = iar) and r = RB (ir = ibr).
%   Vector ordering is column-major, u(:), i.e. ie = ir + Ngr*(ith-1) (par.eqn in setupMesh.m).

c = opts.c;

% ---- Mesh (same names as setupMesh.m) ----
par.Nr  = 5*2^m;                        % number of cells in r
par.Nth = 10*2^m;                       % number of cells in theta
par.dr  = (dm.RB-dm.LB)/par.Nr;         % grid spacing in r
par.dth = (dm.TB-dm.BB)/par.Nth;        % grid spacing in theta
par.rNumGhost  = 0;
par.thNumGhost = 1;
par.iar  = 1+par.rNumGhost;             % index of boundary pt at r=a
par.ibr  = par.iar+par.Nr;              % index of boundary pt at r=b
par.iath = 1+par.thNumGhost;            % index of pt at theta=0
par.ibth = par.iath+par.Nth;            % index of pt at theta=2*pi
par.iaGhost = par.iath-1;               % index of ghost pt at theta=-dth
par.i1r  = par.iar+1;                   % index of first interior pt in r
par.i2r  = par.ibr-1;                   % index of last interior pt in r
par.i1th = par.iath+1;                  % index of first interior pt in theta
par.i2th = par.ibth-1;                  % index of last interior pt in theta
par.I    = par.i1r:par.i2r;             % indices of interior pts in r
par.J    = par.iath:par.i2th;           % indices of unknowns in theta (periodic)
par.Ib   = par.iar:par.ibr;             % indices of interior + boundary pts in r
par.JbL  = par.iath:par.i2th;           % interior + left boundary pts in theta
par.Jb   = par.iath:par.ibth;           % interior + boundary pts in theta
par.Ngr  = par.Nr+1+par.rNumGhost;      % total number of grid pts in r
par.Ngth = par.Nth+1+par.thNumGhost;    % total number of grid pts in theta
par.Ng   = par.Ngr*par.Ngth;            % total number of grid pts
par.Ngi  = (par.Nr-1)*par.Nth;          % number of active pts

% Grid coordinates (same formula as setupMesh.m)
rv  = dm.LB+((1:par.Ngr)'-1)*par.dr;
thv = dm.BB-par.thNumGhost*par.dth+((1:par.Ngth)'-1)*par.dth;
[par.r, par.th] = ndgrid(rv, thv);      % Ngr-by-Ngth arrays
par.x = par.r.*cos(par.th);
par.y = par.r.*sin(par.th);

par.interior = false(par.Ngr,par.Ngth); % mask used by waveSolve
par.interior(par.I,par.J) = true;

% Local copies so the handles below don't carry the whole par struct
I = par.I;  J = par.J;  Ib = par.Ib;  Jb = par.Jb;
dr = par.dr;  dth = par.dth;  Ngr = par.Ngr;  Ngth = par.Ngth;
iar = par.iar;  ibr = par.ibr;  iath = par.iath;  i2th = par.i2th;
ibth = par.ibth;  iaGhost = par.iaGhost;
r = par.r;  th = par.th;
rI = r(I,J);  xI = par.x(I,J);  yI = par.y(I,J);

% ---- Discrete operators (ghost columns must hold the periodic copies) ----
DprDmr     = @(w,I,J) (w(I+1,J)-2*w(I,J)+w(I-1,J))/(dr^2);
D0r        = @(w,I,J) (w(I+1,J)-w(I-1,J))/(2*dr);
DpthDmth   = @(w,I,J) (w(I,J+1)-2*w(I,J)+w(I,J-1))/(dth^2);
par.DprDmr   = DprDmr;
par.D0r      = D0r;
par.DpthDmth = DpthDmth;
%   Lh(w) is a periodic grid function: ghost and theta=2*pi columns hold the copies
%   (as expand_matrix did in main.m), so (v, Lh v)_h is consistent with (v, v)_h.
bndP.Ib = Ib;  bndP.iaGhost = iaGhost;  bndP.i2th = i2th;  bndP.iath = iath;  bndP.ibth = ibth;
par.Lh = @(w) periodicCopy(fillInterior(Ngr, Ngth, I, J, ...
              DprDmr(w,I,J) + (1./rI).*D0r(w,I,J) + (1./rI.^2).*DpthDmth(w,I,J)), bndP);

% ---- Forcing (full grid fn, 0 on boundary and ghost columns) ----
par.forcing = @(fun,t) fillInterior(Ngr, Ngth, I, J, fun.f(xI,yI,t));

% ---- Boundary conditions ----
bnd.iar = iar;  bnd.ibr = ibr;  bnd.Ib = Ib;  bnd.Jb = Jb;
bnd.iaGhost = iaGhost;  bnd.i2th = i2th;  bnd.iath = iath;  bnd.ibth = ibth;
bnd.thL = th(iar,Jb);  bnd.thR = th(ibr,Jb);
% Explicit: Dirichlet in r, then copy the periodic columns (same order as WaveHoltz.m)
par.applyBC = @(u,fun,t) setBC(u, fun, t, bnd);
% Implicit rhs: Dirichlet in r, 0 in the periodic rows (A enforces the copies)
par.applyImplicitBC = @(u,fun,t) setImplicitBC(u, fun, t, bnd);

% ---- Time-stepping (same as end of setupMesh.m, with ceil) ----
%   The smallest spacing in theta is r_min*dth (r_min = dm.LB), not dth.
%   setupMesh.m used dth, which gives an effective CFL number > 1 for cfl = 0.9
%   (explicit scheme unstable on fine grids). Fixed here.
par.dt = (opts.cfl/c)/sqrt(1/(dr^2)+1/((dm.LB*dth)^2));
if strcmp(opts.scheme,'Implicit') && strcmp(opts.mode,'EigenWave')
    par.Nt = 10;
    par.dt = opts.T/par.Nt;
else
    par.Nt = ceil(opts.T/par.dt);
    par.dt = opts.T/par.Nt;             % adjust dt to reach T exactly
end

% ---- Implicit matrix (same as formImplicitMatrix.m) ----
%   interior rows: A = I - (c dt)^2/2 * L_h  (5-point polar stencil)
%   r-boundary rows: identity;  ghost/2*pi columns: u(ir,iaGhost)-u(ir,i2th) = 0,
%   u(ir,ibth)-u(ir,iath) = 0
if strcmp(opts.scheme,'Implicit')
    Ng  = par.Ng;
    eqn = @(ir,ith) ir + Ngr*(ith-1);
    s = (1/2)*(c*par.dt/dr)^2;
    p = (1/2)*(c*par.dt/dth)^2;
    q = (1/4)*(c*par.dt)^2/dr;

    % Interior rows
    [IR, ITH] = ndgrid(I, J);  IR = IR(:);  ITH = ITH(:);
    rr = rv(IR);
    k  = eqn(IR,ITH);
    rows = [k; k; k; k; k];
    cols = [k; eqn(IR-1,ITH); eqn(IR+1,ITH); eqn(IR,ITH-1); eqn(IR,ITH+1)];
    vals = [1+2*s+2*p./rr.^2; -s+q./rr; -s-q./rr; -p./rr.^2; -p./rr.^2];

    % r-boundary rows (theta = 0..2*pi-dth): identity
    kb   = [eqn(iar*ones(numel(J),1),J(:)); eqn(ibr*ones(numel(J),1),J(:))];
    rows = [rows; kb];  cols = [cols; kb];  vals = [vals; ones(numel(kb),1)];

    % Periodic rows for all ir: ghost = copy of i2th, 2*pi = copy of iath
    irv  = Ib(:);  one = ones(numel(irv),1);
    kg   = eqn(irv, iaGhost*one);   kp = eqn(irv, ibth*one);
    rows = [rows; kg; kg; kp; kp];
    cols = [cols; kg; eqn(irv,i2th*one); kp; eqn(irv,iath*one)];
    vals = [vals; one; -one; one; -one];

    A = sparse(rows, cols, vals, Ng, Ng);
    par.dA = decomposition(A);
end

% ---- Trapezoidal weights for the discrete inner product (same as dis_product.m) ----
%   (U,V)_h = sum r*U*V*w_r*w_th ; theta: half weights at 0 and 2*pi, 0 on the ghost
wr  = dr*ones(Ngr,1);   wr([iar ibr]) = dr/2;
wth = dth*ones(Ngth,1); wth([iath ibth]) = dth/2;  wth(iaGhost) = 0;
par.weights = r.*(wr*wth.');            % Ngr-by-Ngth (includes Jacobian r)

% ---- Full grid fn <-> active unknowns (for eigs) ----
%   active ordering: column-major over (I,J), same as par.eqn1 / creatMatrix.m
nI = numel(I);  nJ = numel(J);
par.toActive   = @(u) reshape(u(I,J), [], 1);
par.fromActive = @(v) periodicCopy(fillInterior(Ngr, Ngth, I, J, reshape(v, nI, nJ)), bnd);

% ---- Discrete eigenvalues, for checking (replaces creatMatrix.m + eig) ----
%   L_h on the active unknowns: Dirichlet in r, periodic in theta
ri = rv(I);  ri = ri(:);
iv = (1:nI)';
% row i: u(i-1)*(1/dr^2 - 1/(2 r_i dr)) - 2 u(i)/dr^2 + u(i+1)*(1/dr^2 + 1/(2 r_i dr))
Dr = sparse([iv; iv(2:end); iv(1:end-1)], ...
            [iv; iv(2:end)-1; iv(1:end-1)+1], ...
            [-2/dr^2*ones(nI,1); 1/dr^2-1./(2*ri(2:end)*dr); 1/dr^2+1./(2*ri(1:end-1)*dr)], ...
            nI, nI);
eth = ones(nJ,1);
Dth = spdiags([eth -2*eth eth], -1:1, nJ, nJ);
Dth(1,nJ) = 1;  Dth(nJ,1) = 1;          % periodic
Dth = Dth/dth^2;
Lact = kron(speye(nJ), Dr) + kron(Dth, spdiags(1./ri.^2, 0, nI, nI));
par.Lactive = Lact;

% Eigenvalues by separation of variables (instead of eig of the full matrix):
% the periodic theta-operator has eigenvalues mu_k = -(4/dth^2) sin^2(pi k/Nth),
% so L_h's eigenvalues are those of the small r-problems Dr + mu_k diag(1/r^2).
Nth = par.Nth;
mu  = -(4/dth^2)*sin(pi*(0:Nth-1)/Nth).^2;
Dr_full = full(Dr);  invr2 = diag(1./ri.^2);
ev = zeros(nI, Nth);
for kk = 1:Nth
    ev(:,kk) = eig(Dr_full + mu(kk)*invr2);
end
par.lambdaExact = sort(real(sqrt(-ev(:))));
end

% ===== local functions =====
function u = fillInterior(Ngr, Ngth, I, J, vals)
u = zeros(Ngr,Ngth);
u(I,J) = vals;
end

function u = periodicCopy(u, bnd)
% Ghost column = column at 2*pi-dth; column at 2*pi = column at 0
u(bnd.Ib,bnd.iaGhost) = u(bnd.Ib,bnd.i2th);
u(bnd.Ib,bnd.ibth)    = u(bnd.Ib,bnd.iath);
end

function u = setBC(u, fun, t, bnd)
u(bnd.iar,bnd.Jb) = fun.ga(bnd.thL,t);   % r = a
u(bnd.ibr,bnd.Jb) = fun.gb(bnd.thR,t);   % r = b
u = periodicCopy(u, bnd);
end

function u = setImplicitBC(u, fun, t, bnd)
u(bnd.iar,bnd.Jb) = fun.ga(bnd.thL,t);   % r = a
u(bnd.ibr,bnd.Jb) = fun.gb(bnd.thR,t);   % r = b
u(bnd.Ib,bnd.iaGhost) = 0;               % periodic rows of A have rhs 0
u(bnd.Ib,bnd.ibth)    = 0;
end
