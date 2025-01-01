function par = setupMesh(dm,m)
global c scheme option geometry cfl T

if (strcmp(geometry,'1D'))
    % Mesh in 1D
    par.Nx  = 20*2^m;
    par.NIt = par.Nx+1; % Total number of grid pts
    par.NIa = par.Nx-1; % Number of active pts
    par.NJt = 1;
    par.NJa = 1;
    par.Ng  = par.NIt*par.NJt; 
    par.Ngx = par.Nx+1;
    par.Ngi = par.Nx-1; % only active pts
    par.dx  = (dm.RB-dm.LB)/par.Nx;
    par.x   = linspace(dm.LB,dm.RB,par.Nx+1)';
    par.iax = 1; % index of boundary point at x=a
    par.ibx = par.iax+par.Nx; % index of boundary point at x=b
    par.i1x = par.iax+1; % first interior point
    par.i2x = par.ibx-1; % last interior point
    par.Ib  = par.iax:par.ibx; % interior points & boundary points
    par.Jb  = 1;
    par.I   = par.i1x:par.i2x; % interior points
    par.J   = 1;
    % Define D+D- operator
    par.DpxDmx = @(u,I) (u(I+1)-2*u(I)+u(I-1))/(par.dx^2);
    par.Lh     = @(u,I) par.DpxDmx(u,I);
    % Time-stepping
    par.dt = cfl*par.dx/c;
    
elseif (strcmp(geometry,'2D'))
    % Mesh in 2D
    par.Nx = 20*2^(m-1);
    par.Ny = 20*2^(m-1);
    par.dx = (dm.RB-dm.LB)/par.Nx;
    par.dy = (dm.TB-dm.BB)/par.Ny;
    par.iax = 1;
    par.ibx = par.iax+par.Nx;
    par.iay = 1;
    par.iby = par.iay+par.Ny;
    par.i1x = par.iax+1;
    par.i2x = par.ibx-1;
    par.i1y = par.iay+1;
    par.i2y = par.iby-1;
    par.I  = par.i1x:par.i2x;
    par.Ib = par.iax:par.ibx;
    par.J  = par.i1y:par.i2y;
    par.Jb = par.iay:par.iby;
    par.Ngx = par.Nx+1; % total number or grid points in x
    par.Ngy = par.Ny+1; % total number of grid points in y
    par.NIt = par.Nx+1;
    par.NJt = par.Ny+1;
    par.NIa = par.Nx-1; % number of active pts in x
    par.NJa = par.Ny-1; % number of active pts in y
    par.Ng  =  par.Ngx*par.Ngy;
    par.Ngi = (par.Nx-1)*(par.Ny-1);
    par.x = zeros(par.Ngx,par.Ngy);
    par.y = zeros(par.Ngx,par.Ngy);
    % Fill in x and y
    for iy = 1:par.Ngy
        for ix = 1:par.Ngx
            par.x(ix,iy) = dm.LB+(ix-par.iax)*par.dx;
            par.y(ix,iy) = dm.BB+(iy-par.iay)*par.dy;
        end
    end
    % Grid function to vector index
    par.eqn  = @(ix,iy) 1+ix-par.iax+par.Ngx*(iy-par.iay); % include boundary
    par.eqn1 = @(ix,iy) 1+ix-par.iax+(par.Nx-1)*(iy-par.iay);
    % Discrete operators
    par.DpxDmx = @(w,I,J) (w(I+1,J)-2*w(I,J)+w(I-1,J))/(par.dx^2);
    par.DpyDmy = @(w,I,J) (w(I,J+1)-2*w(I,J)+w(I,J-1))/(par.dy^2);
    par.Lh     = @(w,I,J) par.DpxDmx(w,I,J)+par.DpyDmy(w,I,J);
    % Time-stepping
    par.dt = (cfl/c)/sqrt(1/(par.dx^2)+1/(par.dy^2));

elseif (strcmp(geometry,'Annulus'))
    % Mesh on an annulus
    par.Nr  = 5*2^m;
    par.Nth = 10*2^m;
    par.dr = (dm.RB-dm.LB)/par.Nr; % grid spacing in r
    par.dth = (dm.TB-dm.BB)/par.Nth; % grid spacing in theta
    par.rNumGhost  = 0;
    par.thNumGhost = 1;
    par.iar  = 1+par.rNumGhost; % index of boundary point at r=a
    par.ibr  = par.iar+par.Nr; % index of boundary point at r=b
    
    par.iath = 1+par.thNumGhost; % index of boundary pt at theta=0
    par.ibth = par.iath+par.Nth; % index of boundary pt at theta=2*pi
    par.iaGhost = par.iath-1;
    par.i1r = par.iar+1; % index of first interior pt in r
    par.i2r = par.ibr-1; % index of last interior pt in r
    par.i1th = par.iath+1; % index of first interior pt in theta
    par.i2th = par.ibth-1; % index of last interior pt in theta
    par.I   = par.i1r:par.i2r; % interior pts in r
    par.J   = par.iath:par.i2th; % interior pts in theta
    par.Ib  = par.iar:par.ibr; % interior + boundary pts in r
    par.JbL = par.iath:par.i2th; % interior + left boundary pts in theta
    par.Jb  = par.iath:par.ibth; % interior + boundary pts in theta
    par.Ngr = par.Nr+1+par.rNumGhost; % total number of grid pts in r
    par.Ngth = par.Nth+1+par.thNumGhost; % total number of grid points in theta
    par.NIt = par.Nr+1+par.rNumGhost;
    par.NJt = par.Nth+1+par.thNumGhost;
    par.NIa = par.Nr-1; % number of active pts in r
    par.NJa = par.Nth; % number of active pts in theta
    par.Ng  = par.Ngr*par.Ngth; % total number of mesh pts
    par.Ngi = (par.Nr-1)*(par.Nth); % number of active pts
    par.r  = zeros(par.Ngr,par.Ngth);
    par.th = zeros(par.Ngr,par.Ngth);
    % Fill in r and theta
    for ith = 1:par.Ngth
        for ir = 1:par.Ngr
            par.r (ir,ith) = dm.LB+(ir-1)*par.dr;
            par.th(ir,ith) = dm.BB-par.thNumGhost*par.dth+(ith-1)*par.dth;
        end
    end
    % Convert to x&y coordinates
    par.u = @(r,theta) r.*cos(theta);
    par.v = @(r,theta) r.*sin(theta);
    par.x = par.u(par.r,par.th);
    par.y = par.v(par.r,par.th);
    % Grid function to vector index
    par.eqn  = @(ir,ith) 1+ir-par.iar+par.Ngr*(ith-1); % include boundary pts
    par.eqn1 = @(ir,ith) 1+ir-par.iar+(par.Nr-1)*(ith-1); % oinly interior pts
    % Discrete operators
    par.DprDmr   = @(w,I,J) (w(I+1,J)-2*w(I,J)+w(I-1,J))/(par.dr^2);
    par.D0r      = @(w,I,J) (w(I+1,J)-w(I-1,J))/(2*par.dr);
    par.DpthDmth = @(w,I,J) (w(I,J+1)-2*w(I,J)+w(I,J-1))/(par.dth^2);
    par.Lh       = @(w,I,J) par.DprDmr(w,I,J)+(1./par.r(I,J)).*par.D0r(w,I,J)+...
                   (1./(par.r(I,J).^2)).*par.DpthDmth(w,I,J);
    % Time-stepping
    par.dt = (cfl/c)/sqrt(1/(par.dr^2)+1/(par.dth^2));
    
end

if ((strcmp(scheme,'Implicit'))&&(strcmp(option,'EigenWave')))
    par.Nt = 10;
    par.dt = T/par.Nt;
else
    par.Nt = round(T/par.dt);
    
    par.dt = T/par.Nt;  % adjust dt to reach tFinal exactly
end
end