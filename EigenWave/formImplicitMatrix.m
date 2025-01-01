function [dA,time] = formImplicitMatrix(par)
global c geometry
cpu0 = cputime;
% Simplify notations
dt = par.dt;

if (strcmp(geometry,'2D')||strcmp(geometry,'Annulus'))
    eqn = @(ix,iy) par.eqn(ix,iy);
end

if (strcmp(geometry,'1D'))
    %% Implicit matrix for 1D
    dx = par.dx;
    A = sparse(par.Ngx);
    r = (1/2)*(c*dt/dx)^2;
    for i = 1:par.Ngx
        if((i==1)||(i==par.Ngx))
            A(i,i)=1;
        else
            A(i,i-1) = -r;
            A(i,i  ) = 1+2*r;
            A(i,i+1) = -r;
        end
    end
    
elseif (strcmp(geometry,'2D'))
    %% Implicit matrix for 2D
    % Simplify notations
    dx = par.dx; dy = par.dy;
    i1x = par.i1x; i2x = par.i2x;
    i1y = par.i1y; i2y = par.i2y;
    iax = par.iax; ibx = par.ibx;
    iay = par.iay; iby = par.iby;
    r = (1/2)*(c*dt/dx)^2;
    s = (1/2)*(c*dt/dy)^2;
    A = sparse(par.Ngx,par.Ngy);
    for ix = i1x:i2x
        for iy = i1y:i2y
            ie = eqn(ix,iy);
            A(ie,eqn(ix  ,iy-1)) = -s;
            A(ie,eqn(ix-1,iy  )) = -r;
            A(ie,ie            ) = 1+2*r+2*s;
            A(ie,eqn(ix+1,iy  )) = -r;
            A(ie,eqn(ix  ,iy+1)) = -s;
        end
    end
    for iy = iay:iby
        ie = eqn(iax,iy);
        A(ie,ie) = 1;
        ie = eqn(ibx,iy);
        A(ie,ie) = 1;
    end
    for ix = iax:ibx
        ie = eqn(ix,iay);
        A(ie,ie) = 1;
        ie = eqn(ix,iby);
        A(ie,ie) = 1;
    end
    
elseif (strcmp(geometry,'Annulus'))
    %% Implicit matrix for annulus
    % Simplify notations
    Ngr = par.Ngr; Ngth = par.Ngth;
    dr  = par.dr;  dth  = par.dth;
    iar = par.iar; iath = par.iath;
    ibr = par.ibr; ibth = par.ibth;
    i1r = par.i1r; %i1th = par.i1th;
    i2r = par.i2r; i2th = par.i2th;
    iaGhost = par.iaGhost;
    A = sparse(Ngr,Ngth);
    s = (1/2)*(c*dt/dr)^2;
    p = (1/2)*(c*dt/dth)^2;
    q = (1/4)*(c*dt)^2/dr;
    for ir = i1r:i2r
        for ith = iath:i2th
            rij = par.r(ir,ith);
            ie = eqn(ir,ith);
            A(ie,eqn(ir,ith-1)) = -p/(rij^2);
            A(ie,eqn(ir-1,ith)) = -s+q/rij;
            A(ie,ie           ) = 2*s+2*p/(rij^2)+1;
            A(ie,eqn(ir+1,ith)) = -s-q/rij;
            A(ie,eqn(ir,ith+1)) = -p/(rij^2);
        end
    end
    for ith = iaGhost:ibth
        ie = eqn(iar,ith);
        A(ie,ie) = 1;
        ie = eqn(ibr,ith);
        A(ie,ie) = 1;
    end
    for ir = iar:ibr
        ie = eqn(ir,iaGhost);
        A(ie,ie) = 1;
        ie = eqn(ir,ibth);
        A(ie,ie) = 1;
        ie = eqn(ir,iaGhost);
        A(ie,eqn(ir,i2th)) = -1;
        ie = eqn(ir,ibth);
        A(ie,eqn(ir,iath)) = -1;
    end
end
dA = decomposition(A,'banded');

time = cputime-cpu0;

end