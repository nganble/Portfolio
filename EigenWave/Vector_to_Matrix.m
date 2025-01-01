function [U] = Vector_to_Matrix(V,par)
global geometry
n = size(V,1);
if (strcmp(geometry,'1D'))
    U = V;
elseif (strcmp(geometry,'2D'))
    if (mod(n,par.Ngx)==0)
        eqn = @(i,j) par.eqn(i,j);
        U = zeros(par.Ngx,par.Ngy);
        for iy = par.iay:par.iby
            for ix = par.iax:par.ibx
                ie = eqn(ix,iy);
                U(ix,iy) = V(ie);
            end
        end
    else
        eqn = @(i,j) par.eqn1(i,j);
        U = zeros(par.Nx-1,par.Ny-1);
        for iy = par.iay:par.Ny-1
            for ix = par.iax:par.Nx-1
                ie = eqn(ix,iy);
                U(ix,iy) = V(ie);
            end
        end
        
    end
    
elseif (strcmp(geometry,'Annulus'))
    if (mod(n,par.Ngr)==0)
        eqn = @(i,j) par.eqn(i,j);
        U = zeros(par.Ngr,par.Ngth);
        for ir = par.iar:par.ibr
            for ith = par.iaGhost:par.ibth
                ie = eqn(ir,ith);
                U(ir,ith) = V(ie);
            end
        end
    else
        eqn = @(i,j) par.eqn1(i,j);
        U = zeros(par.Nr-1,par.Nth);
        for ith=1:par.Nth
            for ir=1:par.Nr-1
                ie = eqn(ir,ith);
                U(ir,ith) = V(ie);
            end
        end
    end
end

end