function [U] = expand_matrix(V,par)
global geometry
if (strcmp(geometry,'1D'))
    U = zeros(par.Nx+1,1);
    U(par.I) = V;
else
    U = zeros(par.NIt,par.NJt);
    U(par.I,par.J) = V;
    if (strcmp(geometry,'Annulus'))
        U(:,par.iaGhost) = U(:,par.i2th);
        U(:,par.ibth   ) = U(:,par.iath);
    end
end
end