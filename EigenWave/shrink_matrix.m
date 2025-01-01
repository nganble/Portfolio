function [V] = shrink_matrix(U,par)
global geometry
if (strcmp(geometry,'1D'))
    V = U(par.I);
else
    V = U(par.I,par.J);
end
end