function [prod] = dis_product(U,V,par)
global geometry
if (strcmp(geometry,'1D'))
    dx = par.dx;
    Ngx = par.Ngx;
    f = U.*V;
    prod = getIntegral(f,dx,Ngx);
elseif (strcmp(geometry,'2D'))
    f = U.*V;
    I = zeros(par.NIt,1);
    for i = 1: par.NIt
        I(i) = getIntegral(f(i,:),par.dy,par.NJt);
    end
    prod = getIntegral(I,par.dx,par.NIt);
elseif (strcmp(geometry,'Annulus'))
    f = U.*V.*par.r;
    f = f(par.Ib,par.Jb);
    I = zeros(par.NIt,1);
    for i = 1:par.NIt
        I(i) = getIntegral(f(i,:),par.dth,par.NJt-1);
    end
    prod = getIntegral(I,par.dr,par.NIt);
end
    

end