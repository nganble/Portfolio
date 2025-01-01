function [W0, W1] = getICs(par)
global geometry
if (strcmp(geometry,'1D'))
    vx = @(x) 2*x.*(x-1).^2;
    W0 = vx(par.x);
    W1 = zeros(par.Ngx,1); % W1 of size Ngr*Ngth (but all zeros)
elseif (strcmp(geometry,'2D'))
    %vx = @(x,y) 2*x.*(x-pi).^2.*y.^2.*(1-y).^3;
    vx = @(x,y) cos(2*pi*x).*cos(pi^2*y);
    W0 = vx(par.x,par.y);
    W1 = zeros(par.Ngx,par.Ngy);
elseif (strcmp(geometry,'Annulus'))
    vx = @(x,y) 2*x.*(x-pi).^2.*y.^2.*(1-y).^3;
    %vx = @(x,y) cos(2*pi*x).*cos(pi^2*y);
    W0 = vx(par.x,par.y);
    W1 = zeros(par.Ngr,par.Ngth); % W1 of size Ngr*Ngth (but all zeros)
end
end