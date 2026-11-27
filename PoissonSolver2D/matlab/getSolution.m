% Solving eliptic equation in 2D using iterative method (SOR)
% \partial_x^2 u + \partial_y^2 u = f
% The domain is a square
% Dirichlet boundary conditions


fontSize=16; lineWidth=2; markerSize=10;
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);
% clc;
clear; close all; clc;

%% Variables specification
par.ax = 0;
par.bx = 1;
par.ay = 0;
par.by = 1;
par.mx = 5;
% par.my = 1;
par.my = par.mx;
par.Nx = 10 * 2 ^ (par.mx - 1);
par.Ny = 10 * 2 ^ (par.my - 1);
par.Ng = (par.Nx + 1) * (par.Ny + 1);
par.dx = (par.bx - par.ax) / par.Nx;
par.dy = (par.by - par.ay) / par.Ny;
par.h = par.dx;
par.I = (2 : par.Nx)';
par.J = (2 : par.Ny)';
% par.I_red = par.I(mod(par.I, 2) == 0);
% par.I_black = par.I(mod(par.I, 2) == 1);
par.I_even = par.I(1 : 2 : end);
par.I_odd = par.I(2 : 2 : end);
par.J_even = par.J(1 : 2 : end);
par.J_odd = par.J(2 : 2 : end);
par.x = linspace(par.ax, par.bx, par.Nx + 1)';
par.y = linspace(par.ay, par.by, par.Ny + 1)';
[par.X, par.Y] = meshgrid(par.x, par.y);
% meshgrid returns length(y) * length(x) dimesions
par.X = par.X';
par.Y = par.Y';
par.lambda = 100; % parameter for the forcing
par.kx = 2;
par.ky = 3.5;
par.x0 = 0.5;
par.y0 = 0.5;

% Difference operators
par.DpxDmx = @(U, I, J) U(I + 1, J) - 2 * U(I, J) + U(I - 1, J);
par.DpyDmy = @(U, I, J) U(I, J + 1) - 2 * U(I, J) + U(I, J - 1);

%% Specifying the forcing and exact solution for testing
% fun is a struct of function handles
fun.f = @(x,y,lambda) (10*par.lambda/pi) * exp(-par.lambda^2*((x-0.5).^2 + (y-0.5).^2)) + ...
                  (10*par.lambda/pi) * exp(-par.lambda^2*((x-0.25).^2 + (y-0.25).^2));

% Plot the forcing
Z = fun.f(par.X,par.Y);
figure;
surf(par.X,par.Y, Z);

xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Forcing term f(x,y)');
colormap(getRainbow);
shading interp;   % smooth surface
colorbar;         % show color scale
savePlotFile("../docs/fig/forcing");    

%% Struct to store the iterative solver
iter.max = 30000;
iter.tol = 1e-10;
iter.omega = 1.9;         % omega = 1 is the Gauss-Seidel method
iter.iter_counter = 0;

ukp1 = zeros(par.Nx + 1, par.Ny + 1);
uk = zeros(par.Nx + 1, par.Ny + 1);
R = zeros(par.Nx + 1, par.Ny + 1);

%% Gauss-Seidel method
% Reset to initial guess to zeros
for k = 1 : iter.max
    % Update red points first
    uk(par.I_even, par.J_even) = (1 - iter.omega) * uk(par.I_even, par.J_even) + ...
        (iter.omega / 4) * ( uk(par.I_even + 1, par.J_even) + ...
        uk(par.I_even - 1, par.J_even) + ...
        uk(par.I_even, par.J_even + 1) + ...
        uk(par.I_even, par.J_even - 1)   ...
        - par.h ^ 2 * fun.f(par.X(par.I_even, par.J_even), par.Y(par.I_even, par.J_even)) );
    
    uk(par.I_odd, par.J_odd) =  (1 - iter.omega) * uk(par.I_odd, par.J_odd) + ...
        (iter.omega / 4) * ( uk(par.I_odd + 1, par.J_odd) + ...
        uk(par.I_odd - 1, par.J_odd) + ...
        uk(par.I_odd, par.J_odd + 1) + ...
        uk(par.I_odd, par.J_odd - 1)   ...
        - par.h ^ 2 * fun.f(par.X(par.I_odd, par.J_odd), par.Y(par.I_odd, par.J_odd)) );

    % Update black points
    uk(par.I_even, par.J_odd) =  (1 - iter.omega) * uk(par.I_even, par.J_odd) + ...
        (iter.omega / 4) * ( uk(par.I_even + 1, par.J_odd) + ...
        uk(par.I_even - 1, par.J_odd) + ...
        uk(par.I_even, par.J_odd + 1) + ...
        uk(par.I_even, par.J_odd - 1)   ...
        - par.h ^ 2 * fun.f(par.X(par.I_even, par.J_odd), par.Y(par.I_even, par.J_odd)) );
    
    uk(par.I_odd, par.J_even) =  (1 - iter.omega) * uk(par.I_odd, par.J_even) + ...
        (iter.omega / 4) * ( uk(par.I_odd + 1, par.J_even) + ...
        uk(par.I_odd - 1, par.J_even) + ...
        uk(par.I_odd, par.J_even + 1) + ...
        uk(par.I_odd, par.J_even - 1)   ...
        - par.h ^ 2 * fun.f(par.X(par.I_odd, par.J_even), par.Y(par.I_odd, par.J_even)) );
    
    % Need to add the way to compute the residual
    R(par.I, par.J) = (1 / 4) * ( par.DpxDmx(uk, par.I, par.J) + par.DpyDmy(uk, par.I, par.J) ) - ...
        (par.h ^ 2 / 4) * fun.f(par.X(par.I, par.J), par.Y(par.I, par.J));
    max_res = max(max(abs(R)));
    if max_res < iter.tol
        fprintf("Gauss-Seidel solver converges after %d iterations!\n", k);
        break;
    end
end
if k == iter.max
    fprintf('Gauss-Seidel Iterative solver did not converged! max_res=%4.2e\n', max_res);
end

% Plot the solution
figure;
surf(par.X,par.Y, uk);

xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Computed solution u(x,y)');
colormap(getRainbow);
shading interp;   % smooth surface
colorbar;         % show color scale
savePlotFile("../docs/fig/computed_solution");