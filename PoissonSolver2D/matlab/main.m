% Solving eliptic equation in 2D using iterative method (SOR)
% \partial_x^2 u + \partial_y^2 u = f
% The domain is a square
% Dirichlet boundary conditions

% clc;
clear; close all;

%% Variables specification
par.ax = 0;
par.bx = 1;
par.ay = 0;
par.by = 1;
par.mx = 4;
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
fun.ue = @(x, y) exp(par.kx * (x - par.x0) .^ 2) .* exp(par.ky * (y - par.y0) .^ 2);
fun.ue_xx = @(x, y) (4 * par.kx ^ 2 * (x - par.x0) .^ 2 + 2 * par.kx) .* fun.ue(x, y);
fun.ue_yy = @(x, y) (4 * par.ky ^ 2 * (y - par.y0) .^ 2 + 2 * par.ky) .* fun.ue(x, y);
fun.f = @(x, y) fun.ue_xx(x, y) + fun.ue_yy(x, y); 

% Specifying boundary conditions
fun.gax = @(y) fun.ue(par.ax, y);
fun.gbx = @(y) fun.ue(par.bx, y);
fun.gay = @(x) fun.ue(x, par.ay);
fun.gby = @(x) fun.ue(x, par.by);



%% Struct to store the iterative solver
iter.max = 30000;
iter.tol = 1e-10;
iter.omega = 1.9;         % omega = 1 is the Gauss-Seidel method
iter.iter_counter = 0;

ukp1 = zeros(par.Nx + 1, par.Ny + 1);
uk = zeros(par.Nx + 1, par.Ny + 1);

% Fill in the boundary points
% left
ukp1(1, 1 : par.Ny + 1) = fun.gax(par.y); 
uk(1, 1 : par.Ny + 1) = fun.gax(par.y);
% right
ukp1(par.Nx + 1, 1 : par.Ny + 1) = fun.gbx(par.y); 
uk(par.Nx + 1, 1 : par.Ny + 1) = fun.gbx(par.y);
% bottom
ukp1(1 : par.Nx + 1, 1) = fun.gay(par.x); 
uk(1 : par.Nx + 1, 1) = fun.gay(par.x);

% top
ukp1(1 : par.Nx + 1, par.Ny + 1) = fun.gby(par.x); 
uk(1 : par.Nx + 1, par.Ny + 1) = fun.gby(par.x);
R = zeros(par.Nx + 1, par.Ny + 1);
% Exact solution for testing
ue = fun.ue(par.X, par.Y);

%% Jacobi method
max_res = 1;
for k = 1 : iter.max
    R(par.I, par.J) = (1 / 4) * ( par.DpxDmx(uk, par.I, par.J) + par.DpyDmy(uk, par.I, par.J) ) - ...
        (par.h ^ 2 / 4) * fun.f(par.X(par.I, par.J), par.Y(par.I, par.J));
    ukp1(par.I, par.J) = uk(par.I, par.J) + R(par.I, par.J);
    uk = ukp1;
    max_res = max(max(abs(R)));
    if max_res < iter.tol
        fprintf("Jacobi solver converges after %d iterations!\n", k);
        break;
    end
end
if k == iter.max
    fprintf('Jacobi Iterative solver did not converged! max_res=%4.2e\n', max_res);
end
max_err = max(max(abs(ue - uk)));
fprintf("Max error: %4.2e\n", max_err);

%% Gauss-Seidel method
% Reset to initial guess to zeros
uk(par.I, par.J) = 0;
R = zeros(par.Nx + 1, par.Ny + 1);
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
% uk = zeros(par.Nx + 1, par.Ny + 1);
% % Red points
% uk(par.I_even, par.J_even) = fun.ue(par.X(par.I_even, par.J_even), par.Y(par.I_even, par.J_even));
% uk(par.I_odd, par.J_odd) = fun.ue(par.X(par.I_odd, par.J_odd), par.Y(par.I_odd, par.J_odd));
% % Black points
% uk(par.I_even, par.J_odd) = fun.ue(par.X(par.I_even, par.J_odd), par.Y(par.I_even, par.J_odd));
% uk(par.I_odd, par.J_even) = fun.ue(par.X(par.I_odd, par.J_even), par.Y(par.I_odd, par.J_even));

max_err = max(max(abs(ue - uk)));
fprintf("Max error: %4.2e\n", max_err);