function [U, Ufinal] = waveSolve(V0, V1, par, fun, opts)
% WAVESOLVE: solve w_tt = c^2*Lh(w) + f on [0, T] and apply the time filter.
%   Schemes follow waveSchemes.pdf: explicit eq. (2), implicit eq. (8)-(18).
%
%   [U, Ufinal] = waveSolve(V0, V1, par, fun, opts)
%   V0, V1 : initial conditions (grid functions, same size as par.x)
%   par    : struct from setup1D / setupRectangle / setupAnnulus
%   fun    : boundary data and forcing (fun.f, fun.ga, fun.gb from getMS)
%   opts   : omega, T, c, scheme ('Explicit'|'Implicit')
%   U      : filtered solution (2/T)*dt*sum_n' (cos(omega*t^n) - 1/4)*W^n
%            (trapezoidal sum: half weight at n=0 and n=Nt)
%   Ufinal : solution W^Nt at t = T (used to test the wave solver)

% Simplifying notation
omega  = opts.omega;
T      = opts.T;
c      = opts.c;
scheme = opts.scheme;
dt     = par.dt;
Nt     = par.Nt;
Lh     = par.Lh;                          % Lh(u)    = D+D- u (0 on boundary)
F      = @(t) par.forcing(fun,t);         % F(t)     = f(x_j,t) (0 on boundary)
setBC  = @(u,t) par.applyBC(u,fun,t);     % set boundary values at time t

% Allocate space for solutions at three levels
unm1 = zeros(size(V0));
un   = zeros(size(V0));
unp1 = zeros(size(V0));

% Sol at t = 0, eq. (2b)/(8b)
unm1(par.interior) = V0(par.interior);
unm1 = setBC(unm1,0);

% Apply time filter at t = 0
U = (1/2)*(cos(omega*0*dt)-1/4)*unm1;

% Sol at t = dt
if(strcmp(scheme,'Explicit'))
    % Taylor expansion, eq. (2c)
    un = V0+dt*V1+(dt^2/2)*(c^2*Lh(V0)+F(0));
    un = setBC(un,dt);
elseif(strcmp(scheme,'Implicit'))
    % Eliminate W^{-1} using D0t W^0 = V1, eq. (17)-(18)
    rhs = V0+dt*V1-(c^2*dt^3/2)*Lh(V1)+(dt^2/4)*(F(dt)+F(-dt));
    rhs = par.applyImplicitBC(rhs,fun,dt);
    un  = implicitSolve(par,rhs);
end

%------------------------Start time-stepping loop---------------------%
for n = 1:Nt-1
    t = n*dt; % current time
    if(strcmp(scheme,'Explicit'))
        % Eq. (2a)
        unp1 = 2*un-unm1+(dt^2)*(c^2*Lh(un)+F(t));
        unp1 = setBC(unp1,t+dt);
    elseif(strcmp(scheme,'Implicit'))
        % Eq. (10)-(12): A*W^{n+1} = G^n
        rhs = 2*un-unm1+((c*dt)^2/2)*Lh(unm1)+(dt^2/2)*(F(t+dt)+F(t-dt));
        rhs = par.applyImplicitBC(rhs,fun,t+dt);
        unp1 = implicitSolve(par,rhs);
    end
    % Apply the time filter at t = n*dt
    U = U+(cos(omega*n*dt)-1/4)*un;
    % Update unm1 and un for next step
    unm1 = un;
    un   = unp1;
end
%------------------------End time-stepping loop-----------------------%

% Apply time filter at final time t = Nt*dt (un holds W^Nt here)
U = U+(1/2)*(cos(omega*Nt*dt)-1/4)*un;
U = (2/T)*dt*U;

Ufinal = un;

end

% ===== local functions =====
function u = implicitSolve(par,rhs)
% Solve A*u = rhs; works for vectors (1D) and grid arrays (2D)
u = reshape(par.dA\rhs(:),size(rhs));
end
