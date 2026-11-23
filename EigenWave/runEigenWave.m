%% runEigenWave.m - compute eigenvalues of L_h * u + lambda^2 * u = 0 with EigenWave
%  Replaces main.m (EigenWave part). Results are compared to reference/*.mat.


clc; clear; close all;

%% 1. Paths
addpath(genpath(fileparts(mfilename('fullpath'))));

%% 2. Parameters  (same values as the reference runs)
geometry   = '1D';                   % '1D' | '2D' | 'Annulus'
methodName = 'Arnoldi';                 % 'eigs' | 'Power' | 'Arnoldi'
m = 3;git   Np = 1;  numRequest = 2;
opts.omega = 10;  opts.T = Np*2*pi/opts.omega;
opts.c = 1;  opts.cfl = 0.9;
opts.scheme = 'Explicit';            % 'Explicit' | 'Implicit'
opts.mode = 'EigenWave';  opts.geometry = geometry;


%% 3. Set up domain, grid, boundary data
dm  = getDomain(geometry);
switch geometry
    case '1D',  par = setup1D(dm, m, opts);
    case '2D',  par = setupRectangle(dm, m, opts);
    case 'Annulus', par = setupAnnulus(dm, m, opts);
    otherwise,  error('runEigenWave: unknown geometry %s', geometry);
end
fun = getMS(dm, '', opts);                 % zero BCs, zero forcing


%% 4. Exact discrete eigenvalues (for checking)
lambda_h = par.lambdaExact;


%% 5. Run the eigensolver
counter = SolveCounter();                       % counts the wave solves
S  = @(x) filterOperator(x, par, fun, opts, counter);   % y = S_h*x (one wave solve)
ip = @(a,b) innerProduct(par.fromActive(a), par.fromActive(b), par);  % weighted (a,b)_h
x0 = par.toActive(initialGuess(par, geometry)); % initial guess (Power/Arnoldi)
tic
switch methodName
    case 'eigs'
        [Q, B] = eigs(S, par.Ngi, numRequest);
        comBeta = diag(B);
    case 'Power'
        maxit = 2000;  tolPower = 1e-8;
        [comBeta, Q, info] = powerMethod(S, x0, ip, maxit, tolPower);
        fprintf('Power method: %d iterations\n', info.numIter);
    case 'Arnoldi'
        numSteps = 80;  tolArnoldi = 1e-7;
        reorth = false;              % true -> full reorthogonalization
        [comBeta, Q, info] = arnoldi(S, x0, ip, numSteps, tolArnoldi, reorth);
        fprintf('Arnoldi: %d steps, %d Ritz pairs converged\n', info.numSteps, numel(comBeta));
end
totalTime = toc;
if isempty(comBeta)
    error('runEigenWave: no eigenvalues found');
end


%% 6. Eigenvalues of the original problem (Rayleigh quotient)
%  lambda^2 = -(v, L_h v)_h / (v, v)_h
numComputed = numel(comBeta);
comLambda = zeros(numComputed,1);
for j = 1:numComputed
    v = par.fromActive(Q(:,j));
    comLambda(j) = sqrt(-innerProduct(v, par.Lh(v), par)/normW(v,par)^2);
end


%% 7. Match with the exact discrete eigenvalues and get modified eigenvalues
disLambda = zeros(numComputed,1);
for j = 1:numComputed
    [~, index] = min(abs(comLambda(j)-lambda_h));
    disLambda(j) = lambda_h(index);
end

% Modified eigenvalues (account for time-discretization error)
dt = par.dt;
if (strcmp(opts.scheme,'Explicit'))
    modLambda = (2/dt)*asin(dt*disLambda/2);
elseif (strcmp(opts.scheme,'Implicit'))
    modLambda = (2/dt)*asin((dt*disLambda/2)./sqrt(1+(dt*disLambda).^2/2));
end
modBeta = betaDiscrete(modLambda, par, opts);


%% 8. Plots
if (opts.omega < 21)
    numPoints = 4*round(opts.omega/pi)+1;
else
    numPoints = 3*round(opts.omega/pi)+1;
end

meshPoint = linspace(0, 2 * pi*numPoints, numPoints*50+1);

% Continuous filter function
figure(1);
plot(meshPoint, betaContinuous(meshPoint,opts), '-b', 'LineWidth', 1.5); hold on;
plot(disLambda, betaContinuous(disLambda,opts), '*r', 'LineWidth', 1.2);
plot(comLambda, betaContinuous(comLambda,opts), 'ok', 'LineWidth', 1.2);
legend('\beta', 'discrete \lambda_j (true)', 'computed');
xlabel('\lambda'); ylabel('\beta'); grid on;
title(sprintf('Continuous filter function: \\omega=%g, Np=%d, %d eigs computed', ...
      opts.omega, Np, numComputed));

% Discrete filter function (evaluated at the modified eigenvalues)
figure(2);
plot(meshPoint, betaDiscrete(meshPoint,par,opts), '-b', 'LineWidth', 1.5); hold on;
plot(modLambda, modBeta, '*r', 'LineWidth', 1.2);
plot(modLambda, comBeta, 'ok', 'LineWidth', 1.2);
legend('\beta_h', 'discrete \lambda_j (modified)', 'computed');
xlabel('\lambda'); ylabel('\beta'); grid on;
title(sprintf('Discrete filter function, %d eigs computed', numComputed));


%% 9. Print results
comLambda = sort(comLambda, 'descend');
disLambda = sort(disLambda, 'descend');
errMax = abs(comLambda-disLambda);

switch geometry
    case '1D'
        fprintf('EigenWave algorithm: Nx=%d Nt=%d T=%4.2f dx=%4.2e dt=%4.2e\n', ...
                par.Nx, par.Nt, opts.T, par.dx, par.dt);
    case '2D'
        fprintf('EigenWave algorithm: Nx=Ny=%d Nt=%d T=%4.2f dx=%4.2e dy=%4.2e dt=%4.2e\n', ...
                par.Nx, par.Nt, opts.T, par.dx, par.dy, par.dt);
    case 'Annulus'
        fprintf('EigenWave algorithm: Nr=%d Nth=%d Nt=%d T=%4.2f dr=%4.2e dth=%4.2e dt=%4.2e\n', ...
                par.Nr, par.Nth, par.Nt, opts.T, par.dr, par.dth, par.dt);
end
fprintf('Domain: %s, Scheme: %s, Method: %s\n', geometry, opts.scheme, methodName);
fprintf('numComputed = %d, totalTime = %5.2f(s)\n', numComputed, totalTime);
fprintf('numWaveSolves = %d (%.1f per eigenvalue)\n', counter.n, counter.n/numComputed);
if strcmp(opts.scheme,'Implicit')
    fprintf('timeFactor = %.2e(s) (factorize the implicit matrix, size %d)\n', par.timeFactor, par.sizeA);
end
fprintf('    disLambda           comLambda          err\n');
for j = 1:numComputed
    fprintf('%15.10f     %15.10f      %5.2e\n', disLambda(j), comLambda(j), errMax(j));
end


%% 10. Compare with the reference (saved before refactoring)
refFile = fullfile(fileparts(mfilename('fullpath')), 'reference', ...
                   sprintf('%s_%s_%s.mat', geometry, opts.scheme, methodName));
if isfile(refFile)
    ref = load(refFile);
    if (ref.omega~=opts.omega || ref.m~=m || ref.Np~=Np || ref.cfl~=opts.cfl)
        fprintf('\nReference %s was computed with different parameters; skipped.\n', refFile);
    elseif (numel(ref.comLambda)~=numComputed)
        fprintf('\nReference %s has %d eigenvalues, computed %d; not compared.\n', ...
                refFile, numel(ref.comLambda), numComputed);
    else
        % eigs uses a random start vector -> compare sorted values with a tolerance
        diffCom = max(abs(sort(comLambda) - sort(ref.comLambda)));
        diffDis = max(abs(sort(disLambda) - sort(ref.disLambda)));
        tol = 1e-8;
        if (diffCom < tol && diffDis < tol), status = 'PASS'; else, status = 'FAIL'; end
        fprintf('\nReference %s:\n  max|comLambda-ref| = %.2e, max|disLambda-ref| = %.2e  %s\n', ...
                refFile, diffCom, diffDis, status);
    end
else
    fprintf('\nNo reference file %s; skipped.\n', refFile);
end


% ===== local functions =====
function u = initialGuess(par, geometry)
% Initial guess for Power/Arnoldi (same functions as getICs.m)
switch geometry
    case '1D', u = 2*par.x.*(par.x-1).^2;
    case '2D', u = cos(2*pi*par.x).*cos(pi^2*par.y);
    case 'Annulus', u = 2*par.x.*(par.x-pi).^2.*par.y.^2.*(1-par.y).^3;
end
end
