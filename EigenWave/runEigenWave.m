%% runEigenWave.m - compute eigenvalues of L_h * u + lambda^2 * u = 0 with EigenWave
%  Replaces main.m (EigenWave part). Results are compared to reference/*.mat.


clc; clear; close all;

%% 1. Paths
addpath(genpath(fileparts(mfilename('fullpath'))));

%% 2. Parameters  (same values as the reference runs)
geometry   = '1D';                   % later: '2D', 'Annulus'
methodName = 'eigs';                 % later: 'Power', 'Arnoldi'
m = 2;  Np = 1;  numRequest = 16;
opts.omega = 10;  opts.T = Np*2*pi/opts.omega;
opts.c = 1;  opts.cfl = 0.6;
opts.scheme = 'Implicit';            % 'Explicit' | 'Implicit'
opts.mode = 'EigenWave';  opts.geometry = geometry;


%% 3. Set up domain, grid, boundary data
dm  = getDomain(geometry);
par = setup1D(dm, m, opts);                % later: switch on geometry
fun = getMS(dm, '', opts);                 % zero BCs, zero forcing


%% 4. Exact discrete eigenvalues (for checking)
lambda_h = par.lambdaExact;


%% 5. Run the eigensolver
S = @(x) filterOperator(x, par, fun, opts);    % y = S_h*x (one wave solve)
tic
switch methodName
    case 'eigs'
        [Q, B] = eigs(S, par.Ngi, numRequest);
        comBeta = diag(B);
    case 'Power'
        % later: x0 = par.toActive(2*par.x.*(par.x-1).^2);  (from getICs)
        %        [comBeta, Q] = powerMethod(S, x0, par, ...);   (solvers/)
        error('runEigenWave: Power method not implemented yet');
    case 'Arnoldi'
        % later: [comBeta, Q] = arnoldi(S, x0, par, ...);       (solvers/)
        error('runEigenWave: Arnoldi method not implemented yet');
end
totalTime = toc;


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

meshPoint = linspace(0, pi*numPoints, numPoints*50+1);

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

fprintf('EigenWave algorithm: Nx=%d Nt=%d T=%4.2f dx=%4.2e dt=%4.2e\n', ...
        par.Nx, par.Nt, opts.T, par.dx, par.dt);
fprintf('Domain: %s, Scheme: %s, Method: %s\n', geometry, opts.scheme, methodName);
fprintf('numComputed = %d, totalTime = %5.2f(s)\n', numComputed, totalTime);
fprintf('    disLambda           comLambda          err\n');
for j = 1:numComputed
    fprintf('%15.10f     %15.10f      %5.2e\n', disLambda(j), comLambda(j), errMax(j));
end


%% 10. Compare with the reference (saved before refactoring)
refFile = fullfile(fileparts(mfilename('fullpath')), 'reference', ...
                   sprintf('%s_%s_%s.mat', geometry, opts.scheme, methodName));
if isfile(refFile)
    ref = load(refFile);
    if (ref.omega==opts.omega && ref.m==m && ref.Np==Np && ref.cfl==opts.cfl ...
            && numel(ref.comLambda)==numComputed)
        % eigs uses a random start vector -> compare sorted values with a tolerance
        diffCom = max(abs(sort(comLambda) - sort(ref.comLambda)));
        diffDis = max(abs(sort(disLambda) - sort(ref.disLambda)));
        tol = 1e-8;
        if (diffCom < tol && diffDis < tol), status = 'PASS'; else, status = 'FAIL'; end
        fprintf('\nReference %s:\n  max|comLambda-ref| = %.2e, max|disLambda-ref| = %.2e  %s\n', ...
                refFile, diffCom, diffDis, status);
    else
        fprintf('\nReference %s was computed with different parameters; skipped.\n', refFile);
    end
else
    fprintf('\nNo reference file %s; skipped.\n', refFile);
end
