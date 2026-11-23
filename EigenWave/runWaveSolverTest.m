% Convergence test for the 1D wave solver (core/waveSolve.m):
%   poly  : exact for quadratics -> maxErr ~ machine precision
%   pulse : second order         -> ratio -> 4
%   eigen : second order         -> ratio -> 4

clc; clear; close all;

%% 1. Paths
addpath(genpath(fileparts(mfilename('fullpath'))));

%% 2. Parameters
dm.LB = 0;  dm.RB = 1;              % domain [a,b]

opts.c        = 1;                  % wave speed
opts.cfl      = 0.9;                % CFL number
opts.T        = .65;                  % final time (try 0.7 for 'eigen')
opts.omega    = 10;                 % not used for this test (filter only)
opts.mode     = 'testWaveSolver';
opts.geometry = '1D';

schemes = {'Explicit','Implicit'};
tests   = {'poly','pulse','eigen'};
mList   = 1:6;                      % Nx = 20*2^m = 40,...,1280
nM      = numel(mList);

doPlot       = false;
saveBaseline = true;               % true -> reference/waveSolverBaseline.mat

%% 3. Convergence study
results = struct();
for is = 1:numel(schemes)
    opts.scheme = schemes{is};
    fprintf('\n=========== %s scheme ===========\n', opts.scheme);

    for it = 1:numel(tests)
        ms  = tests{it};
        fun = getMS(dm, ms, opts);          % exact sol, forcing, BCs, ICs

        Nx  = zeros(nM,1);
        err = zeros(nM,1);
        for k = 1:nM
            m   = mList(k);
            par = setup1D(dm, m, opts);     % grid, dt, Nt, (dA if implicit)

            % Initial conditions
            V0 = fun.u0(par.x);
            V1 = fun.u1(par.x);

            % Solve to t = T
            [~, W] = waveSolve(V0, V1, par, fun, opts);

            % Max-norm error at t = T
            err(k) = max(abs(W - fun.ue(par.x, opts.T)));
            Nx(k)  = par.Nx;

            % Same format as waveSchemes.pdf (level k = 1 <-> Nx = 20)
            fprintf('ms=%s: t=%.2e: m=%d Nx=%4d Nt=%4d dt=%.2e maxErr=%.2e', ...
                    ms, opts.T, k, par.Nx, par.Nt, par.dt, err(k));
            if k > 1
                fprintf(' ratio=%5.2f', err(k-1)/err(k));
            end
            fprintf('\n');
        end
        fprintf('\n');
        results.(opts.scheme).(ms) = struct('Nx', Nx, 'err', err);
    end
end

%% 4. Pass/fail checks
fprintf('=========== Checks ===========\n');
allPass = true;
for is = 1:numel(schemes)
    s = schemes{is};
    r = results.(s);

    % poly: exact up to round-off
    ok = max(r.poly.err) < 1e-10;
    fprintf('%-8s poly : max err = %.1e              %s\n', ...
            s, max(r.poly.err), passStr(ok));
    allPass = allPass && ok;

    % pulse/eigen: observed order on the two finest grids >= 1.8
    for ms = {'pulse','eigen'}
        e = r.(ms{1}).err;
        p = log2(e(end-1)/e(end));
        ok = p > 1.8;
        fprintf('%-8s %-5s: observed order = %4.2f        %s\n', ...
                s, ms{1}, p, passStr(ok));
        allPass = allPass && ok;
    end
end
fprintf('Overall: %s\n', passStr(allPass));

%% 5. Plots
if doPlot
    figure('Position', [100 100 1000 400]);
    for is = 1:numel(schemes)
        s = schemes{is};
        subplot(1, numel(schemes), is);
        for ms = {'pulse','eigen'}
            r = results.(s).(ms{1});
            loglog(r.Nx, r.err, 'o-', 'LineWidth', 1.5, ...
                   'DisplayName', ms{1}); hold on;
        end
        % Reference slope Nx^(-2) anchored at the pulse curve
        Nx = results.(s).pulse.Nx;
        e0 = results.(s).pulse.err(1);
        loglog(Nx, e0*(Nx/Nx(1)).^(-2), 'k--', 'DisplayName', 'O(\Delta x^2)');
        grid on; xlabel('N_x'); ylabel('max error at t = T');
        title(sprintf('%s scheme', s)); legend('Location', 'southwest');
    end
end

%% 6. Save baseline
if saveBaseline
    outFile = fullfile(fileparts(mfilename('fullpath')), ...
                       'reference', 'waveSolverBaseline.mat');
    save(outFile, 'results', 'opts', 'mList', 'dm');
    fprintf('Saved baseline to %s\n', outFile);
end

% ===== local functions =====
function s = passStr(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
