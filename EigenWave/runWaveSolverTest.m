% Convergence test for the wave solver (core/waveSolve.m) in 1D and 2D:
%   poly  : exact for quadratics -> maxErr ~ machine precision
%   1D pulse, eigen        : second order -> ratio -> 4
%   2D trig, eigen         : second order -> ratio -> 4
%   1D tables are in waveSchemes.pdf.

clc; clear; close all;

%% 1. Paths
addpath(genpath(fileparts(mfilename('fullpath'))));

%% 2. Parameters
geometry = '2D';                    % '1D' | '2D'  (later: 'Annulus')

opts.c        = 1;                  % wave speed
opts.cfl      = 0.9;                % CFL number
opts.T        = .65;                % final time
opts.omega    = 10;                 % not used for this test (filter only)
opts.mode     = 'testWaveSolver';
opts.geometry = geometry;

dm = getDomain(geometry);           % 1D: [0,1], 2D: [0,1]x[0,1]
schemes = {'Explicit','Implicit'};
switch geometry
    case '1D'
        tests = {'poly','pulse','eigen'};
        mList = 1:6;                % Nx = 20*2^m     = 40,...,1280
        setupGrid = @setup1D;
    case '2D'
        tests = {'poly','trig','eigen'};
        mList = 1:5;                % Nx = Ny = 20*2^(m-1) = 20,...,320
        setupGrid = @setupRectangle;
    otherwise
        error('runWaveSolverTest: geometry %s not implemented yet', geometry);
end
nM = numel(mList);

doPlot       = false;
saveBaseline = false;               % true -> reference/waveSolverBaseline[_2D].mat

%% 3. Convergence study
results = struct();
for is = 1:numel(schemes)
    opts.scheme = schemes{is};
    fprintf('\n=========== %s: %s scheme ===========\n', geometry, opts.scheme);

    for it = 1:numel(tests)
        ms  = tests{it};
        fun = getMS(dm, ms, opts);          % exact sol, forcing, BCs, ICs

        Nx  = zeros(nM,1);
        err = zeros(nM,1);
        for k = 1:nM
            m   = mList(k);
            par = setupGrid(dm, m, opts);   % grid, dt, Nt, (dA if implicit)

            % Initial conditions
            V0 = gridEval(fun.u0, par, geometry);
            V1 = gridEval(fun.u1, par, geometry);

            % Solve to t = T
            [~, W] = waveSolve(V0, V1, par, fun, opts);

            % Max-norm error at t = T
            uexact = gridEval(@(varargin) fun.ue(varargin{:}, opts.T), par, geometry);
            err(k) = max(abs(W - uexact), [], 'all');
            Nx(k)  = par.Nx;

            % Same format as waveSchemes.pdf
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
fprintf('=========== Checks (%s) ===========\n', geometry);
allPass = true;
convTests = setdiff(tests, {'poly'}, 'stable');   % tests that should converge at 2nd order
for is = 1:numel(schemes)
    s = schemes{is};
    r = results.(s);

    % poly: exact up to round-off
    ok = max(r.poly.err) < 1e-10;
    fprintf('%-8s poly : max err = %.1e              %s\n', ...
            s, max(r.poly.err), passStr(ok));
    allPass = allPass && ok;

    % other tests: observed order on the two finest grids >= 1.8
    for ms = convTests
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
        for ms = convTests
            r = results.(s).(ms{1});
            loglog(r.Nx, r.err, 'o-', 'LineWidth', 1.5, ...
                   'DisplayName', ms{1}); hold on;
        end
        % Reference slope Nx^(-2) anchored at the first convergence test
        r  = results.(s).(convTests{1});
        loglog(r.Nx, r.err(1)*(r.Nx/r.Nx(1)).^(-2), 'k--', 'DisplayName', 'O(\Delta x^2)');
        grid on; xlabel('N_x'); ylabel('max error at t = T');
        title(sprintf('%s: %s scheme', geometry, s)); legend('Location', 'southwest');
    end
end

%% 6. Save baseline
if saveBaseline
    if strcmp(geometry,'1D')
        fname = 'waveSolverBaseline.mat';            % keep the existing 1D name
    else
        fname = sprintf('waveSolverBaseline_%s.mat', geometry);
    end
    outFile = fullfile(fileparts(mfilename('fullpath')), 'reference', fname);
    save(outFile, 'results', 'opts', 'mList', 'dm');
    fprintf('Saved baseline to %s\n', outFile);
end

% ===== local functions =====
function s = passStr(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function u = gridEval(f, par, geometry)
% Evaluate f on the grid: f(x) in 1D, f(x,y) in 2D
switch geometry
    case '1D', u = f(par.x);
    case '2D', u = f(par.x, par.y);
end
end
