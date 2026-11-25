%% runBaselineMATLAB.m
% Baseline for anisotropicWave_2D.m: 5 cases with an exact solution,
% Nx = Ny = 20, 40, 80, 160 on (-1,1)^2, tFinal = 0.5, cfl = 0.9.
%   saveRef = true  -> write reference/matlab_<coeffs>_<ms>.mat
%   saveRef = false -> compare with the saved reference (PASS/FAIL)

clc; clear; close all;

%% 1. Settings
saveRef = false;
cases   = ["const","poly";
           "const","trig";
           "const","pulse";
           "const","trigPulse";
           "var",  "poly"];

nCases = size(cases, 1);

pass = false(nCases, 1);

%% 2. Run all cases
for k = 1:nCases
    coeffs = cases(k, 1);
    ms     = cases(k, 2);
    fprintf('\n===== case %d: coeffs=%s, ms=%s ====\n', k, coeffs, ms);
    anisotropicWave_2D;
    pass(k) = ok;
end

%% 3. Summary
fprintf('\n ==== Summary (saveRef=%d) ====\n', saveRef);
for k = 1:nCases
    if pass(k), s='PASS'; else, s='FAIL'; end
    fprintf('%d %-6s %-10s %s\n', k, cases(k,1), cases(k,2), s);
end
if all(pass), fprintf('Overall: PASS\n'); else, fprintf('Overall: FAIL\n'); end

clear cases