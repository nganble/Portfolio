clc;
clear; close all;

fontSize=16; lineWidth=2; markerSize=10;
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);
% ================== WEAK SCALING (TIME PER ITER) ==================

% Data
np_weak = [2 8 32 128];
N_weak  = [256 512 1024 2048];
tpi_weak = [0.6148 0.7022 0.8214 0.9889]; % ms

% Sort (safety)
[np_weak, idx] = sort(np_weak);
tpi_weak = tpi_weak(idx);
N_weak = N_weak(idx);

% Ideal: constant time per iteration
ideal_tpi = tpi_weak(1) * ones(size(tpi_weak));

% --------- PLOT ---------
figure;
% semilogx(np_weak, tpi_weak);
semilogx(np_weak, tpi_weak, '--o');
% hold on;
% semilogx(np_weak, ideal_tpi, '--','LineWidth',1.5);

xticks(np_weak);
xticklabels(string(np_weak));

xlabel('Number of Processes (np)');
ylabel('Time per Iteration (ms)');
title('Weak Scaling (Time per Iteration)');
% legend('Measured','Ideal (constant)','Location','northwest');
ylim([0 2]);
grid on;
savePlotFile("../docs/fig/weak_scaling_runtime");

% --------- Efficiency (optional) ---------
eff_weak = tpi_weak(1) ./ tpi_weak;

figure;
semilogx(np_weak, eff_weak, '--o','LineWidth',2);

xticks(np_weak);
xticklabels(string(np_weak));

xlabel('Number of Processes (np)');
ylabel('Weak Scaling Efficiency');
title('Weak Scaling Efficiency (Time per Iterarion)');
grid on;
savePlotFile("../docs/fig/weak_scaling_efficiency");

% --------- PRINT TABLE ---------
disp('--------------------------------------');
disp('Weak Scaling (Time per Iter)');
disp('np        N        tpi(ms)     efficiency');
disp([np_weak' N_weak' tpi_weak' eff_weak']);

