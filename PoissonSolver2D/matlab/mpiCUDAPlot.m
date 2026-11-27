
clear; clc; close all;
fontSize=16; lineWidth=2; markerSize=10;
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);


% ================== STRONG SCALING (TOTAL RUNTIME) ==================

% Process counts
np = [1 2 4 8 16];

% Grid sizes
N_vals = [256 512 1024 2048 4096];

% Total runtime (seconds)
% Rows = N, Columns = np

runtime = [
    0.126  0.532  0.689  0.721  0.866;   % N = 256
    0.281  1.063  1.258  1.281  1.464;   % N = 512
    0.817  2.065  2.473  2.617  2.962;   % N = 1024
    4.961  5.293  5.459  5.662  7.448;   % N = 2048
    23.320 22.227 18.411 17.949 17.756;  % N = 4096
];

% Plot
figure;
hold on;

for i = 1:length(N_vals)
    semilogx(np, runtime(i,:), '-o', 'LineWidth', 2, ...
        'DisplayName', sprintf('N = %d', N_vals(i)));
end

xlabel('Number of Processes (np)');
ylabel('Total Runtime (s)');
title('Runtime for (MPI + CUDA)');
legend('show', 'location', 'best');
grid on;
xticks(np);
savePlotFile("../docs/fig/mpi_cuda_runtime");


% ================== CUDA ONLY ==================

N_cuda = [256 512 1024 2048 4096];

runtime_cuda = [ ...
    0.126 ...
    0.281 ...
    0.817 ...
    4.961 ...
    23.320 ...
];
% 
% figure;
% semilogx(N_cuda, runtime_cuda, '-o','LineWidth',2);
% 
% xlabel('Grid Size N');
% ylabel('Total Runtime (s)');
% title('CUDA Only (np = 1)');
% grid on;
% xticks(N_cuda);


% Compare CUDA (np=1) vs MPI+CUDA (np=16)

runtime_np16 = runtime(:,5); % np = 16

figure;
semilogx(N_vals, runtime_cuda, '-o','LineWidth',2);
hold on;
plot(N_vals, runtime_np16, '-s','LineWidth',2);

xlabel('Grid Size N');
ylabel('Total Runtime (s)');
title('CUDA vs MPI+CUDA (np=16)');
legend('CUDA only (np=1)', 'MPI+CUDA (np=16)', 'location', 'best');
grid on;
xticks(N_cuda);
savePlotFile("../docs/fig/mpi_cuda_compare");
