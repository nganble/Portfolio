fontSize=16; lineWidth=2; markerSize=10;
set(0,'DefaultLineMarkerSize',markerSize);
set(0,'DefaultLineLineWidth',lineWidth);
set(0,'DefaultAxesFontSize',fontSize);

% ================== STRONG SCALING PLOT SCRIPT ==================
clear; clc; close all;

% --------- USER CHOICE ---------
N = 128;

% --------- DATA ---------
switch N
    case 128
        np   = [128 64 32 16 8 4 2];
        time = [0.0107 0.0131 0.0149 0.0155 0.0352 0.0378 0.0678];
        
    case 256
        np   = [128 64 32 16 8 4 2];
        time = [0.0290 0.0450 0.0688 0.0966 0.1620 0.2818 0.5431];
        
    case 512
        np   = [128 64 32 16 8 4 2];
        time = [0.1556 0.2274 0.3881 0.6973 1.3198 2.2888 4.1806];
        
    case 1024
        np   = [128 64 32 16 8 4 2];
        time = [1.1979 1.8131 2.8987 4.8855 9.2304 17.8128 33.1049];
        
    otherwise
        error('Invalid N value. Choose 128, 256, 512, or 1024.');
end

% --------- SORT (ascending np for nicer plots) ---------
[np, idx] = sort(np);
time = time(idx);

% --------- COMPUTE SPEEDUP & EFFICIENCY ---------
T1 = time(1); % smallest np (baseline)
speedup = T1 ./ time;
efficiency = speedup ./ np;

% --------- PLOT 1: TIME vs NP ---------
figure;
semilogx(np, time, '-o','LineWidth',2);
xlabel('Number of Processes (np)');
ylabel('Time (s)');
xticks(np);
title(['Strong Scaling: Time vs np (N = ', num2str(N), ')']);
grid on;
savePlotFile(sprintf("../docs/fig/strong_scaling_time_N_%d",N));

% --------- PLOT 2: SPEEDUP ---------
figure;
loglog(np, speedup, '-o','LineWidth',2);
hold on;
loglog(np, np, '--','LineWidth', 1.5); % ideal
xlabel('Number of Processes (np)');
ylabel('Speedup');
xticks(np);
title(['Strong Scaling: Speedup (N = ', num2str(N), ')']);
legend('Measured','Ideal','Location','northwest');
grid on;
savePlotFile(sprintf("../docs/fig/strong_scaling_speedup_N_%d",N));

% % --------- PLOT 3: EFFICIENCY ---------
% figure;
% semilogx(np, efficiency, '-o','LineWidth',2);
% xlabel('Number of Processes (np)');
% ylabel('Efficiency');
% title(['Strong Scaling: Efficiency (N = ', num2str(N), ')']);
% xticks(np);
% grid on;

% --------- OPTIONAL: PRINT TABLE ---------
% disp('--------------------------------------');
% disp(['Results for N = ', num2str(N)]);
% disp('np        time(s)     speedup     efficiency');
% disp([np' time' speedup' efficiency']);



