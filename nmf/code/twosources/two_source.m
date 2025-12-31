% Two-Source Mixing Simulation with Lognormal Distributions: A set of synthetic data consisting of 99 specimens "measured" with 100 grain size bins.
% Reference New methods for unmixing sediment grain size data, 2015.
clear; close all; clc;
cd('D:\MVCEMA\nmf');
restoredefaultpath
addpath(genpath('D:\MVCEMA\nmf'))

rng(42); % Set seed for reproducibility
%% Parameters
n_specimens = 99;          % Number of specimens
min_abundance = 0.13;  % Minimum abundance (13%)
max_abundance = 0.87;  % Maximum abundance (87%)
% Grain size vector
T = readtable('D:\MVCEMA\nmf\data\AnalySize-master\AnalySize-master\Example_Data\Example_Data_2.xlsx','Sheet','True_Data');
grain_size = T{:,1};
%% Define Lognormal Distributions for Sources
% Source 1
% Parameters for Source 1
sigma1 = 0.55;
mu1 = 2 + sigma1^2;  

% Parameters for Source 2
sigma2 = 1;
mu2 = 6+sigma2^2;

% Calculate lognormal probability density functions
source1 = (1 ./ (grain_size * sigma1 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu1).^2) / (2 * sigma1^2));
source1 = source1/sum(source1);
source2 = (1 ./ (grain_size * sigma2 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu2).^2) / (2 * sigma2^2));
source2 = source2/sum(source2);
endmember_matrix = zeros(2, length(source1));
endmember_matrix(1, :) = source1;
endmember_matrix(2, :) = source2;

%% Generate True Abundances
% Initialize abundance matrix
abundance_matrix = zeros(n_specimens, 2);
% Generate random abundances for Source 1 in range [0.13, 0.87]
source1_abundance = min_abundance + (max_abundance - min_abundance) * rand(n_specimens, 1);
% Set abundances for Source 2 as 1 - Source 1 abundance
source2_abundance = 1 - source1_abundance;

% Combine into abundance matrix
[~, sortIndex] = sort(source1_abundance, "descend");

abundance_matrix(:, 1) = source1_abundance(sortIndex);
abundance_matrix(:, 2) = source2_abundance(sortIndex);


%% Generate Mixed Distributions for Each Specimen
mixed_distributions = abundance_matrix * endmember_matrix;
% % The row sum of mixed_distributions
% sum(mixed_distributions, 2)
% % The size of mixed_distributions
% size(mixed_distributions)
% Save as CSV with grain size as first column
mixed_data_with_grain = [grain_size, mixed_distributions'];
csvwrite('createddata\twosources\mixed_distributions_with_grain.csv', mixed_data_with_grain);

% Note the code is for X = MH where M is basis/end-member analysis and H is
% coefficient/abundance matrix. Therefore, here we need to transform the
% matrix mixed_distributions.
% Note in the paper, before lambda, there is negative. In the algorithm,
% the sign before lambda is positive. Therefore, when we use maximum volume
% we set lambda to be negative; when we use minimum volume volume we set
% lambda to be negative. There is the sign difference from the paper.

X = mixed_distributions';

r = 2; 
options.maxiter = 1000; % In the book: 1000
options.display = 0;
% min-vol NMF (4) 
options.model=4; 

disp('Running min-vol NMF (4)...'); 
options.lambda = 500;
[Wmin,Hmin,e1min,er11min,er21min,lambdamin] = maxvolEMA(X,r,options); 


disp('Running min-vol NMF (4)...'); 
options.lambda = 0;
[Wno,Hno,e1no,er11no,er21no,lambdano] = maxvolEMA(X,r,options); 


disp('Running max-vol NMF (4)...'); 
options.lambda = -500;
[Wmax,Hmax,e1max,er11max,er21max, lambdamax] = maxvolEMA(X,r,options); 


figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
% Add the y-axis label
ylabel("Volume content", 'FontSize', 16)
legend('True EM 1', 'True EM 2', 'FontSize', 16)
title("True End Members", 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/True End Members two sources', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]);
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(1:n_specimens, abundance_matrix(:, 1), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(1:n_specimens, abundance_matrix(:, 2), 'r', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Specimens", 'FontSize', 16)
% Add the y-axis label
ylabel("Abundances", 'FontSize', 16)
ylim([0 1])
legend('True Abundances for EM 1', 'True Abundances for EM 2', 'FontSize', 16)
title("True Abundances", 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/True Abundances two sources', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
plot(log(grain_size), Wmin(:,1), 'b--', 'LineWidth', 2)  % First column in red
plot(log(grain_size), Wmin(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
% Add the y-axis label
ylabel("Volume content", 'FontSize', 16)
legend('True EM 1', 'True EM 2', "Esti. EM 1", "Esti. EM 2", 'FontSize', 16)
lambdamin_rounded = -round(lambdamin, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdamin_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourceminendmem', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 2)  % Second column in blue
axis equal
hold on
plot(abundance_matrix(:, 1), Hmin(1, :), 'bo', 'LineWidth', 2)  % First column in red
plot(abundance_matrix(:, 2), Hmin(2, :),'ro', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 16)
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 16)
legend("1:1", "Specimens: EM1", "Specimens: EM2", Location="northwest", FontSize=16)
lambdamin_rounded = -round(lambdamin, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdamin_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourceminabundances', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); 
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
plot(log(grain_size), Wno(:,1), 'b--', 'LineWidth', 2)  % First column in red
plot(log(grain_size), Wno(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
% Add the y-axis label
ylabel("Volume content", 'FontSize', 16)
legend('True EM 1', 'True EM 2', "Esti. EM 1", "Esti. EM 2", 'FontSize', 16)
lambdano_rounded = round(lambdano, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdano_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourcenoendmem', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 2)  % Second column in blue
axis equal
hold on
plot(abundance_matrix(:, 1), Hno(1, :), 'bo', 'LineWidth', 2)  % First column in red
plot(abundance_matrix(:, 2), Hno(2, :),'ro', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 16)
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 16)
legend("1:1", "Specimens: EM1", "Specimens: EM2", Location="northwest", FontSize=16)
lambdano_rounded = round(lambdano, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdano_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourcenoabundances', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
plot(log(grain_size), Wmax(:,1), 'b--', 'LineWidth', 2)  % First column in red
plot(log(grain_size), Wmax(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
% Add the y-axis label
ylabel("Volume content", 'FontSize', 16)
legend('True EM 1', 'True EM 2', "Esti. EM 1", "Esti. EM 2", 'FontSize', 16)
lambdamax_rounded = -round(lambdamax, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdamax_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourcemaxendmem', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]);
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 2)  % Second column in blue
axis equal
hold on
plot(abundance_matrix(:, 1), Hmax(1, :), 'bo', 'LineWidth', 2)  % First column in red
plot(abundance_matrix(:, 2), Hmax(2, :),'ro', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 16)
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 16)
legend("1:1", "Specimens: EM1", "Specimens: EM2", Location="northwest", FontSize = 16)
lambdamax_rounded = -round(lambdamax, 3);
title(sprintf('APFGM: \\lambda = %.3f', lambdamax_rounded), 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/QPtwosourcemaxabundances', '-dpdf', '-bestfit');

%%% EMMA %%%
% Before load the data, please generate data using
% D:\MVCEMA\lbaemma\code\twosourceEMMA.R
load("D:/MVCEMA/nmf/createddata/twosources/EMMAresultsfromR/EMMAtwosourceendmem.mat")
load("D:/MVCEMA/nmf/createddata/twosources/EMMAresultsfromR/EMMAtwosourceabundances.mat")

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
hold on
plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
plot(log(grain_size), EMMAtwosourceendmem(:,1), 'b--', 'LineWidth', 2)  % First column in red
plot(log(grain_size), EMMAtwosourceendmem(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
% Add the x-axis label with the mu symbol
xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
% Add the y-axis label
ylabel("Volume content", 'FontSize', 16)
legend('True EM 1', 'True EM 2', "Esti. EM 1", "Esti. EM 2", 'FontSize', 16)
title("EMMA", 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/EMMAtwosourceendmem', '-dpdf', '-bestfit');

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
plot(abundance_matrix(:, 1), abundance_matrix(:, 1), 'Color', 'black', 'LineWidth', 2)  % Second column in blue
axis equal
hold on
plot(abundance_matrix(:, 1), EMMAtwosourceabundances(1, :), 'bo', 'LineWidth', 2)  % First column in red
plot(abundance_matrix(:, 2), EMMAtwosourceabundances(2, :), 'ro', 'LineWidth', 2)  % First column in red
% Add the x-axis label with the mu symbol
xlabel("True Abundances", 'FontSize', 16)
% Add the y-axis label
ylabel("Estimated Abundances", 'FontSize', 16)
legend("1:1", "Specimens: EM1", "Specimens: EM2", Location="northwest", FontSize = 16)
title("EMMA", 'FontSize', 16)
set(gca, 'FontSize', 16)
print('plots/EMMAtwosourceabundances', '-dpdf', '-bestfit');

det(Wmin'*Wmin)
det(Wno'*Wno)
det(Wmax'*Wmax)
det(EMMAtwosourceendmem'*EMMAtwosourceendmem)