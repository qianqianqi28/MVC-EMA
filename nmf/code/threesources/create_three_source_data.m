% Create data: Three-Source Mixing Simulation with Lognormal Distributions by fixing peak location.
% Reference New methods for unmixing sediment grain size data, 2015.
% For lognormal distribution, the distribution is f(x) = (1/(xσ√(2π))) * exp(-(ln(x)-μ)²/(2σ²)) for x > 0; 
% the location peak is in exp(mu-delta^2);
% peak_height = exp(-(mu - sigma^2/2)) / (sigma * sqrt(2*pi));
clear; close all; clc;
cd('D:\MVCEMA\nmf');
restoredefaultpath
addpath(genpath('D:\MVCEMA\nmf'))
rng(42); % Set seed for reproducibility
%% Parameters
% Grain size vector
T = readtable('D:\MVCEMA\nmf\data\AnalySize-master\AnalySize-master\Example_Data\Example_Data_2.xlsx','Sheet','True_Data');
grain_size = T{:,1};

%% Define Lognormal Distributions for Sources (Figure 1d): 

sigma1 = 0.55;    % Standard deviation of log(grain_size)
mu1 = 2 + sigma1^2; % when we plot Figure 1d where log(x) as peak location, we hope the peak location is 2. Thus, peak location is log(x) = 2, i.e., x = exp(2). Therefore, mu = 2+delta^2
sigma2 = 0.55;      % Standard deviation of log(grain_size)
mu2 = 4+sigma2^2;         % Mean of log(grain_size)
sigma3 = 0.55;      % Standard deviation of log(grain_size)
mu3 = 6+sigma3^2;         % Mean of log(grain_size)

% Calculate lognormal probability density functions
source1 = (1 ./ (grain_size * sigma1 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu1).^2) / (2 * sigma1^2));
source1 = source1/sum(source1);
source2 = (1 ./ (grain_size * sigma2 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu2).^2) / (2 * sigma2^2));
source2 = source2/sum(source2);
source3 = (1 ./ (grain_size * sigma3 * sqrt(2 * pi))) .* ...
          exp(-((log(grain_size) - mu3).^2) / (2 * sigma3^2));
source3 = source3/sum(source3);
endmember_matrix = zeros(3, length(source1));
endmember_matrix(1, :) = source1;
endmember_matrix(2, :) = source2;
endmember_matrix(3, :) = source3;

% Sort rows by the position of the maximum value in each row.
[~, maxPositions] = max(endmember_matrix, [], 2);
[~, sortOrder] = sort(maxPositions);
endmember_matrix = endmember_matrix(sortOrder, :);

filename = char("createddata\threesources\endmember_matrix.csv");
csvwrite(filename, endmember_matrix);

%% Generate True Abundances (Figure 3b)
% Initialize abundance matrix
n_specimens = 200;          % Number of specimens
numDims = 3;
% Generate random matrix in [0,1]
matrix = rand(n_specimens, numDims);
% Normalize each row to sum to 1
rowSums = sum(matrix, 2);
normalizedMatrix = matrix ./ rowSums;


% Define the abundance ranges
abundance_ranges = {
    0.00, 1.00, "00";
    0.05, 0.95, "05";
    0.10, 0.90, "10";
    0.15, 0.85, "15";
    0.20, 0.80, "20";
    0.25, 0.75, "25"
};

% Loop through each abundance range
for i = 1:size(abundance_ranges, 1)
    min_abundance = abundance_ranges{i, 1};  % Minimum abundance
    max_abundance = abundance_ranges{i, 2};  % Maximum abundance
    range = abundance_ranges{i, 3};          % Range identifier
    
    fprintf('Processing range %s: min=%.2f, max=%.2f\n', min_abundance, max_abundance, range);
    
    % Generate the data
    [abundance_matrix, mixed_data_with_grain] = simulate_Mixed_abundance_data(...
        normalizedMatrix, endmember_matrix, grain_size, min_abundance, max_abundance);
    
    % Save abundance matrix
    filename = char("createddata\threesources\abundance_matrix_" + range + ".csv");
    csvwrite(filename, abundance_matrix);
    fprintf('Saved: %s\n', filename);
    
    % Save mixed data with grain
    filename = char("createddata\threesources\data_set_" + range + ".csv");
    csvwrite(filename, mixed_data_with_grain);
    fprintf('Saved: %s\n', filename);
end
