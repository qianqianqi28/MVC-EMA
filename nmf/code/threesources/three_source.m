clear; close all; clc;
rng(42); 
% set parameters
numDims = 3;
r = numDims;
options.maxiter = 1000;
options.display = 0;
options.model = 4;

% range
ranges = {'00', '05', '10', '15', '20', '25'};
for i = 1:length(ranges)
    range = ranges{i};
    data_filename = char("createddata\threesources\data_set_" + range + ".csv");
    data_set = csvread(data_filename);
    X = data_set(:, 2:end);

    % set lambda minvol
    options.lambda = 1;
    [Wmin, Hmin, e1min, er11min, er21min, lambdamin] = maxvolEMA(X, r, options);
    
    [~, maxPositions] = max(Wmin);
    [~, sortOrder] = sort(maxPositions);
    Hmin = Hmin(sortOrder, :);
    Wmin = Wmin(:, sortOrder);
    
    filename = ['createddata\threesources\Hmin_', range, '.csv'];
    csvwrite(filename, Hmin);
    fprintf('Saved: %s\n', filename);
    
    filename = ['createddata\threesources\Wmin_', range, '.csv'];
    csvwrite(filename, Wmin);
    fprintf('Saved: %s\n', filename);

    filename = ['createddata\threesources\lambdamin_', range, '.csv'];
    csvwrite(filename, lambdamin);
    fprintf('Saved: %s\n', filename);
    
    % set lambda novol
    options.lambda = 0;
    [Wno, Hno, e1no, er11no, er21no, lambdano] = maxvolEMA(X, r, options);
    
    [~, maxPositions] = max(Wno);
    [~, sortOrder] = sort(maxPositions);
    Hno = Hno(sortOrder, :);
    Wno = Wno(:, sortOrder);
  
    filename = ['createddata\threesources\Hno_', range, '.csv'];
    csvwrite(filename, Hno);
    fprintf('Saved: %s\n', filename);
    
    filename = ['createddata\threesources\Wno_', range, '.csv'];
    csvwrite(filename, Wno);
    fprintf('Saved: %s\n', filename);

    filename = ['createddata\threesources\lambdano_', range, '.csv'];
    csvwrite(filename, lambdano);
    fprintf('Saved: %s\n', filename);

    % set lambda maxvol
    options.lambda = -1;
    [Wmax, Hmax, e1max, er11max, er21max, lambdamax] = maxvolEMA(X, r, options);
    
    [~, maxPositions] = max(Wmax);
    [~, sortOrder] = sort(maxPositions);
    Hmax = Hmax(sortOrder, :);
    Wmax = Wmax(:, sortOrder);
    
    filename = ['createddata\threesources\Hmax_', range, '.csv'];
    csvwrite(filename, Hmax);
    fprintf('Saved: %s\n', filename);
    
    filename = ['createddata\threesources\Wmax_', range, '.csv'];
    csvwrite(filename, Wmax);
    fprintf('Saved: %s\n', filename);

    filename = ['createddata\threesources\lambdamax_', range, '.csv'];
    csvwrite(filename, lambdamax);
    fprintf('Saved: %s\n', filename);
end
