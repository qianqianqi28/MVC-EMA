clear; close all; clc;

data_set_00 = csvread('createddata\threesources\data_set_00.csv');
grain_size = data_set_00(:, 1);
n_specimens = 200;
ranges = {'00', '05', '10', '15', '20', '25'};

for i = 1:length(ranges)
    range = ranges{i};

    endmember_matrix = csvread("createddata\threesources\endmember_matrix.csv");
    abundance_matrix = csvread("createddata\threesources\abundance_matrix_"+range+".csv");
    data_set = csvread("createddata\threesources\data_set_"+range+".csv");

    lambdamin = csvread("createddata\threesources\lambdamin_"+range+".csv");
    lambdano = csvread("createddata\threesources\lambdano_"+range+".csv");
    lambdamax = csvread("createddata\threesources\lambdamax_"+range+".csv");

    Hmin = csvread("createddata\threesources\Hmin_"+range+".csv");
    Wmin = csvread("createddata\threesources\Wmin_"+range+".csv");

    Hno = csvread("createddata\threesources\Hno_"+range+".csv");
    Wno = csvread("createddata\threesources\Wno_"+range+".csv");

    Hmax = csvread("createddata\threesources\Hmax_"+range+".csv");
    Wmax = csvread("createddata\threesources\Wmax_"+range+".csv");


    load("D:/MVCEMA/nmf/createddata/threesources/EMMAresultsfromR/EMMA_mixed_"+range+"_endmem.mat");
    load("D:/MVCEMA/nmf/createddata/threesources/EMMAresultsfromR/EMMA_mixed_"+range+"_abundances.mat");
    [~, maxPositions] = max(EMMA_mixed_endmem);
    [~, sortOrder] = sort(maxPositions);
    EMMA_mixed_endmem = EMMA_mixed_endmem(:, sortOrder);
    EMMA_mixed_abundances = EMMA_mixed_abundances(sortOrder, :);


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), endmember_matrix(3, :), 'c', 'LineWidth', 2)  % Second column in blue

    % Add the x-axis label with the mu symbol
    xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Volume content", 'FontSize', 16)
    legend('True EM 1', 'True EM 2', 'True EM 3', Location="northwest", FontSize= 16)
    title("True End Members", 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/True End Members three sources", "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); 
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    plot(1:n_specimens, abundance_matrix(:, 1), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(1:n_specimens, abundance_matrix(:, 2), 'r', 'LineWidth', 2)  % Second column in blue
    plot(1:n_specimens, abundance_matrix(:, 3), 'c', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("Specimens", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Abundances", 'FontSize', 16)
    ylim([0 1])
    legend('True Abun. 1', 'True Abun. 2', 'True Abun. 3', Location="northwest", FontSize = 16)
    title("True Abundances", 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/True Abundances"+range, "-dpdf", "-bestfit");

    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), endmember_matrix(3, :), 'c', 'LineWidth', 2)  % Second column in blue
    
    plot(log(grain_size), Wmin(:,1), 'b--', 'LineWidth', 2)  % First column in red
    plot(log(grain_size), Wmin(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), Wmin(:,3), 'c--', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Volume content", 'FontSize', 16)
    legend('True EM 1', 'True EM 2', 'True EM 3', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", Location="northwest", FontSize = 16)
    lambdamin_rounded = -round(lambdamin, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdamin_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/QPthreesourceminendmem"+range, "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 2)  % Second column in blue
    axis equal
    hold on
    plot(abundance_matrix(:, 1), Hmin(1, :), 'bo', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 2), Hmin(2, :),'ro', 'LineWidth', 2)  % Second column in blue
    plot(abundance_matrix(:, 3), Hmin(3, :),'co', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("True Abundances", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Estimated Abundances", 'FontSize', 16)
    legend("1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", Location="northwest", FontSize=16)
    % title("APFGM: \lambda = "+lambdamin, 'FontSize', 16)
    lambdamin_rounded = -round(lambdamin, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdamin_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    % 使用 -bestfit 选项保存为PDF
    print("plots/QPthreesourceminabundances"+range, "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), endmember_matrix(3, :), 'c', 'LineWidth', 2)  % Second column in blue
    
    plot(log(grain_size), Wno(:,1), 'b--', 'LineWidth', 2)  % First column in red
    plot(log(grain_size), Wno(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), Wno(:,3), 'c--', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Volume content", 'FontSize', 16)
    legend('True EM 1', 'True EM 2', 'True EM 3', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", Location="northwest", FontSize = 16)
    lambdano_rounded = round(lambdano, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdano_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/QPthreesourcenoendmem"+range, "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(abundance_matrix(:, 1), abundance_matrix(:, 1),'Color', 'black', 'LineWidth', 2)  % Second column in blue
    axis equal
    hold on
    plot(abundance_matrix(:, 1), Hno(1, :), 'bo', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 2), Hno(2, :),'ro', 'LineWidth', 2)  % Second column in blue
    plot(abundance_matrix(:, 3), Hno(3, :),'co', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("True Abundances", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Estimated Abundances", 'FontSize', 16)
    legend("1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", Location="northwest", FontSize=16)
    % title("APFGM: \lambda = "+lambdano, 'FontSize', 16)
    lambdano_rounded = round(lambdano, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdano_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    % 使用 -bestfit 选项保存为PDF
    print("plots/QPthreesourcenoabundances"+range, "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), endmember_matrix(3, :), 'c', 'LineWidth', 2)  % Second column in blue
    
    plot(log(grain_size), Wmax(:,1), 'b--', 'LineWidth', 2)  % First column in red
    plot(log(grain_size), Wmax(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), Wmax(:,3), 'c--', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Volume content", 'FontSize', 16)
    legend('True EM 1', 'True EM 2', 'True EM 3', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", Location="northwest", FontSize=16)
    
    lambdamax_rounded = -round(lambdamax, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdamax_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/QPthreesourcemaxendmem"+range, "-dpdf", "-bestfit");

    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    
    plot(abundance_matrix(:, 1), abundance_matrix(:, 1), 'Color', 'black', 'LineWidth', 2)  % Second column in blue
    axis equal
    hold on
    plot(abundance_matrix(:, 1), Hmax(1, :), 'bo', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 2), Hmax(2, :), 'ro', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 3), Hmax(3, :), 'co', 'LineWidth', 2)  % First column in red
    
    % Add the x-axis label with the mu symbol
    xlabel("True Abundances", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Estimated Abundances", 'FontSize', 16)
    legend("1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", Location="northwest", FontSize = 16)
    
    lambdamax_rounded = -round(lambdamax, 3);
    title(sprintf('APFGM: \\lambda = %.3f', lambdamax_rounded), 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/QPthreesourcemaxabundances"+range, "-dpdf", "-bestfit");


    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(log(grain_size), endmember_matrix(1, :), 'b', 'LineWidth', 2)  % First column in red
    hold on
    plot(log(grain_size), endmember_matrix(2, :), 'r', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), endmember_matrix(3, :), 'c', 'LineWidth', 2)  % Second column in blue
    
    plot(log(grain_size), EMMA_mixed_endmem(:,1), 'b--', 'LineWidth', 2)  % First column in red
    plot(log(grain_size), EMMA_mixed_endmem(:,2), 'r--', 'LineWidth', 2)  % Second column in blue
    plot(log(grain_size), EMMA_mixed_endmem(:,3), 'c--', 'LineWidth', 2)  % Second column in blue
    
    % Add the x-axis label with the mu symbol
    xlabel("Ln(Grain size in \mum)", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Volume content", 'FontSize', 16)
    legend('True EM 1', 'True EM 2', 'True EM 3', "Esti. EM 1", "Esti. EM 2", "Esti. EM 3", Location="northwest", FontSize = 16)
    title("EMMA", 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/EMMAthreesourceendmem"+range, "-dpdf", "-bestfit");

    figure('Position', [100, 100, 800, 800]);
    set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
    set(gcf, 'PaperPosition', [0 0 21.0 21.0]);
    
    plot(abundance_matrix(:, 1), abundance_matrix(:, 1), 'Color', 'black', 'LineWidth', 2)  % Second column in blue
    axis equal
    hold on
    plot(abundance_matrix(:, 1), EMMA_mixed_abundances(1, :), 'bo', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 2), EMMA_mixed_abundances(2, :), 'ro', 'LineWidth', 2)  % First column in red
    plot(abundance_matrix(:, 3), EMMA_mixed_abundances(3, :), 'co', 'LineWidth', 2)  % First column in red
    
    % Add the x-axis label with the mu symbol
    xlabel("True Abundances", 'FontSize', 16)
    % Add the y-axis label
    ylabel("Estimated Abundances", 'FontSize', 16)
    legend("1:1", "Specimens: EM1", "Specimens: EM2", "Specimens: EM3", Location="northwest", FontSize = 16)
    title("EMMA", 'FontSize', 16)
    set(gca, 'FontSize', 16)
    print("plots/EMMAthreesourceabundances"+range, "-dpdf", "-bestfit");

end

