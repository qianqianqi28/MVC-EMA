clear; close all; clc;

n_specimens = 200;
ranges = {'00', '05', '10', '15', '20', '25'};
mixed_range = 0:5:25;

EMMA_EMs_Angle = NaN(length(ranges), 1);
EMMA_Spec_Angle = NaN(length(ranges), 1);
EMMA_DataSet_det = NaN(length(ranges), 1);
EMMA_DataSet_data_fit = NaN(length(ranges), 1);

QP_EMs_Angle_min = NaN(length(ranges), 1);
QP_EMs_Angle_no = NaN(length(ranges), 1);
QP_EMs_Angle_max = NaN(length(ranges), 1);

QP_Spec_Angle_min = NaN(length(ranges), 1);
QP_Spec_Angle_no = NaN(length(ranges), 1);
QP_Spec_Angle_max = NaN(length(ranges), 1);

QP_DataSet_det_min = NaN(length(ranges), 1);
QP_DataSet_det_no = NaN(length(ranges), 1);
QP_DataSet_det_max = NaN(length(ranges), 1);

QP_DataSet_data_fit_min = NaN(length(ranges), 1);
QP_DataSet_data_fit_no = NaN(length(ranges), 1);
QP_DataSet_data_fit_max = NaN(length(ranges), 1);

endmember_matrix = csvread("createddata\threesources\endmember_matrix.csv");

for i = 1:length(ranges)
    range = ranges{i};
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

    V1 = abundance_matrix;
    V2 = Hmin';
    QP_Spec_Angle_min(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));

    V1 = endmember_matrix;
    V2 = Wmin';
    QP_EMs_Angle_min(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));


    V1 = data_set(:, 2:end)';
    V2 = Wmin * Hmin;
    V2 = V2';
    QP_DataSet_det_min(i, 1) = det(Wmin' * Wmin);
    QP_DataSet_data_fit_min(i, 1) = norm(V1 - V2, 'fro');

    V1 = abundance_matrix;
    V2 = Hno';
    QP_Spec_Angle_no(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));

    V1 = endmember_matrix;
    V2 = Wno';
    QP_EMs_Angle_no(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));

    V1 = data_set(:, 2:end)';
    V2 = Wno * Hno;
    V2 = V2';
    QP_DataSet_det_no(i, 1) = det(Wno' * Wno);
    QP_DataSet_data_fit_no(i, 1) = norm(V1 - V2, 'fro');


    V1 = endmember_matrix;
    V2 = Wmax';
    QP_EMs_Angle_max(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));
    
    
    V1 = abundance_matrix;
    V2 = Hmax';
    QP_Spec_Angle_max(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));

    V1 = data_set(:, 2:end)';
    V2 = Wmax * Hmax;
    V2 = V2';
    QP_DataSet_det_max(i, 1) = det(Wmax' * Wmax);
    QP_DataSet_data_fit_max(i, 1) = norm(V1 - V2, 'fro');

    V1 = endmember_matrix;
    V2 = EMMA_mixed_endmem';
    EMMA_EMs_Angle(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));
    
    V1 = abundance_matrix;
    V2 = EMMA_mixed_abundances';
    EMMA_Spec_Angle(i, 1) = mean(rad2deg(acos((sum(V1.*V2,2)./sqrt(sum(V1.^2,2).*sum((V2).^2, 2))))));


    V1 = data_set(:, 2:end)';
    V2 = (EMMA_mixed_abundances') * (EMMA_mixed_endmem');
    EMMA_DataSet_det(i, 1) = det(EMMA_mixed_endmem' * EMMA_mixed_endmem);
    EMMA_DataSet_data_fit(i, 1) = norm(V1 - V2, 'fro');

end


figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);

plot(mixed_range, EMMA_EMs_Angle, 'bo-', 'LineWidth', 2)  % First column in red
hold on
plot(mixed_range, QP_EMs_Angle_min, 'ro-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_EMs_Angle_no, 'go-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_EMs_Angle_max, 'ko-', 'LineWidth', 2)  % Second column in blue

% Add the x-axis label with the mu symbol
xlabel("Minimum abundance [%]", 'FontSize', 16)
% Add the y-axis label
ylabel("MAEM", 'FontSize', 16)
legend('EMMA', "APFGM: Minimum volume", "APFGM: No volume", "APFGM: Maximum volume", 'Location', 'northwest', 'FontSize', 16)
title("MAEM", 'FontSize', 16)
set(gca, 'FontSize', 16)
print("plots/End Member Misfit", "-dpdf", "-bestfit");


figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);

plot(mixed_range, EMMA_Spec_Angle, 'bo-', 'LineWidth', 2)  % First column in red
hold on
plot(mixed_range, QP_Spec_Angle_min, 'ro-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_Spec_Angle_no, 'go-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_Spec_Angle_max, 'ko-', 'LineWidth', 2)  % Second column in blue

% Add the x-axis label with the mu symbol
xlabel("Minimum abundance [%]", 'FontSize', 16)
% Add the y-axis label
ylabel("MAAB", 'FontSize', 16)
legend('EMMA', "APFGM: Minimum volume", "APFGM: No volume", "APFGM: Maximum volume", 'Location', 'northwest', 'FontSize', 16)
title("MAAB", 'FontSize', 16)
set(gca, 'FontSize', 16)
print("plots/Abundance Misfit", "-dpdf", "-bestfit");


figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);

plot(mixed_range, EMMA_DataSet_det, 'bo-', 'LineWidth', 2)  % First column in red
hold on
plot(mixed_range, QP_DataSet_det_min, 'ro-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_DataSet_det_no, 'go-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_DataSet_det_max, 'ko-', 'LineWidth', 2)  % Second column in blue

% Add the x-axis label with the mu symbol
xlabel("Minimum abundance [%]", 'FontSize', 16)
% Add the y-axis label
ylabel("Determinant", 'FontSize', 16)
legend('EMMA', "APFGM: Minimum volume", "APFGM: No volume", "APFGM: Maximum volume", 'Location', 'east', 'FontSize', 16)
title("Determinant", 'FontSize', 16)
set(gca, 'FontSize', 16)
print("plots/Determinant for Basis Matrix", "-dpdf", "-bestfit");

figure('Position', [100, 100, 800, 800]);
set(gcf, 'PaperSize', [21.0 21.0]); % A4 landscape in centimeters
set(gcf, 'PaperPosition', [0 0 21.0 21.0]);

plot(mixed_range, EMMA_DataSet_data_fit, 'bo-', 'LineWidth', 2)  % First column in red
hold on
plot(mixed_range, QP_DataSet_data_fit_min, 'ro-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_DataSet_data_fit_no, 'go-', 'LineWidth', 2)  % Second column in blue
plot(mixed_range, QP_DataSet_data_fit_max, 'ko-', 'LineWidth', 2)  % Second column in blue

% Add the x-axis label with the mu symbol
xlabel("Minimum abundance [%]", 'FontSize', 16)
% Add the y-axis label
ylabel("Data Fit", 'FontSize', 16)
legend('EMMA', "APFGM: Minimum volume", "APFGM: No volume", "APFGM: Maximum volume", 'Location', 'northwest', 'FontSize', 16)
title("Data Fit for Data Set", 'FontSize', 16)
set(gca, 'FontSize', 16)
print("plots/Data Fit for Data Set", "-dpdf", "-bestfit");

MAEM_apfgd = [EMMA_EMs_Angle, QP_EMs_Angle_min, QP_EMs_Angle_no, QP_EMs_Angle_max];
disp(round(MAEM_apfgd, 4))

MAAB_apfgd = [EMMA_Spec_Angle, QP_Spec_Angle_min, QP_Spec_Angle_no, QP_Spec_Angle_max];
disp(round(MAAB_apfgd, 4))

det_apfgd = [EMMA_DataSet_det, QP_DataSet_det_min, QP_DataSet_det_no, QP_DataSet_det_max];
disp(round(det_apfgd * 1e3, 4))