close all
clc
clear
directory = 'F:\Olivier';
filename = fullfile(directory, 'dataset_thesis_complete.csv');
df = readtable(filename);
df = df(df.survival_time < 1280, :);
Npat = height(df);
N = height(df);
lastMRI = df.LAST_MR/30; % True survival days
% load('unpdf.mat')
% load('pdf_sim.mat')
load('unpdf_KDE_thesis.mat')
unpdf=unpdfLOOCV;
load('pdf_sim_Norm.mat')
% pdf = pdf;

survival_time = df.survival_time/30;
BaM3 = cell(1, N);
pdf_disc = 100;  % Number of discretization steps for the PDF
target_range = linspace(0, 40, pdf_disc);
for i = 1:N
    % Calculate the corrected probability distribution and normalize it
    temp3 = unpdf{i} .* pdf{i}';
    BaM3{i} = temp3 ./ trapz(target_range, temp3);
end

for i = 1:N
    % For pdf
    pdf_means(i) = trapz(target_range, target_range .* pdf{i});
    [~, idx] = max(pdf{i});
    pdf_modes(i) = target_range(idx);
    cdf = cumtrapz(target_range, pdf{i});
    [~, idx] = min(abs(cdf - 0.5));
    pdf_medians(i) = target_range(idx);
    
%     % For unpdf
    unpdf_means(i) = trapz(target_range, target_range .* unpdf{i}');
    [~, idx] = max(unpdf{i});
    unpdf_modes(i) = target_range(idx);
    cdf = cumtrapz(target_range, unpdf{i}');
    [~, idx] = min(abs(cdf - 0.5));
    unpdf_medians(i) = target_range(idx);
%     
%     % For BaM3
    BaM3_means(i) = trapz(target_range, target_range .* BaM3{i}');
    [~, idx] = max(BaM3{i});
    BaM3_modes(i) = target_range(idx);
    cdf = cumtrapz(target_range, BaM3{i}');
    [~, idx] = min(abs(cdf - 0.5));
    BaM3_medians(i) = target_range(idx);
end

% %% Calculate the error for the expected value, mode and median
MSE_mod = 1/N * ( sum( (lastMRI' - pdf_modes).^2 ) );
MSE_bam3 = 1/N * ( sum( (lastMRI' - BaM3_modes).^2 ) );
MSE_unmod = 1/N * ( sum( (lastMRI' - unpdf_modes).^2 ) );

fprintf('\n\nMSE_mod = %.3f',MSE_mod);
fprintf('\n\nMSE_bam = %.3f ', MSE_bam3);
fprintf('\n\nMSE_unmod = %.3f\n', MSE_unmod);


% Enhanced visualization for Medians Estimation
figure;
scatter(lastMRI, pdf_modes, 100, 'r.','linewidth', 3,'DisplayName', 'Modelable'); hold on;
scatter(lastMRI, unpdf_modes, 100, 'g+','linewidth', 2, 'DisplayName', 'Unmodelable');
scatter(lastMRI, BaM3_modes, 100, 'bx','linewidth', 2, 'DisplayName', 'BaM3');
% Plotting the perfect fit line
max_val = max([max(lastMRI), max(pdf_medians), max(unpdf_medians), max(BaM3_medians)]);
min_val = min([min(lastMRI), min(pdf_medians), min(unpdf_medians), min(BaM3_medians)]);
plot([min_val, 30], [min_val, 30], 'k--', 'DisplayName', 'Perfect Fit');
axis square; % To make x and y axis equal for better visual comparison
xlabel('Actual Time to Relapse', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Estimated Time to Relapse', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12, 'FontWeight', 'bold');
%title('Mod Estimation', 'FontSize', 16, 'FontWeight', 'bold');
grid on;  % Optional: Adding a grid for better visualization
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
figure;
numPatientsToPlot = 16; % Number of patients to visualize
for i = 1:numPatientsToPlot
    subplot(4, 4, i);  % Arrange plots in 4 rows and 4 columns 
    % Plot each PDF
    plot(target_range, pdf{i}, 'b', 'LineWidth', 2); hold on;
    plot(target_range, unpdf{i}, 'r', 'LineWidth', 2);
    plot(target_range, BaM3{i}, 'g', 'LineWidth', 2); hold off;
    line([lastMRI(i), lastMRI(i)], [0, 0.30], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
    line([survival_time(i), survival_time(i)], [0, 0.30], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--');
    title(['Patient ' num2str(i)], 'FontSize', 12);
    xlabel('Growth Time (in months)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
    if i == 1  % Only add legend to the first subplot
        legend('PDF', 'UnPDF', 'BaM3', 'lastMRI', 'Survival', 'Location', 'best');
    end
    grid on; % Optional: Add grid for clarity

    % Set font size for axis numbers
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
end
% 
% sgtitle('Probability Density Functions for the Lux Patients');  % Super title for the entire figure
