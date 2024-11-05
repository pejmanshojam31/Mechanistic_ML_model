close all
clc
clear
% Load the dataset
directory = 'F:\Olivier';
filename = fullfile(directory, 'dataset_thesis_complete.csv');
df = readtable(filename);
df = df(df.survival_time < 1280, :);
%Survival_days = df.survival_time / 30;
Survival_days= df.LAST_MR / 30; % True survival days

% Parameters
Npat = height(df);
% K = 4.287e10;
K= 100;
lamda = 1;
% Initial conditions and time span
% initial_conditions = 1.7511e+09;
initial_conditions = 3.19;
t = linspace(0, 1300, 1000);
pmvar = 0.4;
num_samples = 100;
rng(123)

r1_mean = zeros(1, Npat); 
time_results = NaN(Npat, 1);
time_results_rmin = NaN(Npat, 1);
time_results_rmax = NaN(Npat, 1);
time_results_avg = NaN(Npat, 1);
patient_std = NaN(Npat, 1);
pdf_disc = 100;
target_range = linspace(0, 40, pdf_disc);
pdf = cell(1, Npat);

for ipat = 1:Npat % Loop over the patient ensemble
    current_time = Survival_days(ipat);
    r1_mean(ipat) = df.growth(ipat);

    rmin = r1_mean(ipat) * (1 - pmvar);
    rmax = r1_mean(ipat) * (1 + pmvar);
    
    % Find results for rmin
    [T_rmin, Y_rmin] = ode45(@(t, y) tumor_system(t, y, rmin, lamda, K), t, initial_conditions);
    idx_rmin = find(Y_rmin(:, 1) >= 0.80 * K, 1, 'first');
    if ~isempty(idx_rmin)
        time_results_rmin(ipat) = T_rmin(idx_rmin) / 30;
    else
        time_results_rmin(ipat) = 24; % Set to 24 if the condition isn't met
    end
    
    % Find results for rmax
    [T_rmax, Y_rmax] = ode45(@(t, y) tumor_system(t, y, rmax, lamda, K), t, initial_conditions);
    idx_rmax = find(Y_rmax(:, 1) >= 0.80 * K, 1, 'first');
    if ~isempty(idx_rmax)
        time_results_rmax(ipat) = T_rmax(idx_rmax) / 30;
    else
        time_results_rmax(ipat) = 24; % Set to 24 if the condition isn't met
    end
    
    % Find results for r1_mean
    [T, Y] = ode45(@(t, y) tumor_system(t, y, r1_mean(ipat), lamda, K), t, initial_conditions);
    idx = find(Y(:, 1) >= 0.99 * K, 1, 'first');
    if ~isempty(idx)
        time_results(ipat) = T(idx) / 30;
    else
        time_results(ipat) = 15; % Set to 15 if the condition isn't met
    end
    
    IQR = time_results_rmin(ipat) + time_results_rmax(ipat);
    time_results_avg(ipat) = IQR / 2;

    % Compute the standard deviation for the current patient
    patient_std(ipat) =  IQR / 1.349;
end

pdf_means = NaN(1, Npat);
pdf_modes = NaN(1, Npat);
pdf_medians = NaN(1, Npat);
bin_centers = target_range;

for ipat = 1:Npat
    pdf{ipat} = normpdf(target_range, time_results_avg(ipat), patient_std(ipat));
    % Calculations for pdf mean, mode, median
    pdf_means(ipat) = trapz(bin_centers, bin_centers .* pdf{ipat});
    [~, idx] = max(pdf{ipat});
    pdf_modes(ipat) = bin_centers(idx);
    cdf = cumtrapz(bin_centers, pdf{ipat});
    [~, idx] = min(abs(cdf - 0.5));
    pdf_medians(ipat) = bin_centers(idx);
end

% Calculate the error for the expected value, mode, and median
MSE_mod = 1 / Npat * sum((Survival_days' - pdf_modes).^2);
fprintf('\n\nMSE_mod = %.3f\n', MSE_mod);

for i = 1:16
    subplot(4, 4, i); % Arrange plots in 5 rows and 4 columns
    % Plot each PDF
    plot(bin_centers, pdf{i}, 'b', 'LineWidth', 1); hold on;
    line([Survival_days(i), Survival_days(i)], [0, 0.45], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--')
    title(['Patient ' num2str(i)]);
    xlabel('Growth Time (in months)');
    ylabel('Probability Density');
    grid on; % Optional: Add grid for clarity
end

figure;
scatter(Survival_days, pdf_modes, 'filled');
xlabel('True Survival Time (months)');
ylabel('Predicted Survival Time (months)');
title('Validation: True vs. Predicted Survival Time');
axis equal;
grid on;
hold on;
x = linspace(0, 35, 100);
plot(x, x, 'r--');
hold off;

% Box Plot
figure;
boxplot([time_results, time_results_avg]);
title('Box Plot of time\_results and time\_results\_avg');
set(gca, 'XTickLabel', {'time\_results', 'time\_results\_avg'});
ylabel('Values');
grid on;

% Bar Plot
figure;
bar(1:length(time_results), [time_results, time_results_avg]);
title('Bar Plot of time\_results and time\_results\_avg');
xlabel('Patient');
ylabel('Values');
legend('time\_results', 'time\_results\_avg');
grid on;

save('pdf_sim_Norm.mat', 'pdf_means', 'pdf_modes', 'pdf_medians', 'pdf')

function dy = tumor_system(~, y, r1, lamda, K)
    c = y(1);
    dcdt = r1 * c^lamda * (1 - (c / K));
    dy = dcdt;
end
