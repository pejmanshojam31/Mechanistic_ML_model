close all
clc
clear
directory = 'F:\Olivier';
filename = fullfile(directory, 'dataset_thesis_complete.csv');
df = readtable(filename);
df = df(df.survival_time <= 1280, :);
Time_to_progression = df.LAST_MR;
Time_to_progression = Time_to_progression/30;
r = df.growth;
threshold = 10;%mean(Time_to_progression);
pmvar = 0.4;
num_samples = 500;
sigma = 0.01;
pdf_disc = 100;
rng(123)
target_range = linspace(0, 45, pdf_disc+1);
Npat = height(df);
% Iterate over patients
r1_values = zeros(num_samples, 1);
pdf = cell(1, Npat);
for ipat = 1:Npat
    r1_mean=r(ipat);
    for i = 1:num_samples
        temp = normrnd(r1_mean, sigma);
        while temp <= 0
            temp = normrnd(r1_mean, sigma);
        end
        r1_values(i) = temp;
    end
    teta = linspace(1,10,num_samples); %This range referes to the range of ln(teta/C_0)
    tau = teta ./ r1_values; %predicting the time to relapse
    tau = tau/30;
    [temp_pdf, ~] = histcounts(tau, target_range);  % Calculate temp_pdf using histcounts
    bin_centers = 0.5 * (target_range(1:end-1) + target_range(2:end));
    pdf{ipat}= temp_pdf / trapz(bin_centers, temp_pdf);
end

for i = 1:Npat
    % For pdf
    pdf_means(i) = trapz(bin_centers, bin_centers .* pdf{i});
    [~, idx] = max(pdf{i});
    pdf_modes(i) = bin_centers(idx);
    cdf = cumtrapz(bin_centers, pdf{i});
    [~, idx] = min(abs(cdf - 0.5));
    pdf_medians(i) = bin_centers(idx);
end
save('MC_thesis.mat')
% %% Calculate the error for the expected value, mode and median
MSE_mod = 1/Npat* ( sum( (Time_to_progression' - pdf_modes).^2 ) );
fprintf('\n\nMSE_mod = %.3f',MSE_mod);

figure;
scatter(Time_to_progression, pdf_modes, 'filled');
xlabel('True Survival Time (months)');
ylabel('Predicted Survival Time (months)');
title('Validation: True vs. Predicted Survival Time');
axis equal;
grid on;
hold on;
x = linspace(0, 35, 100);
plot(x, x, 'r--');
hold off;
figure;
for i = 1:14
    subplot(5, 4, i);  % Arrange plots in 5 rows and 4 columns
    % Plot each PDF
    plot(bin_centers, pdf{i}, 'b', 'LineWidth', 1); hold on;
    line([Time_to_progression(i),Time_to_progression(i)],[0,0.25],'Color','red','LineWidth',2,'LineStyle','--')
    title(['Patient ' num2str(i)]);
    xlabel('Growth Time (in months)');
    ylabel('Probability Density');
    grid on; % Optional: Add grid for clarity
end