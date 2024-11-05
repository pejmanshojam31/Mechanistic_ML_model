close all
clear
% Define the directory and filename
directory = 'F:\Olivier';
filename = fullfile(directory, 'dataset_thesis_complete.csv');
df = readtable(filename);
df = df(df.survival_time < 1280, :); % Added growth condition
Npat = height(df);
Time_to_progression = df.LAST_MR/30; % True survival days
Age = df.Age;
Total_resection = df.Total_resection;
PC3 = df.PC3;
Vol_rat = df.Vol_rat;
mu2mu1 = df.mu2mu1;
Vol_nec= df.Vol_Nec;
sex = df.sex;
pdf_disc = 100; % Discretization of the modelable pdf
LAST_MRI_RANGE = linspace(0, 45, pdf_disc);
data = [Time_to_progression, Age, Total_resection, Vol_nec, PC3];

% Initialize variables for LOOCV
means = zeros(1, Npat);
medians = zeros(1, Npat);
modes = zeros(1, Npat);
unpdfLOOCV = cell(1, Npat); % Renamed to unpdfLOOCV

% Loop over the patients for LOOCV
for i = 1:Npat
    % Split the data into training and test sets
    train_data = data([1:i-1, i+1:Npat], :);
    test_data = data(i, :);
    
    % Calculate bandwidths dynamically for the training data
    bandwidths = zeros(1, size(train_data, 2));
    for j = 1:size(train_data, 2)
        bandwidths(j) = std(train_data(:, j)) * (4 / ((Npat - 1) * (size(train_data, 2) + 2))) ^ (1 / (size(train_data, 2) + 4));
    end

    % Ensure set_points is properly constructed
    set_points = [LAST_MRI_RANGE', repmat(test_data(2:end), pdf_disc, 1)];
    
    % Estimate the density using training data
    temp = mvksdensity(train_data, set_points, 'Bandwidth', bandwidths);
    
    % Normalize the distribution
    un_pdf = temp / trapz(LAST_MRI_RANGE, temp);
    unpdfLOOCV{i} = un_pdf; % Assign to unpdfLOOCV
    
    % Mean
    means(i) = trapz(LAST_MRI_RANGE,LAST_MRI_RANGE .* un_pdf');
    % Median
    cdf = cumtrapz(LAST_MRI_RANGE, un_pdf);
    [~, idx] = min(abs(cdf - 0.5));
    medians(i) = LAST_MRI_RANGE(idx);
    % Mode (Peak)
    [~, idx] = max(un_pdf);
    modes(i) = LAST_MRI_RANGE(idx);
end
MSE_mod = 1/Npat * ( sum( (Time_to_progression' - modes).^2 ) );
fprintf('\n\nMSE_mod = %.3f',MSE_mod);
true_Survival_days = Time_to_progression; % Since you already have this
predicted_Survival_days = modes; % Using modes as the predicted value
% Save the unpdfLOOCV data
save('unpdf_KDE_thesis.mat', 'unpdfLOOCV');
