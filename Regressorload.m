% Load the .mat file
data = load('Reg_thesis.mat');

% Get the original data array
pdf_cells = data.pdf_cells;

% Initialize a cell array of size 1x16
regunpdf_STlasso = cell(1, size(pdf_cells, 1));

% Loop over each row of the original array
for i = 1:size(pdf_cells, 1)
    % Assign each row (as a column vector) to the corresponding cell
    regunpdf_STlasso{i} = pdf_cells(i, :)';
end

% Save the new cell array to a .mat file
save('regunpdf_thesis.mat',"regunpdf_STlasso");