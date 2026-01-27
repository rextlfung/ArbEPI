function mask = rand_sample(dims, R, pdf_in)
% RAND_SAMPLE Randomly samples a 2D grid based on weights
%
% Inputs:
%   dims:   [Nx, Ny] vector specifying grid size
%   R:      Acceleration factor (Scalar)
%   pdf_in: Matrix of sampling probabilities (weights).
%           Must be size [Nx, Ny]. If empty, uses uniform weights.
%
% Output:
%   mask:   Logical matrix of size [Nx, Ny]

    % 1. Unpack dimensions
    Nx = dims(1);
    Ny = dims(2);
    
    % 2. Calculate required sample count
    total_points = Nx * Ny;
    num_samples = round(total_points / R);

    % 3. Validate and Assign Weights
    if isempty(pdf_in)
        weights = ones(Nx, Ny);
    else
        weights = pdf_in;
        if ~isequal(size(weights), [Nx, Ny])
            error('PDF dimensions must match the input [Nx, Ny] vector.');
        end
    end

    % 4. Weighted Random Sampling (Key-Sorting Method)
    % Generate random keys u^(1/w). Heavier weights = larger keys.
    u = rand(Nx, Ny);
    
    % Small epsilon added to weights to avoid division by zero if PDF has 0s
    safe_weights = weights + eps; 
    keys = u .^ (1 ./ safe_weights);

    % 5. Sort and Select
    [~, sort_idx] = sort(keys(:), 'descend');
    
    mask = false(Nx, Ny);
    chosen_indices = sort_idx(1:num_samples);
    mask(chosen_indices) = true;

end