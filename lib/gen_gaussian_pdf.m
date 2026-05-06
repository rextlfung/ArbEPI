function pdf = gen_gaussian_pdf(dims, sigma)
% GEN_GAUSSIAN_PDF Generates a 2D Gaussian density map
%
% Inputs:
%   dims:   [Nx, Ny] vector (e.g., [256, 256])
%   sigma:  Standard deviation. Can be a scalar (isotropic) 
%           or a [sig_x, sig_y] vector (anisotropic).
%           (Units are in pixels)
%
% Output:
%   pdf:    [Nx, Ny] matrix with values between 0 and 1

    Nx = dims(1);
    Ny = dims(2);
    
    % Handle scalar vs vector sigma
    if isscalar(sigma)
        sig_x = sigma;
        sig_y = sigma;
    else
        sig_x = sigma(1);
        sig_y = sigma(2);
    end

    % 1. Create centered grid coordinates
    % Range from -Nx/2 to Nx/2
    cx = ceil(Nx/2);
    cy = ceil(Ny/2);
    [Y, X] = meshgrid((1:Ny) - cy, (1:Nx) - cx);

    % 2. Calculate Gaussian
    % Formula: exp( - (x^2/(2*sx^2) + y^2/(2*sy^2)) )
    pdf = exp(-((X.^2) ./ (2 * sig_x^2) + (Y.^2) ./ (2 * sig_y^2)));

    % 3. (Optional) Ensure peak is exactly 1.0 (it usually is by math, but safe to enforce)
    pdf = pdf ./ max(pdf(:));

end