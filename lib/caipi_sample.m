%% caipi_sample
%
% Generate regular CAIPI shifted 2D sampling mask
%
% Last modified Oct 23rd, 2025. Rex Fung

function omega = caipi_sample(N, R, caipi_z)
    % Unpack input arguments
    Ny = N(1); Nz = N(2);
    Ry = R(1); Rz = R(2);
    
    % Validate input arguments
    assert(Ny >= 1 && Nz >= 1, ...
           'Dimensions must be >= 1');
    assert(Ry >= 1 && Rz >= 1 && round(Ry) == Ry && round(Rz) == Rz, ...
           'Acceleration factors must be an integer >= 1');
        assert(round(caipi_z) == caipi_z && caipi_z >= 1, ...
           'caipi_z must be an integer >= 1');

    % Create sampling mask
    omega = zeros(Ny, Nz);
    omega(1:Ry:Ny, 1:Rz:Nz) = 1;

    % Apply CAIPI shifts
    for z = 1:caipi_z
        omega(:, (1+Rz*(z-1)):caipi_z*Rz:Nz) = circshift(omega(:, (1+Rz*(z-1)):caipi_z*Rz:Nz), z-1, 1); 
    end
end 