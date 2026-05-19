%% caipi_sample
%
% Generate regular CAIPI shifted 2D sampling mask
%
% Last modified May 14th, 2026. Rex Fung

function omega = caipi_sample(N, R, shift_offset)
    if nargin < 3, shift_offset = 0; end

    Ny = N(1); Nz = N(2);

    assert(Ny >= 1 && Nz >= 1, 'Dimensions must be >= 1');
    assert(R >= 1 && round(R) == R, 'R must be a positive integer');
    assert(round(shift_offset) == shift_offset && shift_offset >= 0, ...
           'shift_offset must be a non-negative integer');

    % Most balanced factorization Ry * Rz = R, with Rz <= sqrt(R)
    Rz = floor(sqrt(R));
    while mod(R, Rz) ~= 0
        Rz = Rz - 1;
    end
    Ry = R / Rz;
    caipi_z = Rz;

    omega = zeros(Ny, Nz);
    omega(1:Ry:Ny, 1:Rz:Nz) = 1;

    for z = 1:caipi_z
        shift_amount = mod(z - 1 + shift_offset, caipi_z);
        cols = (1 + Rz*(z-1)) : caipi_z*Rz : Nz;
        if isempty(cols), continue; end
        omega(:, cols) = circshift(omega(:, cols), shift_amount, 1);
    end
end
