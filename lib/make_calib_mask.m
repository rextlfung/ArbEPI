%% make_calib_mask  Build a centered rectangular calibration region mask.
%
% Returns a Ny×Nz logical mask with a (calib × calib) block at the center
% set to true.  If calib == 0, returns an all-false mask.
%
% N:     [Ny, Nz]
% calib: half-size of the calibration region in samples (scalar, 0 = disabled)

function mask = make_calib_mask(N, calib)
    Ny = N(1); Nz = N(2);
    mask = false(Ny, Nz);
    if calib > 0
        y0 = max(1, floor(Ny/2 - calib/2) + 1);
        y1 = min(Ny, floor(Ny/2 + calib/2));
        z0 = max(1, floor(Nz/2 - calib/2) + 1);
        z1 = min(Nz, floor(Nz/2 + calib/2));
        mask(y0:y1, z0:z1) = true;
    end
end
