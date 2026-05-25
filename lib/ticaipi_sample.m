%% ticaipi_sample
%
% Generate temporally-interleaved CAIPI sampling mask.
%
% Cycles through all R = Ry * Rz unique (y_shift, kz_offset) combinations
% so that every frame gets a distinct mask and the union over R consecutive
% frames covers k-space exactly once.
%
% The CAIPI pattern is fixed (shift_offset = 0).  Full coverage is achieved
% by two independent global circshifts applied to that base pattern:
%   kz_offset (0..Rz-1) shifts the sampled kz columns (dim 2)
%   y_shift   (0..Ry-1) shifts the sampled ky rows    (dim 1)
%
% N:     [Ny, Nz] grid dimensions
% R:     net integer acceleration factor (Ry * Rz = R)
% frame: 1-based frame index

function omega = ticaipi_sample(N, R, frame)
    Rsmall = floor(sqrt(R));
    while mod(R, Rsmall) ~= 0, Rsmall = Rsmall - 1; end
    if N(1) >= N(2)
        Rz = Rsmall;
    else
        Rz = R / Rsmall;
    end
    frame_idx = mod(frame - 1, R);
    kz_offset = mod(frame_idx, Rz);
    y_shift   = floor(frame_idx / Rz);
    omega = caipi_sample(N, R, 0);
    omega = circshift(omega, y_shift,   1);
    omega = circshift(omega, kz_offset, 2);
end
