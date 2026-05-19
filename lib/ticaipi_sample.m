%% ticaipi_sample
%
% Generate temporally-interleaved CAIPI sampling mask.
%
% Thin wrapper around caipi_sample that derives shift_offset from frame,
% rotating which z-group receives which ky-shift across frames.
%
% N:     [Ny, Nz] grid dimensions
% R:     net integer acceleration factor (Ry * Rz = R)
% frame: 1-based frame index

function omega = ticaipi_sample(N, R, frame)
    Rz = floor(sqrt(R));
    while mod(R, Rz) ~= 0, Rz = Rz - 1; end
    shift_offset = mod(frame - 1, Rz);
    omega = caipi_sample(N, R, shift_offset);
end
