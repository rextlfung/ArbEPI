%% ticaipi_sample
%
% Generate temporally-interleaved CAIPI sampling mask.
%
% Thin wrapper around caipi_sample that passes a per-frame shift_offset,
% rotating which z-group receives which ky-shift across frames.
% Caller computes: shift_offset = mod(frame - 1, Rz), where Rz is the
% largest factor of R that is <= sqrt(R).
%
% N:            [Ny, Nz] grid dimensions
% R:            net integer acceleration factor (Ry * Rz = R)
% shift_offset: non-negative integer offset applied to all z-group shifts

function omega = ticaipi_sample(N, R, shift_offset)
    omega = caipi_sample(N, R, shift_offset);
end
