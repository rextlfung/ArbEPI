function omegas = gen_sampling_masks(R)
%% gen_sampling_masks  Generate spatiotemporally incoherent 2D sampling masks.
%   R:      target acceleration factor
%   omegas: Ny x Nz x Nframes logical sampling mask
%
% The sampling method is controlled by samplingMethod in params.m.
% Typical usage:
%   omegas  = gen_sampling_masks(R);
%   ArbEPI(omegas);

projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;

omegas = false(Ny, Nz, Nframes);
for frame = 1:Nframes
    rng(frame); % deterministic per-frame seed for reproducibility
    switch samplingMethod
        case 'caipi'
            omegas(:,:,frame) = caipi_sample([Ny, Nz], R);
        case 'ticaipi'
            Rz = floor(sqrt(R));
            while mod(R, Rz) ~= 0, Rz = Rz - 1; end
            shift_offset = mod(frame - 1, Rz);
            omegas(:,:,frame) = ticaipi_sample([Ny, Nz], R, shift_offset);
        case 'pd'
            omegas(:,:,frame) = pd_sample([Ny, Nz], R);
        case 'rand'
            weights = gen_gaussian_pdf([Ny, Nz], [Ny, Nz]./10);
            omegas(:,:,frame) = rand_sample([Ny, Nz], R, weights);
        case 'radial'
            omegas(:,:,frame) = radial_sample([Ny, Nz], R, frame, radialCalib);
        case 'jitter'
            omegas(:,:,frame) = jitter_sample([Ny, Nz], R, jitterGlobalShift);
        otherwise
            error('Unknown samplingMethod "%s". Choose ''caipi'', ''ticaipi'', ''pd'', ''rand'', ''radial'', or ''jitter''.', samplingMethod);
    end
end
end
