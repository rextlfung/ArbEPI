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
            omegas(:,:,frame) = caipi_sample([Ny, Nz], caipiR, caipiShift);
        case 'pd'
            omegas(:,:,frame) = pd_sample([Ny, Nz], R);
        case 'rand'
            weights = gen_gaussian_pdf([Ny, Nz], [Ny, Nz]./10);
            omegas(:,:,frame) = rand_sample([Ny, Nz], R, weights);
        otherwise
            error('Unknown samplingMethod "%s". Choose ''caipi'', ''pd'', or ''rand''.', samplingMethod);
    end
end
end
