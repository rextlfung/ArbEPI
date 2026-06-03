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

if isempty(randGaussianSigma)
    randGaussianSigma = [Ny, Nz] ./ 6;
end

omegas = false(Ny, Nz, Nframes);
for frame = 1:Nframes
    rng(frame); % deterministic per-frame seed for reproducibility
    switch samplingMethod
        case 'caipi'
            omegas(:,:,frame) = caipi_sample([Ny, Nz], R);
        case 'ticaipi'
            omegas(:,:,frame) = ticaipi_sample([Ny, Nz], R, frame);
        case 'pd'
            omegas(:,:,frame) = pd_sample([Ny, Nz], R, ...
                'calib', pdCalib, 'crop_corner', pdCropCorner, 'decay', pdDecay);
        case 'rand'
            weights = gen_gaussian_pdf([Ny, Nz], randGaussianSigma);
            omegas(:,:,frame) = rand_sample([Ny, Nz], R, weights);
        otherwise
            error('Unknown samplingMethod "%s". Choose ''caipi'', ''ticaipi'', ''pd'', or ''rand''.', samplingMethod);
    end
end
end
