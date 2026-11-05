function omegas = gen_sampling_masks(R)
%% gen_sampling_masks  Generate spatiotemporally incoherent 2D sampling masks.
%   R:      target acceleration factor
%   omegas: Ny x Nz x Nframes logical sampling mask
%
% The sampling method is controlled by samplingMethod in params.m.
% Typical usage:
%   omegas  = gen_sampling_masks(R);
%   ArbEPI(omegas);

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'lib'));
params;

if isempty(randGaussianSigma)
    randGaussianSigma = [Ny, Nz] ./ 6;
end

omegas = false(Ny, Nz, Nframes);
if strcmp(samplingMethod, 'gocart')
    % GOCART is a continuous stream cut into frames, so all frames are generated at once
    omegas = gocart_sample([Ny, Nz], R, Nframes, 'P', gocartP, 'central', gocartCentral, ...
        'W', gocartW, 'seed', gocartSeed);
    return
end
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
        case 'circus'
            % frame index doubles as the time-interleave index (Eq. 4 of the CIRCUS paper)
            omegas(:,:,frame) = circus_sample([Ny, Nz], R, frame, ...
                'variant', circusVariant, 'b', circusB, ...
                'calib', circusCalib, 'crop_corner', circusCropCorner);
        otherwise
            error('Unknown samplingMethod "%s". Choose ''caipi'', ''ticaipi'', ''pd'', ''rand'', ''circus'', or ''gocart''.', samplingMethod);
    end
end
end
