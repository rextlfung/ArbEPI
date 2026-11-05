function export_ball_masks(outfile, Nframes)
% EXPORT_BALL_MASKS  Export the three SI-demo ky-kz-t masks as a v6 .mat for ArbEPI-python.
%   Ny x Nz = 90 x 60, R = 6 (900 PEs per frame = 12 shots x ETL 75), Nframes frames each.
%   omegas_radialpd : Poisson disc exactly as gen_sampling_masks builds it (rng(frame) per frame,
%                     params.m's pdCalib / pdCropCorner / pdDecay). Its pd_sample.m is identical to
%                     the manuscript's ArbEPI SHA 02984ac. The "radial" in Radial-PD is the
%                     echo-train ordering (ArbEPI-python mask2epi_radial), not the mask.
%   omegas_circus   : CIRCUS, params.m's circus* options.
%   omegas_gocart   : GOCART, params.m's gocart* options.
if nargin < 2, nf_req = 75; else, nf_req = Nframes; end
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(repo); addpath(fullfile(repo, 'lib')); params;   % params.m defines Ny, Nz, R, Nframes, ...: set ours after
Ny = 90; Nz = 60; R = 6; Nper = Ny*Nz/R; Nframes = nf_req;
omegas_radialpd = false(Ny, Nz, Nframes); omegas_circus = omegas_radialpd;
for f = 1:Nframes
    rng(f);
    omegas_radialpd(:,:,f) = pd_sample([Ny Nz], R, 'calib', pdCalib, 'crop_corner', pdCropCorner, 'decay', pdDecay);
    omegas_circus(:,:,f) = circus_sample([Ny Nz], R, f, 'variant', circusVariant, 'b', circusB, ...
        'calib', circusCalib, 'crop_corner', circusCropCorner);
end
omegas_gocart = gocart_sample([Ny Nz], R, Nframes, 'P', gocartP, 'central', gocartCentral, 'W', gocartW, 'seed', gocartSeed);
for nm = {'radialpd', 'circus', 'gocart'}
    om = eval(['omegas_' nm{1}]);
    c = squeeze(sum(sum(om, 1), 2));
    assert(all(c == Nper), '%s: frame counts not all %d', nm{1}, Nper);
    d = 0; for f = 2:Nframes, d = d + isequal(om(:,:,f), om(:,:,f-1)); end
    fprintf('%-9s %d frames x %d PEs, identical consecutive frames: %d, distinct frames: %d\n', nm{1}, Nframes, Nper, d, size(unique(reshape(om, [], Nframes)', 'rows'), 1));
end
git_sha = strtrim(evalc('system(sprintf(''git -C %s rev-parse --short HEAD'', repo));'));
meta = struct('Ny', Ny, 'Nz', Nz, 'R', R, 'Nframes', Nframes, 'pdCalib', pdCalib, 'pdDecay', pdDecay, ...
    'pdCropCorner', pdCropCorner, 'circusVariant', circusVariant, 'circusB', circusB, 'circusCalib', circusCalib, ...
    'gocartP', gocartP, 'gocartCentral', gocartCentral, 'gocartW', gocartW, 'gocartSeed', gocartSeed, 'ArbEPI_revised_base_commit', git_sha); %#ok<NASGU>
save(outfile, 'omegas_radialpd', 'omegas_circus', 'omegas_gocart', 'meta', '-v6');
end
