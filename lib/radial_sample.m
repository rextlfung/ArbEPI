%% radial_sample
%
% Generate a Cartesian k-space mask approximating a stack-of-stars radial
% trajectory. Each frame uses Nspokes spokes whose base angles are rotated
% by the golden angle between frames, providing temporal incoherence.
%
% N:          [Ny, Nz] grid dimensions
% R:          target acceleration factor (determines Nspokes)
% frame:      1-based frame index (drives the golden-angle rotation)
% radialCalib: calibration region half-size in samples (0 = disabled).
%              The center block is always fully sampled, and the DC point
%              is always included regardless of this setting.

function mask = radial_sample(N, R, frame, radialCalib)
    Ny = N(1); Nz = N(2);
    target_samples = floor(Ny * Nz / R);

    phi = (1 + sqrt(5)) / 2;

    % Number of spokes per frame (Nyquist radial criterion)
    Nspokes = max(1, round(pi/2 * max(Ny, Nz) / R));

    % Calibration region (always-sampled center block + DC point)
    cy = ceil(Ny / 2);
    cz = ceil(Nz / 2);

    calib_mask = make_calib_mask([Ny, Nz], radialCalib);
    calib_mask(cy, cz) = true;  % DC point always included; atan2(0,0) is degenerate

    % Centered coordinate grids
    [IY, IZ] = ndgrid(1:Ny, 1:Nz);
    Yc = IY - cy;
    Zc = IZ - cz;

    % Spoke angles for this frame: evenly spaced in [0, pi), rotated by
    % (frame-1) golden angles
    theta_base   = (0 : Nspokes - 1) * pi / Nspokes;
    theta_offset = mod((frame - 1) * pi / phi, pi);
    thetas       = mod(theta_base + theta_offset, pi);

    half_dtheta = pi / (2 * Nspokes);

    % Polar angle of every grid point in [0, pi)
    angles = mod(atan2(Yc, Zc), pi);

    % Mark grid points within half_dtheta of any spoke (vectorized over spokes)
    diff_angles = abs(angles - reshape(thetas, 1, 1, []));
    diff_angles = min(diff_angles, pi - diff_angles);
    radial_mask = any(diff_angles <= half_dtheta, 3);
    radial_mask(cy, cz) = true;

    combined    = radial_mask | calib_mask;
    all_sampled = find(combined);
    n_sampled   = numel(all_sampled);

    if n_sampled > target_samples
        % Prune outermost (largest-radius) non-protected samples first
        removable  = all_sampled(~calib_mask(all_sampled));
        radius_sq  = Yc(removable).^2 + Zc(removable).^2;
        [~, ord]   = sort(radius_sq, 'descend');
        n_to_remove = min(n_sampled - target_samples, numel(removable));
        mask = combined;
        mask(removable(ord(1 : n_to_remove))) = false;
    elseif n_sampled < target_samples
        % Fill randomly from unsampled points (consumes rng(frame) state)
        candidates = find(~combined);
        n_to_add   = target_samples - n_sampled;
        perm       = randperm(numel(candidates));
        mask = combined;
        mask(candidates(perm(1 : n_to_add))) = true;
    else
        mask = combined;
    end

    assert(sum(mask, 'all') == target_samples, ...
        'radial_sample: sample count mismatch (%d vs %d).', sum(mask,'all'), target_samples);
end
