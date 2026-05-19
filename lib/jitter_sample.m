%% jitter_sample
%
% Generate a jittered-grid Cartesian sampling mask. Starts from a regular
% undersampled grid with strides [Ry, Rz] chosen to preserve the Ny:Nz
% aspect ratio, then applies a random integer shift via circshift.
%
% In global-shift mode (jitterGlobalShift = true) the same [dy, dz] offset
% is applied to every cell, so inter-sample spacing is unchanged while the
% grid origin shifts each frame — low blip demand, moderate incoherence.
%
% In per-cell mode (jitterGlobalShift = false) each base sample gets an
% independent random jitter within its Ry x Rz cell, giving higher
% incoherence at the cost of variable inter-sample spacing.
%
% Caller sets rng(frame) before calling to ensure reproducibility.
%
% N:                [Ny, Nz] grid dimensions
% R:                target acceleration factor
% jitterGlobalShift: logical (true = global shift, false = per-cell)

function mask = jitter_sample(N, R, jitterGlobalShift)
    Ny = N(1); Nz = N(2);
    target_samples = floor(Ny * Nz / R);

    % Strides that approximate aspect ratio while keeping Ry*Rz ≈ R
    Ry = max(1, round(sqrt(R * Nz / Ny)));
    Rz = max(1, round(R / Ry));

    % Base regular grid
    base = false(Ny, Nz);
    base(1:Ry:Ny, 1:Rz:Nz) = true;

    if jitterGlobalShift
        % Single random shift applied to the entire grid via circshift
        dy   = randi([0, Ry - 1]) - 1;
        dz   = randi([0, Rz - 1]) - 1;
        mask = circshift(base, dy, 1);
        mask = circshift(mask, dz, 2);
    else
        % Independent per-cell jitter: each base sample moves within its cell
        [base_iy, base_iz] = find(base);
        n_base = numel(base_iy);
        dy_vec = randi([-floor(Ry/2), floor(Ry/2)], n_base, 1);
        dz_vec = randi([-floor(Rz/2), floor(Rz/2)], n_base, 1);
        new_iy = mod(base_iy + dy_vec - 1, Ny) + 1;
        new_iz = mod(base_iz + dz_vec - 1, Nz) + 1;
        lin    = unique(sub2ind([Ny, Nz], new_iy, new_iz));
        mask   = false(Ny, Nz);
        mask(lin) = true;
    end

    % Prune or fill to enforce exact target sample count
    n_sampled = sum(mask, 'all');
    if n_sampled > target_samples
        candidates  = find(mask);
        n_to_remove = n_sampled - target_samples;
        perm        = randperm(numel(candidates));
        mask(candidates(perm(1 : n_to_remove))) = false;
    elseif n_sampled < target_samples
        candidates = find(~mask);
        n_to_add   = target_samples - n_sampled;
        perm       = randperm(numel(candidates));
        mask(candidates(perm(1 : n_to_add))) = true;
    end

    assert(sum(mask, 'all') == target_samples, ...
        'jitter_sample: sample count mismatch (%d vs %d).', sum(mask,'all'), target_samples);
end
