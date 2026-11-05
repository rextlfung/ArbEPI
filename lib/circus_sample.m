function [mask, info] = circus_sample(N, R, frame, varargin)
% CIRCUS_SAMPLE  CIRCUS (CIRcular Cartesian UnderSampling) ky-kz mask, one frame.
%
%   mask = CIRCUS_SAMPLE([Ny Nz], R, frame, ...) returns a logical Ny x Nz mask with
%   EXACTLY floor(Ny*Nz/R) true entries, for the 1-based time-frame index 'frame'
%   (frame = 1 is t = 0 in the notation of the paper).  Fully deterministic given
%   (N, R, frame, options): no global RNG is used.
%
%   Reference: Liu J, Saloner D. Accelerated MRI with CIRcular Cartesian
%   UnderSampling (CIRCUS). Quant Imaging Med Surg 2014;4(1):57-67.
%
%   Algorithm
%     * The grid is decomposed into nested concentric "rings" (squares of side
%       J = 2,4,...,Ny in the paper's symmetric case, K = 4J-4 points each; the
%       K's sum to Ny^2).  Ring j holds M_j of its K_j points.
%     * Points on a ring are indexed clockwise from the bottom-left corner and the
%       selected ones are i = floor(mod((m + b*J + t*M_j)*r - 1, 1) * K), with
%       m = 0..M_j-1, r = (1+sqrt(5))/2, shift b*J (Eq. 2) and interleave offset
%       t*M_j (Eq. 4).  Golden-ratio collisions (equal floor()) are skipped by
%       advancing m until M_j distinct points are found.
%     * M_j is constant (= variable density) except that rings with K_j <= M are
%       fully sampled and the remaining budget is re-distributed ("water-filling"):
%       M_j = min(K_j, mu) with mu chosen so the total hits the target; the
%       fractional remainder is given as +1 to evenly spaced unsaturated rings.
%
%   [mask, info] = ... also returns diagnostics: info.M (points per ring), info.Tsym
%   (symmetric-grid budget, rowselect only), info.keep (kept kz columns, rowselect
%   only) and info.nfix (signed number of points trimmed(-)/added(+) by the exact
%   count fix-up, rowselect only; 0 for 'rect').
%
%   Name-value options
%     'variant'     'rect' (default) | 'rowselect'
%        'rect'      (NOT from the paper) nested rectangles in index space with
%                    aspect ratio Nz/Ny, i.e. squares in physical k-space when
%                    the voxels are isotropic.  Each grid point is assigned to
%                    ring ceil(max(dy, dz*Ny/Nz)) where dy, dz are the distances
%                    (in pixels, +0.5) from the k-space centre; this partitions
%                    the grid exactly (no gaps/collisions).
%        'rowselect' (the paper's recipe for Nz < Ny) build the pattern on a
%                    symmetric Ny x Ny matrix, keep Nz of the Ny kz-columns
%                    selected with the golden-ratio profile, and COMPACT the kept
%                    columns into an Ny x Nz grid.  The kept column set is the
%                    same for every frame.
%     'b'           integer shift multiplier (Eq. 2). 0 = pure radial lines. Default 30.
%     'calib'       [cy cz] fully-sampled central block (0 = none, default [0 0]);
%                   placed like pd_sample's calibration region.
%     'rect_alloc'  'smooth' (default) | 'const'. Only for variant 'rect': 'const' is the
%                   paper's literal constant M per ring, which on the 3:2 index grid
%                   leaves rings of alternating thickness (K = 28,36,16,48,56,24,...)
%                   and hence a period-3 density modulation along ky (coherent PSF
%                   sidelobe at dy = Ny/3).  'smooth' scales M_j by K_j/movmean(K,3),
%                   i.e. keeps the constant-M density law but evaluated on the
%                   locally averaged ring size.
%     'crop_corner' logical, skip points outside the inscribed ellipse
%                   (default false).
%
%   Example:
%       mask = circus_sample([90 60], 6, 1);            % rect, b = 30
%       mask = circus_sample([90 60], 6, 3, 'variant', 'rowselect', 'b', 10);

    if nargin < 3 || isempty(frame), frame = 1; end
    p = inputParser;
    p.addParameter('variant', 'rect', @(x) any(strcmpi(x, {'rect', 'rowselect'})));
    p.addParameter('b', 30, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    p.addParameter('calib', [0 0], @(x) isnumeric(x) && numel(x) == 2);
    p.addParameter('crop_corner', false, @islogical);
    p.addParameter('rect_alloc', 'smooth', @(x) any(strcmpi(x, {'smooth', 'const'})));
    p.parse(varargin{:});
    variant = lower(p.Results.variant);
    b       = p.Results.b;
    calib   = p.Results.calib(:).';
    crop    = p.Results.crop_corner;
    rectAlloc = p.Results.rect_alloc;

    Ny = N(1); Nz = N(2);
    assert(R >= 1, 'R must be >= 1');
    assert(frame >= 1 && round(frame) == frame, 'frame must be a positive integer (1-based)');
    target = floor(Ny*Nz/R);
    t = frame - 1;

    % fully sampled calibration block on the final grid (same indexing as pd_sample)
    calibMask = false(Ny, Nz);
    if any(calib > 0)
        ys = max(1, floor(Ny/2 - calib(1)/2) + 1); ye = min(Ny, floor(Ny/2 + calib(1)/2));
        zs = max(1, floor(Nz/2 - calib(2)/2) + 1); ze = min(Nz, floor(Nz/2 + calib(2)/2));
        calibMask(ys:ye, zs:ze) = true;
    end
    nCalib = nnz(calibMask);
    assert(nCalib <= target, 'calibration block larger than the sample budget');
    budget = target - nCalib;

    switch variant
        case 'rect'
            rings = build_rings(Ny, Nz, calibMask, crop);
            [mask, info] = sample_rings(rings, budget, t, b, Ny, Nz, [], 0, strcmpi(rectAlloc, 'smooth'));
        case 'rowselect'
            assert(Nz <= Ny, 'rowselect requires Nz <= Ny');
            keep = golden_rows(Ny, Nz);                 % kept kz columns of the symmetric matrix
            blockedSym = false(Ny, Ny);                 % calib block expressed in symmetric coords
            blockedSym(:, keep) = calibMask;
            rings = build_rings(Ny, Ny, blockedSym, crop);
            % The symmetric pattern must contain more than 'budget' points because
            % ~(Ny-Nz)/Ny of them fall on discarded columns.  Calibrate the symmetric
            % budget (frame-independent, fixed 4-frame average) so the mean kept
            % count matches the target; remaining per-frame deviation is fixed exactly.
            Tsym = round(budget * Ny / Nz);
            for it = 1:6
                cnt = 0;
                for tt = 0:3
                    cnt = cnt + count_kept(rings, Tsym, tt, b, keep, Ny, false);
                end
                cnt = cnt / 4;
                if cnt == 0, break; end
                Tnew = round(Tsym * budget / cnt);
                if Tnew == Tsym, break; end
                Tsym = Tnew;
            end
            [mask, info] = sample_rings(rings, Tsym, t, b, Ny, Nz, keep, budget, false);
            info.Tsym = Tsym; info.keep = keep;
    end
    mask = mask | calibMask;
    assert(nnz(mask) == target, 'circus_sample: failed to hit exact sample count (%d vs %d)', nnz(mask), target);
end

% ------------------------------------------------------------------------------
function rings = build_rings(Ny, Nz, blocked, crop)
% Partition the Ny x Nz grid into nested rings; return cell array (outer = last) of
% [y z] lists ordered clockwise from the bottom-left corner.
    cy = (Ny + 1)/2; cz = (Nz + 1)/2;          % k-space centre sits between pixels for even N
    [Z, Y] = meshgrid(1:Nz, 1:Ny);
    dyS = Y - cy; dzS = Z - cz;
    level = ceil(max(abs(dyS) + 0.5, (abs(dzS) + 0.5) * (Ny/Nz)) - 1e-9);   % 1..ceil(Ny/2)
    valid = ~blocked;
    if crop
        valid = valid & ((dyS/(Ny/2)).^2 + (dzS/(Nz/2)).^2 <= 1);
    end
    % angle measured on physically normalised coordinates, starting at the corner (-,-)
    ang = mod(atan2(dzS/Nz, dyS/Ny) + 3*pi/4, 2*pi);
    nlev = max(level(:));
    rings = cell(1, nlev);
    for j = 1:nlev
        idx = find(level == j & valid);
        if isempty(idx), rings{j} = zeros(0, 2); continue; end
        rr = hypot(dyS(idx)/Ny, dzS(idx)/Nz);
        [~, o] = sortrows([ang(idx), rr]);
        idx = idx(o);
        rings{j} = [Y(idx), Z(idx)];
    end
    rings = rings(~cellfun(@isempty, rings));
    if Ny == Nz && ~any(blocked(:)) && ~crop && mod(Ny, 2) == 0
        K = cellfun(@(c) size(c, 1), rings);
        J = 2*(1:numel(rings));
        assert(isequal(K, 4*J - 4) && sum(K) == Ny^2, 'ring decomposition does not match K = 4J-4');
    end
end

function M = allocate(K, T, w)
% Water-filling allocation of T points to rings: M_j = min(K_j, mu*w_j), integers with
% sum(M) = min(T, sum(K)).  w = ones gives the paper's constant M (rings with
% K <= M fully sampled); w = K./movmean(K,3) is used for 'rect' (see header).
% Fractional parts are rounded by error diffusion so extras are spread evenly.
    K = K(:).'; w = w(:).'; T = min(T, sum(K));
    lo = 0; hi = max(K ./ w);
    for it = 1:80
        mu = (lo + hi)/2;
        if sum(min(K, mu*w)) < T, lo = mu; else, hi = mu; end
    end
    Mreal = min(K, hi*w);
    q = floor(cumsum(Mreal) + 0.5);
    q(end) = T;
    M = diff([0, q]);
    assert(sum(M) == T && all(M >= 0) && all(M <= K), 'allocation failure');
end

function w = ring_weights(K, smooth)
% smooth = false: w = 1 (paper).  smooth = true: w = K./(3-ring moving mean of K), which
% removes the period-3 density modulation caused by rect rings whose thickness
% (1 ky row, 2/3 kz column) does not match the integer grid.
    K = K(:).';
    if smooth
        Kbar = movmean(K, 3, 'Endpoints', 'shrink');
        w = K ./ Kbar;
    else
        w = ones(size(K));
    end
end

function seq = golden_seq(K, J, Mj, t, b)
% Ordered list of distinct ring-point indices (1..K) generated by Eqs. 1,2,4:
% i = floor(mod((q + b*J + t*Mj)*r - 1, 1)*K), q = 0,1,2,...; duplicates skipped.
    r = (1 + sqrt(5))/2;
    Q = 12*K + 50;
    q = (0:Q-1).' + b*J + t*Mj;
    i = floor(mod(q*r - 1, 1)*K) + 1;
    i = min(i, K);
    seq = unique(i, 'stable');
    if numel(seq) < K                           % append never-hit points (very rare) in index order
        seq = [seq; setdiff((1:K).', seq)];
    end
end

function [sel, ext] = ring_pick(ring, Mj, J, t, b, keepCol)
% sel: selected ring-point indices; ext: following candidates (for top-up).
    K = size(ring, 1);
    if Mj >= K
        sel = (1:K).'; ext = zeros(0, 1);
    else
        seq = golden_seq(K, J, Mj, t, b);
        sel = seq(1:Mj); ext = seq(Mj+1:end);
    end
    sel = sel(keepCol(ring(sel, 2)));
    ext = ext(keepCol(ring(ext, 2)));
end

function n = count_kept(rings, T, t, b, keep, Ny, smooth)
    K = cellfun(@(c) size(c, 1), rings);
    M = allocate(K, T, ring_weights(K, smooth));
    keepCol = false(1, Ny); keepCol(keep) = true;
    n = 0;
    for j = 1:numel(rings)
        sel = ring_pick(rings{j}, M(j), 2*j, t, b, keepCol);
        n = n + numel(sel);
    end
end

function [mask, info] = sample_rings(rings, T, t, b, Ny, Nz, keep, target, smooth)
% Select points on every ring, map to the final Ny x Nz grid and (if target > 0,
% i.e. in the rowselect case where the discarded-column losses are unknown in
% advance) trim/top-up to exactly 'target' from the outermost rings inward.
    K = cellfun(@(c) size(c, 1), rings);
    M = allocate(K, T, ring_weights(K, smooth));
    nr = numel(rings);
    isSym = ~isempty(keep);                       % rowselect: rings live on the symmetric grid
    if isSym
        keepCol = false(1, Ny); keepCol(keep) = true;
    else
        keepCol = true(1, Nz);
    end
    sels = cell(1, nr); exts = cell(1, nr);
    for j = 1:nr
        [sels{j}, exts{j}] = ring_pick(rings{j}, M(j), 2*j, t, b, keepCol);
    end
    info.M = M; info.nfix = 0;
    if target > 0
        n = sum(cellfun(@numel, sels)); n0 = n;
        j = nr; guard = 0;
        while n > target && guard < 1e5               % trim: drop last-picked point, outer rings first
            if ~isempty(sels{j}), sels{j}(end) = []; n = n - 1; end
            j = j - 1; if j < 1, j = nr; end; guard = guard + 1;
        end
        while n < target && guard < 1e5               % top up with the next golden candidates
            if ~isempty(exts{j})
                sels{j}(end+1, 1) = exts{j}(1); exts{j}(1) = []; n = n + 1;
            end
            j = j - 1; if j < 1, j = nr; end; guard = guard + 1;
        end
    end
    if target > 0, info.nfix = n - n0; info.n_before_fix = n0; end
    mask = false(Ny, Nz);
    newcol = zeros(1, numel(keepCol)); newcol(keep) = 1:numel(keep);   % compaction map (rowselect only)
    for j = 1:nr
        yz = rings{j}(sels{j}, :);
        zc = yz(:, 2);
        if isSym, zc = newcol(zc).'; end
        mask(sub2ind([Ny, Nz], yz(:, 1), zc)) = true;
    end
end

function keep = golden_rows(Ny, Nz)
% Choose Nz of Ny kz columns with the golden-ratio profile i = floor(mod(m*r-1,1)*Ny),
% skipping repeats, starting at the smallest sequence offset m0 that leaves exactly
% floor(Nz/2) kept columns in the lower half (so the compacted grid stays centred
% on the k-space centre).  Returned sorted ascending.
    r = (1 + sqrt(5))/2;
    for m0 = 0:Ny-1
        i = floor(mod(((0:12*Ny).' + m0)*r - 1, 1)*Ny) + 1;
        u = unique(i, 'stable');
        keep = sort(u(1:Nz)).';
        if nnz(keep <= Ny/2) == floor(Nz/2), return; end
    end
    error('golden_rows: no centred offset found');
end
