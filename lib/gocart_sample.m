function [omegas, info] = gocart_sample(N, R, Nframes, varargin)
% GOCART_SAMPLE  GOCART (GOlden-angle CArtesian Randomized Time-resolved) ky-kz-t masks.
%
%   omegas = GOCART_SAMPLE([Ny Nz], R, Nframes, ...) returns a logical Ny x Nz x Nframes
%   array in which every frame has EXACTLY floor(Ny*Nz/R) true entries, all distinct.
%   Unlike the other generators this one is a stream: phase encodes (PEs) are produced
%   continuously and cut into frames of floor(Ny*Nz/R) PEs, so it returns every frame at
%   once (frame f depends on frames 1..f-1 through the repeat window W).
%
%   Reference: Zhu Y, Guo Y, Lebel RM, Law M, Nayak KS. GOCART: GOlden-angle CArtesian
%   randomized time-resolved 3D MRI. Magn Reson Imaging 2016;34(7):940-950.
%
%   Algorithm (paper, "Data acquisition phase encode ordering")
%     * Cartesian golden-angle (GA) radial spokes in the ky-kz plane: spoke n runs from the
%       k-space center to the edge at angle n*111.246 deg (the golden angle 180*(sqrt(5)-1)/2),
%       rasterized onto the PE grid.
%     * Each PE a spoke passes through is acquired with probability P (P in (0,1]; P = 1 is
%       plain GA), so part of every spoke is skipped.
%     * P = 1 inside a central region (about 15% of the PEs of a frame), because every spoke
%       passes through the center.
%     * A PE is not repeated within a temporal window W (in acquired PEs).
%     * Corners of the ky-kz plane (outside the inscribed circle of physical k-space) are
%       not acquired.
%     * The stream is binned into frames afterwards. A frame may start or end mid-spoke.
%
%   Adaptations for ArbEPI (not in the paper)
%     * Frame length is fixed at floor(Ny*Nz/R) PEs so that mask2epi can partition every
%       frame into Nshots x ETL; the paper bins retrospectively.
%     * A PE already in the current frame is also rejected, so frames contain no duplicates
%       even when W is shorter than a frame.
%     * Physical k-space is taken as isotropic: ky = (iy-cy)/(Ny/2), kz = (iz-cz)/(Nz/2) in
%       units of kmax, valid when the in-plane and through-plane resolutions are equal. Spokes
%       are therefore isotropic in physical k-space, and the PE spacing differs on the two axes.
%     * No global RNG state is touched (a private RandStream is used), so results do not depend
%       on, and do not disturb, the caller's rng.
%
%   Name-value options
%     'P'        sampling probability outside the central region (default 0.5)
%     'central'  central P = 1 region, as a fraction of the PEs per frame (default 0.15);
%                the region is the circle containing round(central*Nperframe) grid PEs
%     'W'        repeat window, as a fraction of the PEs per frame (default 0.7, Fig. 1 of the paper)
%     'seed'     RandStream seed; also sets the initial spoke angle (default 0)
%
%   [omegas, info] = ... also returns info.rho_c (central radius, units of kmax), info.nspokes,
%   info.Nperframe, info.W (in PEs).

p = inputParser;
addParameter(p, 'P', 0.5, @(x) isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'central', 0.15, @(x) isscalar(x) && x >= 0 && x < 1);
addParameter(p, 'W', 0.7, @(x) isscalar(x) && x >= 0);
addParameter(p, 'seed', 0, @isscalar);
parse(p, varargin{:});
P = p.Results.P;

Ny = N(1); Nz = N(2);
Nper = floor(Ny*Nz/R);
W = round(p.Results.W * Nper);
stream = RandStream('mt19937ar', 'Seed', p.Results.seed);

cy = Ny/2 + 1; cz = Nz/2 + 1;
[IY, IZ] = ndgrid(1:Ny, 1:Nz);
rho = hypot((IY - cy)/(Ny/2), (IZ - cz)/(Nz/2));       % radius in units of kmax
inside = rho <= 1;                                       % corners are not acquired
if nnz(inside) < Nper
    error('gocart_sample:tooFewPEs', 'Only %d PEs inside the k-space circle, need %d per frame.', nnz(inside), Nper);
end
rs = sort(rho(:));
ncen = round(p.Results.central * Nper);
if ncen > 0, rho_c = rs(ncen); else, rho_c = -1; end
central = rho <= rho_c;

golden = 180 * (sqrt(5) - 1) / 2;                        % 111.246 deg
theta0 = 360 * rand(stream, 1) * (p.Results.seed ~= 0);  % seed 0: start at 0 deg
dr = 0.5 / max(Ny, Nz);                                  % fine enough to touch every PE on a spoke
r = 0:dr:1;

omegas = false(Ny, Nz, Nframes);
last = -inf(Ny, Nz);          % stream position (1-based) at which each PE was last acquired
frameOf = zeros(Ny, Nz);      % frame in which each PE was last acquired
pos = 0; f = 1; nin = 0; n = 0;
maxSpokes = 1000 * Nframes;
while f <= Nframes
    n = n + 1;
    if n > maxSpokes
        error('gocart_sample:stalled', 'Could not fill frame %d (parameters leave too few free PEs).', f);
    end
    th = deg2rad(theta0 + (n - 1) * golden);
    iy = round(cy + r*cos(th)*Ny/2);
    iz = round(cz + r*sin(th)*Nz/2);
    ok = iy >= 1 & iy <= Ny & iz >= 1 & iz <= Nz;
    lin = unique(sub2ind([Ny Nz], iy(ok), iz(ok)), 'stable');   % PEs along the spoke, center -> edge
    lin = lin(inside(lin));
    for q = lin
        if pos - last(q) <= W && last(q) > 0, continue; end        % inside the repeat window
        if frameOf(q) == f, continue; end                          % already in this frame
        if ~central(q) && rand(stream) >= P, continue; end         % skipped with probability 1-P
        pos = pos + 1; last(q) = pos; frameOf(q) = f;
        omegas(IY(q), IZ(q), f) = true; nin = nin + 1;
        if nin == Nper
            f = f + 1; nin = 0;
            if f > Nframes, break; end
        end
    end
end

info = struct('rho_c', rho_c, 'nspokes', n, 'Nperframe', Nper, 'W', W, 'P', P);
end
