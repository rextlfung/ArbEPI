function [rf, gzSS, gzSSR] = make_excitation_pulse(alpha, rfDur, rfTB, fov, sys, CRT)
%% make_excitation_pulse  Create slab-selective sinc excitation pulse.
%   alpha: flip angle (degrees)
%   rfDur: RF pulse duration (s)
%   rfTB:  RF time-bandwidth product
%   fov:   field of view [x y z] (m); slab thickness = 0.9*fov(3)
%   sys:   Pulseq system struct (mr.opts)
%   CRT:   common raster time for GE compatibility (s)

% Target a slightly thinner slab to alleviate aliasing
[rf, gzSS, gzSSR, ~] = mr.makeSincPulse(alpha/180*pi, ...
    'duration', rfDur, ...
    'sliceThickness', 0.9*fov(3), ...
    'timeBwProduct', rfTB, ...
    'system', sys, ...
    'use', 'excitation');
gzSS = trap4ge(gzSS, CRT, sys);
gzSS.delay = rf.delay - gzSS.riseTime; % sync RF pulse onset with slice-select gradient
gzSSR = trap4ge(gzSSR, CRT, sys);
end
