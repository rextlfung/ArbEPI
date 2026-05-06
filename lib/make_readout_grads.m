function rg = make_readout_grads(max_ky_step, max_kz_step, Nx, fov, dwell, sys, CRT)
%% make_readout_grads  Create EPI readout gradients, blips, and ADC event.
%   max_ky_step: largest ky blip step (in k-space index units)
%   max_kz_step: largest kz blip step (in k-space index units)
%   Nx:          readout matrix size
%   fov:         field of view [x y z] (m)
%   dwell:       ADC sample dwell time (s)
%   sys:         Pulseq system struct (mr.opts)
%   CRT:         common raster time for GE compatibility (s)
%
% Returns struct rg with fields:
%   gro, gro1, gro2  - readout trapezoids (full, first piece, last piece)
%   gyBlip, gzBlip   - unit blip trapezoids (scale at call site)
%   adc              - ADC event
%   Tread            - flat-top readout duration (s)
%   Nfid             - number of ADC samples
%   blipDuration     - blip duration = time between echo lines (s)
%   deltak           - k-space spacing [kx ky kz] (1/m)
%   maxBlipArea      - area of the dominant blip axis (1/m)

deltak = 1./fov;

% Size blips to support the largest steps across all frames/shots.
% Blips are stored at unit amplitude and scaled at assembly time.
gyBlip = mr.makeTrapezoid('y', sys, 'area', max_ky_step*deltak(2));
gyBlip = mr.scaleGrad(gyBlip, 1/max_ky_step, sys);
gyBlip = trap4ge(gyBlip, CRT, sys);
gzBlip = mr.makeTrapezoid('z', sys, 'area', max_kz_step*deltak(3));
gzBlip = mr.scaleGrad(gzBlip, 1/max_kz_step, sys);
gzBlip = trap4ge(gzBlip, CRT, sys);

% Match blip durations so they always fit within one readout block
if mr.calcDuration(gyBlip) > mr.calcDuration(gzBlip) % y blip is longer
    maxBlipArea  = max_ky_step * deltak(2);
    blipDuration = mr.calcDuration(gyBlip);
    gzBlip = mr.makeTrapezoid('z', sys, 'area', max_kz_step*deltak(3), 'duration', blipDuration);
    gzBlip = mr.scaleGrad(gzBlip, 1/max_kz_step, sys);
else % z blip is longer
    maxBlipArea  = max_kz_step * deltak(3);
    blipDuration = mr.calcDuration(gzBlip);
    gyBlip = mr.makeTrapezoid('y', sys, 'area', max_ky_step*deltak(2), 'duration', blipDuration);
    gyBlip = mr.scaleGrad(gyBlip, 1/max_ky_step, sys);
end

% Readout trapezoid sized for ramp-sampling and to contain the blip area
systmp = sys;
systmp.maxGrad = deltak(1)/dwell; % enforce Nyquist sampling on flat top
gro = trap4ge(mr.makeTrapezoid('x', systmp, 'area', Nx*deltak(1) + maxBlipArea), CRT, sys);

% Circularly shift gro so blips are contained within each Pulseq block:
% the waveform is split and reassembled as [second half, -first half].
[gro1, gro2] = mr.splitGradientAt(gro, blipDuration/2);
gro2.delay = 0;
gro1.delay = gro2.shape_dur;
gro = mr.addGradients({gro2, mr.scaleGrad(gro1, -1)}, sys);
gro1.delay = 0; % leading piece played before the first echo line

% ADC event covering the flat-top portion
Tread = mr.calcDuration(gro) - blipDuration;
Nfid  = round(Tread/dwell/4)*4; % round to multiple of 4 for GE
adc   = mr.makeAdc(Nfid, 'Dwell', dwell);

% Delay blips to play after the ADC window closes
gyBlip.delay = Tread;
gzBlip.delay = Tread;

rg.gro          = gro;
rg.gro1         = gro1;
rg.gro2         = gro2;
rg.gyBlip       = gyBlip;
rg.gzBlip       = gzBlip;
rg.adc          = adc;
rg.Tread        = Tread;
rg.Nfid         = Nfid;
rg.blipDuration = blipDuration;
rg.deltak       = deltak;
rg.maxBlipArea  = maxBlipArea;
end
