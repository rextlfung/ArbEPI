function GRE()
%% GRE  T1-weighted gradient echo sequence for sensitivity map estimation.
%
% Acquires a slightly larger FOV than EPI, then crop during reconstruction.
% Outputs: GRE.seq  (Pulseq format)
%          GRE.pge  (GE TOPPE format)

%% Parameters
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;

seqname = 'GRE';
seq = mr.Sequence(sys);

%% Fat-sat (longer pulse than EPI for better fat suppression)
fatsat.flip    = 90;   % degrees
fatsat.slThick = 1e5;  % dummy large value (reduces dead time; no slice selection needed)
fatsat.tbw     = 3;    % time-bandwidth product
fatsat.dur     = 8.0;  % pulse duration (ms)
rfsat = make_fatsat_rf(fatsat, sys, fatOffresFreq);

%% Slab-selective excitation (same pulse as EPI)
[rf, gzSS, gzSSR, ~] = mr.makeSincPulse(fa/180*pi, ...
    'duration', rfDur, ...
    'sliceThickness', 0.9*fov(3), ...
    'timeBwProduct', rfTB, ...
    'system', sys, ...
    'use', 'excitation');
gzSS = trap4ge(gzSS, CRT, sys);
gzSSR = trap4ge(gzSSR, CRT, sys);

%% Readout and phase-encode gradients
deltak = 1./fov_gre;
Tread = Nx_gre*dwell;

gyPre = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', Ny_gre*deltak(2)/2, ...
    'Duration', Tpre), CRT, sys);
gzPre = trap4ge(mr.makeTrapezoid('z', sys, ...
    'Area', Nz_gre*deltak(3)/2, ...
    'Duration', Tpre), CRT, sys);

gxtmp = mr.makeTrapezoid('x', sys, ...
    'Amplitude', Nx_gre*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', -gxtmp.area/2, ...
    'Duration', Tpre), CRT, sys);

adc = mr.makeAdc(Nx_gre, sys, ...
    'Duration', Tread, ...
    'Delay', gxtmp.riseTime);

% Extend flat time to split at end of ADC dead time
gx = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Amplitude', Nx_gre*deltak(1)/Tread, ...
    'FlatTime', Tread + adc.deadTime), CRT, sys);

gzSpoil = mr.makeTrapezoid('z', sys, 'Area', Nx_gre*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeTrapezoid('x', sys, 'Area', Nx_gre*deltak(1)*nCyclesSpoil);

%% Phase-encode step vectors
pe1Steps = ((0:Ny_gre-1) - Ny_gre/2) / Ny_gre * 2;
pe2Steps = ((0:Nz_gre-1) - Nz_gre/2) / Nz_gre * 2;

%% Calculate TE and TR delays
TEmin = max(mr.calcDuration(rf), mr.calcDuration(gzSS))/2 + mr.calcDuration(gzSSR) ...
      + mr.calcDuration(gxPre) + adc.delay + Nx_gre/2*dwell;
delayTE = ceil((TE_gre - TEmin) / seq.gradRasterTime) * seq.gradRasterTime;
TRmin = max(mr.calcDuration(rf), mr.calcDuration(gzSS)) + mr.calcDuration(gzSSR) ...
      + delayTE + mr.calcDuration(gxPre) ...
      + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = ceil((TR_gre - TRmin) / seq.gradRasterTime) * seq.gradRasterTime;

%% Assemble sequence
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC on for receive gain calibration (GE auto-prescan)
% iZ > 0: Image acquisition

rf_count = 1;
lastmsg = [];
for iZ = -NdummyZloops:Nz_gre
    isDummyTR = iZ < 0;

    for ii = 1:length(lastmsg), fprintf('\b'); end
    msg = sprintf('z encode %d of %d ', iZ, Nz_gre);
    fprintf(msg);
    lastmsg = msg;

    % Fat-sat
    TRID = 2 - isDummyTR;
    seq.addBlock(rfsat, mr.makeLabel('SET', 'TRID', TRID));
    seq.addBlock(gzSpoil);

    for iY = 1:Ny_gre
        % Phase-encode steps off during dummy and gain-cal shots
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1, iZ));

        % RF spoiling
        rf_phase = mod(0.5 * rf_phase_0 * rf_count^2, 360.0);
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_count = rf_count + 1;

        seq.addBlock(rf, gzSS);
        seq.addBlock(gzSSR);
        seq.addBlock(mr.makeDelay(delayTE));
        seq.addBlock(gxPre, mr.scaleGrad(gyPre, yStep), mr.scaleGrad(gzPre, zStep));
        if isDummyTR
            seq.addBlock(gx);
        else
            seq.addBlock(gx, adc);
        end
        seq.addBlock(gxSpoil, mr.scaleGrad(gyPre, -yStep), mr.scaleGrad(gzPre, -zStep));
        seq.addBlock(mr.makeDelay(delayTR));
    end
end

%% Check sequence timing
[ok, error_report] = seq.checkTiming;
if ok
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Write Pulseq .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
fn_seq = fullfile(outputDir, [seqname '.seq']);
seq.write(fn_seq);

%% Write GE TOPPE .pge file
sysPGE2 = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, 'xrm');
pge2.seq2ge(fn_seq, sysPGE2, pislquant, PNSwt);


end % GRE
