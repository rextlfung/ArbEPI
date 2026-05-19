function EPIcal()
%% EPIcal  Calibration sequence for EPI ghost correction and receiver gain tuning.
%
% Structurally mirrors ArbEPI but with no ky/kz blips, producing unencoded
% readout lines at k-space center. Acquisition flow:
%   1. Dummy shots (no ADC) — drive magnetization to steady state
%   2. Real shots (ADC on)  — acquire for EPI ghost correction
%
% Requires samp_locs.mat from a prior ArbEPI run (used to match readout
% gradient design so the calibration trajectory matches the imaging trajectory).
%
% Outputs: EPIcal.seq  (Pulseq format)
%          EPIcal.pge  (GE TOPPE format)
%          kxoe<Nx>.mat (odd/even echo k-space trajectories for ghost correction)

%% Parameters
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;
seqname = 'EPIcal';

%% Load EPI schedule to match readout gradient design
load(fullfile(outputDir, 'samp_locs.mat'), 'schedules');

%% Excitation pulse (identical to ArbEPI)
[rf, gzSS, gzSSR] = make_excitation_pulse(fa, rfDur, rfTB, fov, sys, CRT);

%% Fat-sat pulse (identical to ArbEPI)
rfsat = make_fatsat_rf(fatsat, sys, fatOffresFreq);

%% Readout gradients — sized to match ArbEPI so calibration trajectory applies
max_ky_step = max(abs(diff(schedules(:,:,:,1), 1, 3)), [], 'all');
max_kz_step = max(abs(diff(schedules(:,:,:,2), 1, 3)), [], 'all');
rg = make_readout_grads(max_ky_step, max_kz_step, Nx, fov, dwell, sys, CRT);

%% Prephasers and spoilers (identical to ArbEPI)
[gxPre, gyPre, gzPre]       = make_prephasers(Nx, Ny, Nz, fov, sys, CRT);
[gxSpoil, gySpoil, gzSpoil] = make_spoilers(Nx, Ny, Nz, fov, NcyclesSpoil, sys, CRT);

%% Calculate TE and TR delays (identical to ArbEPI)
[TEdelay, TRdelay, minTE, minTR] = calc_te_tr_delays( ...
    rf, rfsat, gzSS, gzSSR, gxPre, gyPre, gzPre, rg.gro, ...
    gxSpoil, gySpoil, gzSpoil, ETL, TE, TR, sys);

%% Assemble sequence
sys.adcDeadTime = 0; % suppress warnings; no back-to-back ADC blocks
seq = mr.Sequence(sys);

rf_count = 1;

for shot = -Ndummyshots+1:Nshots
    isDummy = shot < 1;
    TRID = 1 + (~isDummy); % TRID 1 = dummy, TRID 2 = real (see Pulseq on GE manual)

    % Fat-sat
    seq.addBlock(rfsat, mr.makeLabel('SET', 'TRID', TRID));
    seq.addBlock(gxSpoil, gzSpoil);

    % RF spoiling (quadratic phase cycling)
    rf_phase = mod(0.5 * rf_phase_0 * rf_count^2, 360.0);
    rf.phaseOffset = rf_phase/180*pi;
    rg.adc.phaseOffset = rf_phase/180*pi;
    rf_count = rf_count + 1;

    % Slab-selective excitation + slice-select rephaser
    seq.addBlock(rf, gzSS);
    seq.addBlock(gzSSR);

    % TE padding delay
    if TE > minTE
        seq.addBlock(mr.makeDelay(TEdelay));
    end

    % Prephase to k-space center (no ky/kz encoding — blips will be zero)
    seq.addBlock(gxPre, mr.scaleGrad(gyPre, 0), mr.scaleGrad(gzPre, 0));

    % EPI readout — same waveform as ArbEPI, blips scaled to zero
    seq.addBlock(rg.gro1);
    echo = 0; % initialize so (-1)^echo is valid when ETL=1 (empty loop)
    for echo = 1:(ETL - 1)
        if ~isDummy
            seq.addBlock(rg.adc, mr.scaleGrad(rg.gro, (-1)^(echo-1)), ...
                mr.scaleGrad(rg.gyBlip, 0), mr.scaleGrad(rg.gzBlip, 0));
        else
            seq.addBlock(mr.scaleGrad(rg.gro, (-1)^(echo-1)), ...
                mr.scaleGrad(rg.gyBlip, 0), mr.scaleGrad(rg.gzBlip, 0));
        end
    end
    % Last echo line
    if ~isDummy
        seq.addBlock(rg.adc, mr.scaleGrad(rg.gro2, (-1)^echo));
    else
        seq.addBlock(mr.scaleGrad(rg.gro2, (-1)^echo));
    end

    % Spoilers: x dephase, y/z no net phase (no encoding was applied)
    seq.addBlock(gxSpoil, mr.scaleGrad(gySpoil, 0), mr.scaleGrad(gzSpoil, 1));

    % TR padding delay
    if TR > minTR
        seq.addBlock(mr.makeDelay(TRdelay));
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

%% Save odd/even echo k-space trajectories for ghost correction
[ktraj_adc, ~, ~, ~, ~, ~] = seq.calculateKspacePP();
kxo = ktraj_adc(1, 1:rg.Nfid);
kxe = ktraj_adc(1, rg.Nfid+1:rg.Nfid*2);
save(fullfile(outputDir, sprintf('kxoe%d.mat', Nx)), 'kxo', 'kxe', '-v7.3');

%% Write Pulseq .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
seq.write(fullfile(outputDir, [seqname '.seq']));

%% Write GE TOPPE .pge file
sysPGE2 = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, 'xrm');
write_to_ge(seq, fullfile(outputDir, seqname), sysPGE2, PNSwt, pislquant);

end % EPIcal
