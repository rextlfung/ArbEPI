function noise()
%% noise  Noise prescan for receiver noise covariance estimation.
%
% Acquires ADC data without RF or gradients to measure the noise
% covariance matrix used for noise prewhitening during reconstruction.
%
% Requires samp_locs.mat from a prior ArbEPI run to match readout geometry.
%
% Outputs: noise.seq  (Pulseq format)
%          noise.pge  (GE TOPPE format)

%% Parameters
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;

seqname = 'noise';

%% Rebuild readout gradient objects matching the EPI sequence
load(fullfile(outputDir, 'samp_locs.mat'), 'schedules');
max_ky_step = max(abs(diff(schedules(:,:,:,1), 1, 3)), [], 'all');
max_kz_step = max(abs(diff(schedules(:,:,:,2), 1, 3)), [], 'all');
rg = make_readout_grads(max_ky_step, max_kz_step, Nx, fov, dwell, sys, CRT);
[rf, ~, ~] = make_excitation_pulse(fa, rfDur, rfTB, fov, sys, CRT); 

%% Calculate number of ADC repetitions needed
Nsamples_noise = 20 * Ncoils^2;
Nreps = ceil(Nsamples_noise / rg.Nfid);

%% Calculate delay to pad each block to EPI readout duration
sys.adcDeadTime = 0; % suppress warnings; no back-to-back ADC blocks
adc_total_dur = sys.adcDeadTime + mr.calcDuration(rg.adc);
pad_duration  = mr.calcDuration(rg.gro) - adc_total_dur;

if pad_duration > 1e-9
    delay_block = mr.makeDelay(pad_duration);
else
    delay_block = [];
end

%% Assemble sequence
seq = mr.Sequence(sys);
for rep = 1:Nreps
    seq.addBlock(rg.adc, mr.makeLabel('SET', 'TRID', 1));
    if ~isempty(delay_block)
        seq.addBlock(delay_block);
    end
end

% Dummy RF pulse at end so the scanner recognises a second block type
seq.addBlock(rf, mr.makeLabel('SET', 'TRID', 2));

%% Check sequence timing
fprintf('Validating sequence...\n');
[ok, errorReport] = seq.checkTiming();
if ok
    fprintf('  Timing check PASSED.\n\n');
else
    fprintf('  Timing check WARNINGS/ERRORS:\n');
    fprintf('    %s\n', errorReport{:});
    fprintf('\n');
end

%% Write Pulseq .seq file
fn_seq = fullfile(outputDir, [seqname '.seq']);
seq.write(fn_seq);
fprintf('Sequence written to: %s\n', fn_seq);

%% Write GE TOPPE .pge file
sysPGE2 = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, 'xrm');
pge2.seq2ge(fn_seq, sysPGE2, pislquant, PNSwt);


end % noise
