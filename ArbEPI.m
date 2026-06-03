function ArbEPI(omegas)
%% ArbEPI  Generate 3D-EPI sequence with an arbitrary 2D sampling mask.
%   omegas: Ny x Nz x Nframes logical sampling mask
%           Generate with: omegas = gen_sampling_masks(R);
%
% Outputs: ArbEPI.seq  (Pulseq format)
%          ArbEPI.pge  (GE TOPPE format)
%          samp_locs.mat (schedules and partition map for reconstruction)

%% Parameters
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;
seqname = 'ArbEPI';

%% Excitation pulse
[rf, gzSS, gzSSR] = make_excitation_pulse(fa, rfDur, rfTB, fov, sys, CRT);

%% Fat-sat pulse
rfsat = make_fatsat_rf(fatsat, sys, fatOffresFreq);

%% Generate EPI sampling schedule from mask
[schedules, parts] = compute_schedules(omegas, ETL, Nshots);

% Infer maximum ky and kz blip steps across all frames and shots
max_ky_step = max(abs(diff(schedules(:,:,:,1), 1, 3)), [], 'all');
max_kz_step = max(abs(diff(schedules(:,:,:,2), 1, 3)), [], 'all');

%% Define readout gradients and ADC event
rg = make_readout_grads(max_ky_step, max_kz_step, Nx, fov, dwell, sys, CRT);

%% Prephasers and spoilers
[gxPre, gyPre, gzPre]       = make_prephasers(Nx, Ny, Nz, fov, sys, CRT);
[gxSpoil, gySpoil, gzSpoil] = make_spoilers(Nx, Ny, Nz, fov, NcyclesSpoil, sys, CRT);

%% Calculate delays to achieve desired TE and TR
[TEdelay, TRdelay, minTE, minTR] = calc_te_tr_delays( ...
    rf, rfsat, gzSS, gzSSR, gxPre, gyPre, gzPre, rg.gro, ...
    gxSpoil, gySpoil, gzSpoil, ETL, TE, TR, sys);

%% Assemble sequence
sys.adcDeadTime = 0; % suppress warnings; no back-to-back ADC blocks
seq = mr.Sequence(sys);

rf_count = 1;

for frame = 1:Nframes
    fprintf('Writing frame %d\n', frame);

    for shot = 1:Nshots
        % Fat-sat (label first block in each unique section with TRID for GE)
        seq.addBlock(rfsat, mr.makeLabel('SET', 'TRID', 1));
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

        % Load this shot's k-space locations
        schedule = squeeze(schedules(frame, shot, :, :));
        y_locs = schedule(:, 1);
        z_locs = schedule(:, 2);

        % Move to first k-space location
        gyPreTmp = mr.scaleGrad(gyPre, (y_locs(1) - Ny/2 - 1)/(-Ny/2));
        gzPreTmp = mr.scaleGrad(gzPre, (z_locs(1) - Nz/2 - 1)/(-Nz/2));
        seq.addBlock(gxPre, gyPreTmp, gzPreTmp);

        % EPI readout — zip through k-space with alternating readout polarity
        seq.addBlock(rg.gro1);
        echo = 0; % initialize so (-1)^echo is valid when ETL=1 (empty loop)
        for echo = 1:(length(y_locs) - 1)
            seq.addBlock(rg.adc, mr.scaleGrad(rg.gro, (-1)^(echo-1)), ...
                mr.scaleGrad(rg.gyBlip, y_locs(echo+1) - y_locs(echo)), ...
                mr.scaleGrad(rg.gzBlip, z_locs(echo+1) - z_locs(echo)));
        end
        % Last echo line (no blip needed)
        seq.addBlock(rg.adc, mr.scaleGrad(rg.gro2, (-1)^echo));

        % Spoilers: x/z dephase, y rewind to center
        seq.addBlock(gxSpoil, ...
            mr.scaleGrad(gySpoil, -((y_locs(end) - Ny/2)*rg.deltak(2))/gySpoil.area), ...
            mr.scaleGrad(gzSpoil,  (gzSpoil.area - (z_locs(end) - Nz/2)*rg.deltak(3))/gzSpoil.area));

        % TR padding delay
        if TR > minTR
            seq.addBlock(mr.makeDelay(TRdelay));
        end
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

%% Save sampling locations for reconstruction
save(fullfile(outputDir, 'samp_locs.mat'), 'schedules', 'parts', '-v7.3');

%% Write Pulseq .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
seq.write(fullfile(outputDir, [seqname '.seq']));

%% Write GE TOPPE .pge file
sysPGE2 = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, 'xrm');
write_to_ge(seq, fullfile(outputDir, seqname), sysPGE2, PNSwt, pislquant);

%% K-space trajectory and sampling mask for frame 1
blocks_per_frame = round(length(seq.blockDurations) / Nframes);
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation] = ...
    seq.calculateKspacePP('blockRange', [1, blocks_per_frame]);

deltak = 1./fov;

% Sampling mask
[kys, kzs] = meshgrid(-Ny*deltak(2)/2:deltak(2):Ny*deltak(2)/2, ...
    -Nz*deltak(3)/2:deltak(3):Nz*deltak(3)/2);
figure; hold on;
hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom'); uistack(hy, 'bottom');
plot(kys(:), kzs(:), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 16);
plot(ktraj_adc(2,:), ktraj_adc(3,:), 'r.', 'MarkerSize', 16);
axis equal; set(gca, 'FontSize', 20);
title(sprintf('2D sampling mask. R = %d', round(R)), 'FontSize', 40);
xlabel('k_y (cycles/m)', 'FontSize', 30); ylabel('k_z (cycles/m)', 'FontSize', 30);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);

% Trajectory
figure; hold on;
for i = 1:Nshots
    t_start = t_excitation(i);
    if i < length(t_excitation)
        t_end = t_excitation(i+1);
    else
        t_end = t_adc(end) + 1e-6;
    end
    adc_mask = t_adc >= t_start & t_adc < t_end;
    t_adc_segment = t_adc(adc_mask);
    adc_indices = arrayfun(@(t) find(abs(t_ktraj - t) < 1e-9, 1, 'first'), t_adc_segment);
    ktraj_segment = ktraj(:, adc_indices);
    if ~isempty(ktraj_segment)
        plot(ktraj_segment(2,:), ktraj_segment(3,:), 'b', 'LineWidth', 1.5);
    end
end
plot(ktraj_adc(2,:), ktraj_adc(3,:), 'r.', 'MarkerSize', 16);
hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom'); uistack(hy, 'bottom');
axis equal; set(gca, 'FontSize', 20);
title(sprintf('3D-EPI trajectory. ETL = %d', round(ETL)), 'FontSize', 40);
xlabel('k_y (cycles/m)', 'FontSize', 30); ylabel('k_z (cycles/m)', 'FontSize', 30);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);

%% Pulse sequence diagram
seq.plot('timeRange', [0, TR], 'stacked', 1, 'showBlocks', 1);
set(findall(gcf, 'Type', 'line'), 'LineWidth', 2);

end % ArbEPI

%% Local helper: run mask2epi for each frame
function [schedules, parts] = compute_schedules(omegas, ETL, Nshots)
[Ny, Nz, Nframes] = size(omegas);
schedules = zeros(Nframes, Nshots, ETL, 2);
parts     = zeros(Ny, Nz, Nframes);
for frame = 1:Nframes
    [schedules(frame,:,:,:), parts(:,:,frame)] = mask2epi(omegas(:,:,frame), ETL, Nshots);
end
end
