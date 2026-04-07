%% 3D-EPI sequence for arbitrary sampling along the 2 PE directions
% Rex Fung

%% Define experiment parameters
run('params.m');

%% Path and options
seqname = 'anyEPI_bad';

%% Excitation pulse
% Target a slightly thinner slice to alleviate aliasing
[rf, gzSS, gzSSR, delay] = mr.makeSincPulse(alpha/180*pi,...
                                     'duration',rfDur,...
                                     'sliceThickness',0.9*fov(3),...
                                     'timeBwProduct', rfTB, ...
                                     'system',sys,...
                                     'use','excitation');
gzSS = trap4ge(gzSS,CRT,sys);
gzSS.delay = rf.delay - gzSS.riseTime; % Sync up rf pulse and slice select gradient
gzSSR = trap4ge(gzSSR,CRT,sys);

%% Fat-sat
% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ... % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);

% Create pulseq object
% Try to account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip_ang = fatsat.flip/180*pi;
flipAssumed = abs(sum(rfp));
rfsat = mr.makeArbitraryRf(rfp, ...
    flip_ang*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'system', sys, ...
    'use', 'saturation');
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = -fatOffresFreq; % Hz

%% Generate temporally incoherent sampling masks, partitions, and schedule
omegas = zeros(Ny, Nz, Nframes);
parts = zeros(Ny, Nz, Nframes);
schedules = zeros(Nframes, Nshots, ETL, 2);
for frame = 1:Nframes
    rng(frame);
    % omega = caipi_sample([Ny, Nz], [3, 2], 3);
    omega = pd_sample([Ny, Nz], R);
    % weights = gen_gaussian_pdf([Ny, Nz], [Ny, Nz] ./ 6);
    % omega = rand_sample([Ny, Nz], R, weights);

    omegas(:,:,frame) = omega; 
    [schedules(frame,:,:,:), parts(:,:,frame)] = mask2epi(omega, ETL, Nshots);
end

%% Infer maximum ky, kz jumps from partitions
max_ky_step = 0; max_kz_step = 0;
for frame = 1:Nframes
    for shot = 1:Nshots
        max_ky_step = max(max_ky_step, max(diff(schedules(frame,shot,:,1))));
        max_kz_step = max(max_kz_step, max(diff(schedules(frame,shot,:,2))));
    end
end

%% Define readout gradients and ADC event
% The Pulseq toolbox really shines here!

% Define k-space spacing for fully-sampled data
deltak = 1./fov;

% Start with the blips. Ensure long enough to support the largest blips
gyBlip = mr.makeTrapezoid('y', sys, 'area', max_ky_step*deltak(2));
gyBlip = mr.scaleGrad(gyBlip, 1/max_ky_step, sys);
gyBlip = trap4ge(gyBlip,CRT,sys);
gzBlip = mr.makeTrapezoid('z', sys, 'area', max_kz_step*deltak(3));
gzBlip = mr.scaleGrad(gzBlip, 1/max_kz_step, sys);
gzBlip = trap4ge(gzBlip,CRT,sys);

% Match blip durations
if mr.calcDuration(gyBlip) > mr.calcDuration(gzBlip) % biggest blip in y
    maxBlipArea = max_ky_step*deltak(2);
    blipDuration = mr.calcDuration(gyBlip);
    gzBlip = mr.makeTrapezoid('z', sys, 'area', max_kz_step*deltak(3), 'duration', blipDuration);
    gzBlip = mr.scaleGrad(gzBlip, 1/max_kz_step, sys);
else % biggest blip in z
    maxBlipArea = max_kz_step*deltak(3);
    blipDuration = mr.calcDuration(gzBlip);
    gyBlip = mr.makeTrapezoid('y', sys, 'area', max_ky_step*deltak(2));
    gyBlip = mr.scaleGrad(gyBlip, 1/max_ky_step, sys);
end
gyBlip = trap4ge(gyBlip, CRT, sys);
gzBlip = trap4ge(gzBlip, CRT, sys);

% Readout trapezoid
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;  % to ensure Nyquist sampling
gro = mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea);
gro = trap4ge(gro,CRT,sys);

% Circularly shift gro waveform to contain blips within each block
% [gro1, gro2] = mr.splitGradientAt(gro, blipDuration/2);
% gro2.delay = 0;
% gro1.delay = gro2.shape_dur;
% gro = mr.addGradients({gro2, mr.scaleGrad(gro1, -1)}, sys);
% gro1.delay = 0; % This piece is necessary at the very beginning of the readout

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
Nfid = floor(Tread/dwell/4)*4;
adc = mr.makeAdc(Nfid, 'Dwell', dwell);
adc.delay = blipDuration/2;

% Delay blips so they play after adc stops
% gyBlip.delay = Tread;
% gzBlip.delay = Tread;

% Split blips into up and down parts
[gyBU, gyBD] = mr.splitGradientAt(gyBlip, mr.calcDuration(gyBlip)/2);
gyBU.shape_dur = mr.calcDuration(gyBlip)/2;
gyBD.shape_dur = mr.calcDuration(gyBlip)/2;
gyBU.tt(2) = mr.calcDuration(gyBlip)/2;
gyBD.tt(2) = mr.calcDuration(gyBlip)/2;
gyBU.delay = mr.calcDuration(gro) - mr.calcDuration(gyBU);
gyBD.delay = 0;
[gzBU, gzBD] = mr.splitGradientAt(gzBlip, mr.calcDuration(gzBlip)/2);
gzBU.shape_dur = mr.calcDuration(gzBlip)/2;
gzBD.shape_dur = mr.calcDuration(gzBlip)/2;
gzBU.tt(2) = mr.calcDuration(gzBlip)/2;
gzBD.tt(2) = mr.calcDuration(gzBlip)/2;
gzBU.delay = mr.calcDuration(gro) - mr.calcDuration(gzBU);
gzBD.delay = 0;

% Prephasers (Make duration long enough to support all 3 directions)
gxPre = trap4ge(mr.makeTrapezoid('x',sys,'Area',-(Nx*deltak(1) + maxBlipArea)/2),CRT,sys);
gyPre = trap4ge(mr.makeTrapezoid('y',sys,'Area',-Ny/2*deltak(2)),CRT,sys);
gzPre = trap4ge(mr.makeTrapezoid('z',sys,'Area',-Nz/2*deltak(3)),CRT,sys);

% Spoilers (conventionally only in x and z because ??, might as well do so in y)
gxSpoil = trap4ge(mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*NcyclesSpoil),CRT,sys);
gySpoil = trap4ge(mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)*NcyclesSpoil),CRT,sys);
gzSpoil = trap4ge(mr.scaleGrad(...
    mr.makeTrapezoid('z', sys, 'Area', Nz*deltak(3)*(NcyclesSpoil + 0.5)),...
    NcyclesSpoil/(NcyclesSpoil + 0.5)),CRT,sys);

%% Calculate delay to achieve desired TE
minTE = 0.5*mr.calcDuration(rf)...
      + mr.calcDuration(gzSSR)...
      + max([mr.calcDuration(gxPre), mr.calcDuration(gyPre), mr.calcDuration(gzPre)])...
      + (ETL/2 - 0.5) * mr.calcDuration(gro);
if TE >= minTE
    TEdelay = floor((TE - minTE)/sys.blockDurationRaster) * sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TE (%d) exceeds prescribed TE (%d)',...
                    minTE, TE))
    TEdelay = 0;
end

%% Calculate delay to achieve desired TR
minTR = mr.calcDuration(rfsat)...
      + max([mr.calcDuration(gxSpoil), mr.calcDuration(gzSpoil)])...
      + max([mr.calcDuration(rf), mr.calcDuration(gzSS)])...
      + mr.calcDuration(gzSSR)...
      + TEdelay...
      + max([mr.calcDuration(gxPre), mr.calcDuration(gyPre), mr.calcDuration(gzPre)])...
      + ETL * mr.calcDuration(gro)...
      + max([mr.calcDuration(gxSpoil), mr.calcDuration(gySpoil), mr.calcDuration(gzSpoil)]);
if TR >= minTR
    TRdelay = floor((TR - minTR)/sys.blockDurationRaster)*sys.blockDurationRaster;
else
    warning(sprintf('Minimum achievable TR (%d) exceeds prescribed TR (%d)',...
                    minTR, TR))
    TRdelay = 0;
end

%% Assemble sequence
% manually set to 0 to avoid annoying warnings. 
% Shouldn't be a problem since I don't have back-to-back blocks with adc.
sys.adcDeadTime = 0;

seq = mr.Sequence(sys);

% log the sequence of k-space locations sampled (ky and kz)
samp_log = zeros(Nframes, Nshots*ETL, 2);

% RF spoiling trackers
rf_count = 1;
rf_phase = rf_phase_0;

for frame = 1:Nframes
    fprintf('Writing frame %d\n', frame);

    for shot = 1:Nshots
        % Label the first block in each "unique" section with TRID (see Pulseq on GE manual)
        TRID = 1;

        % Fat-sat
        seq.addBlock(rfsat, mr.makeLabel('SET','TRID',TRID));
        seq.addBlock(gxSpoil, gySpoil, gzSpoil);

        % RF spoiling
        rf_phase = mod(0.5 * rf_phase_0 * rf_count^2, 360.0);
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_count = rf_count + 1;

        % Slab-selective RF excitation + rephase
        seq.addBlock(rf, gzSS);
        seq.addBlock(gzSSR);

        % TE delay
        if TE > minTE
            seq.addBlock(mr.makeDelay(TEdelay));
        end

        % Load in sampling schedule
        schedule = squeeze(schedules(frame,shot,:,:));
        y_locs = schedule(:,1);
        z_locs = schedule(:,2);

        % Move to first location
        gzPreTmp = mr.scaleGrad(gzPre, (z_locs(1) - Nz/2 - 1)/(-Nz/2));
        gyPreTmp = mr.scaleGrad(gyPre, (y_locs(1) - Ny/2 - 1)/(-Ny/2));
        seq.addBlock(gxPre, gyPreTmp, gzPreTmp);

        % Begin ky encoding
        % first line
        gyBlipTmp = mr.scaleGrad(gyBU, y_locs(2) - y_locs(1));
        gzBlipTmp = mr.scaleGrad(gzBU, z_locs(2) - z_locs(1));
        seq.addBlock(adc, gro, gyBlipTmp, gzBlipTmp);

        % echo train
        for echo = 2:(length(y_locs) - 1)
            gyBlipTmp = mr.addGradients(...
                        {mr.scaleGrad(gyBD, y_locs(echo) - y_locs(echo-1)),...
                        mr.scaleGrad(gyBU, y_locs(echo+1) - y_locs(echo))},...
                        sys);
            gzBlipTmp = mr.addGradients(...
                        {mr.scaleGrad(gzBD, z_locs(echo) - z_locs(echo-1)),...
                        mr.scaleGrad(gzBU, z_locs(echo+1) - z_locs(echo))},...
                        sys);

            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(echo-1)),...
                gyBlipTmp,...
                gzBlipTmp...
                );
        end

        % last line
        gyBlipTmp = mr.scaleGrad(gyBD, y_locs(echo+1) - y_locs(echo));
        gzBlipTmp = mr.scaleGrad(gzBD, z_locs(echo+1) - z_locs(echo));
        seq.addBlock(adc, mr.scaleGrad(gro, (-1)^echo), gyBlipTmp, gzBlipTmp);

        % End ky encoding

        % Gx, Gz spoilers. Gy rewinder gradients
        seq.addBlock(gxSpoil, ...
            mr.scaleGrad(gySpoil, (gySpoil.area - (y_locs(end) - Ny/2)*deltak(2))/gySpoil.area), ...
            mr.scaleGrad(gzSpoil, (gzSpoil.area - (z_locs(end) - Nz/2)*deltak(3))/gzSpoil.area));

        % Achieve desired TR
        if TR > minTR
            seq.addBlock(mr.makeDelay(TRdelay));
        end
    end
end

%% Check sequence timing
[ok, error_report] = seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else        
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n'); 
end

%% Save sampling locations for gridding
save('samp_locs.mat', 'schedules', 'parts', '-v7.3');

%% Plot in pulseq (no figure call needed)
seq.plot('timeRange', [0 2*max(minTR, TR)], 'stacked', 1);

%% Write to .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', seqname);
fn_seq = strcat(seqname, '.seq');
seq.write(fn_seq);

%% Interpret to GE via TOPPE
% Define hardware parameters for MR750 scanner
psd_rf_wait = 150e-6;  % RF-gradient delay (s)
psd_grd_wait = 120e-6; % ADC-gradient delay (s)
b1_max = 0.25;         % Gauss
g_max = 5;             % Gauss/cm
slew_max = 20;         % Gauss/cm/ms
sysPGE2 = pge2.opts(psd_rf_wait, psd_grd_wait, b1_max, g_max, slew_max, 'xrm');

% Write to GE compatible file
pislquant = 10; % Number of ADC events at start of scan for receive gain calibration
% PNSwt = [0.8 1 0.7]; % PNS channel/direction weights, less for L/R & S/I directions
PNSwt = [0 0 0]; % for phantom
pge2.seq2ge(fn_seq, sysPGE2, pislquant, PNSwt);

%% Plot in pge2
ceq = seq2ceq(seq);
figure;
S = pge2.plot(ceq, sysPGE2, 'blockRange', [1 10], 'rotate', false, 'interpolate', true, 'wt', [0.8 1 0.7]);

%% Plot only gradients in pge2
figure; plotPGE2grads(ceq, sysPGE2, 'blockRange', [1 10], 'showBlocks', true, 'rotate', false, 'interpolate', true, 'wt', [0.8 1 0.7]);
return;

%% Plot trajectories stringing together samples (takes a while)
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

figure;
hold on;

tic;
% Loop over each excitation
for i = 1:Nshots
    t_start = t_excitation(i);
    
    % Use end of current excitation or end of ADC window
    if i < length(t_excitation)
        t_end = t_excitation(i+1);
    else
        t_end = t_adc(end) + 1e-6;  % Small buffer
    end

    % Find ADC times within this excitation window
    adc_mask = t_adc >= t_start & t_adc < t_end;
    t_adc_segment = t_adc(adc_mask);

    % Map ADC times to indices in t_ktraj
    adc_indices = arrayfun(@(t) find(abs(t_ktraj - t) < 1e-9, 1, 'first'), t_adc_segment);

    % Extract corresponding k-space trajectory points
    ktraj_segment = ktraj(:, adc_indices);

    % Plot this segment as a separate line
    if ~isempty(ktraj_segment)
        plot(ktraj_segment(2,:), ktraj_segment(3,:), 'b', 'LineWidth', 1.5);
    end
end

samps = 1:(length(ktraj_adc)/Nframes);
plot(ktraj_adc(2,samps), ktraj_adc(3,samps),'r.', 'MarkerSize', 12); % plot the sampling points

hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom');
uistack(hy, 'bottom');

axis equal;
title(sprintf('3D-EPI trajectory. R = %d', round(R)), 'FontSize', 18);
xlabel('k_y (m^{-1})', 'FontSize', 18); ylabel('k_z (m^{-1})', 'FontSize', 18);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);

toc;

%% Plot sampling masks on grid
[kys, kzs] = meshgrid([-Ny*deltak(2)/2:deltak(2):Ny*deltak(2)/2], [-Nz*deltak(3)/2:deltak(3):Nz*deltak(3)/2]);

figure; hold on;
hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom');
uistack(hy, 'bottom');
plot(kys(:), kzs(:), 'k.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 12);
plot(ktraj_adc(2,samps), ktraj_adc(3,samps),'r.', 'MarkerSize', 12); % plot the sampling points

axis equal;
title(sprintf('2D sampling mask. R = %d', round(R)), 'FontSize', 18);
xlabel('k_y (m^{-1})', 'FontSize', 18); ylabel('k_z (m^{-1})', 'FontSize', 18);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);

return;

%% Check for forbidden gradient frequencies
check_grad_acoustics(reshape([S.gx.signal S.gy.signal S.gz.signal], [length(S.gx.signal), 1, 3])/100, 'xrm');
return;

%% Detailed sequence report
% Slow but useful for testing during development,
% e.g., for the real TE, TR or for staying within slew rate limits
rep = seq.testReport;
fprintf([rep{[1:9, 11:end]}]); % print report except warnings
return;

%% Plot point spread function in y-z space
omega = omegas(:,:,1);
psf = ifftshift(ifft2(fftshift(omega)));
[Y, Z] = meshgrid(res(2)*(-Ny/2:Ny/2-1), res(3)*(-Nz/2:Nz/2-1));
figure; surf(Y.', Z.', abs(psf), 'FaceColor', 'interp');
xlim([-Ny*res(2)/2, Ny*res(2)/2]); ylim([-Nz*res(3)/2, Nz*res(3)/2]);
zlim([0, max(abs(psf(:)))])
xlabel('y (m)', 'FontSize', 18); ylabel('z (m)', 'FontSize', 18);
zlabel('magnitude (a.u.)', 'FontSize', 18);
title('Corresponding point spread function in y-z space', 'FontSize', 18);
zticks(linspace(0,max(abs(psf(:))),5)); zticklabels(linspace(0,1,5));
return;