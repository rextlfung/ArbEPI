% Calibration scan sequence for EPI
%
% This short sequence first excites the volume to steady state, then 
% acquires many readout lines without Gy and Gz blips to:
% 1. Allow the scanner to tune receiver gains
% 2. Collect data used for EPI ghost correciton

%% Define experiment parameters
run('params.m');

%% Path and options
seqname = 'EPIcal';

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

%% Load in 3D EPI partitions and schedules
load('samp_locs.mat');

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

% Readout trapezoid
systmp = sys;
systmp.maxGrad = deltak(1)/dwell;  % to ensure Nyquist sampling
gro = mr.makeTrapezoid('x', systmp, 'Area', Nx*deltak(1) + maxBlipArea);
gro = trap4ge(gro,CRT,sys);

% Circularly shift gro waveform to contain blips within each block
[gro1, gro2] = mr.splitGradientAt(gro, blipDuration/2);
gro2.delay = 0;
gro1.delay = gro2.shape_dur;
gro = mr.addGradients({gro2, mr.scaleGrad(gro1, -1)}, sys);
gro1.delay = 0; % This piece is necessary at the very beginning of the readout

% ADC event
Tread = mr.calcDuration(gro) - blipDuration;
Nfid = floor(Tread/dwell/4)*4;
adc = mr.makeAdc(Nfid, 'Dwell', dwell);

% Delay blips so they play after adc stops
gyBlip.delay = Tread;
gzBlip.delay = Tread;

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

% RF spoiling trackers
rf_count = 1;
rf_phase = rf_phase_0;

for shot = -Ndummyshots+1:Nshots
    % Label the first block in each "unique" section with TRID (see Pulseq on GE manual)
    if shot < 1
        TRID = 1;
    else
        TRID = 2;
    end

    % Fat-sat
    seq.addBlock(rfsat, mr.makeLabel('SET','TRID',TRID));
    seq.addBlock(gxSpoil, gzSpoil);

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

    % Move to first location
    gzPreTmp = mr.scaleGrad(gzPre, 0);
    gyPreTmp = mr.scaleGrad(gyPre, 0);
    seq.addBlock(gxPre, gyPreTmp, gzPreTmp);

    % Begin ky encoding
    % Zip through k-space with EPI trajectory
    seq.addBlock(gro1);
    for echo = 1:(ETL - 1)
        if shot > 0
            seq.addBlock(adc, mr.scaleGrad(gro, (-1)^(echo-1)),...
                mr.scaleGrad(gyBlip, 0),...
                mr.scaleGrad(gzBlip, 0)...
                );
        else
            seq.addBlock(mr.scaleGrad(gro, (-1)^(echo-1)),...
                mr.scaleGrad(gyBlip, 0),...
                mr.scaleGrad(gzBlip, 0)...
                );
        end
    end

    % Last line
    if shot > 0
        seq.addBlock(adc, mr.scaleGrad(gro2, (-1)^echo));
    else
        seq.addBlock(mr.scaleGrad(gro2, (-1)^echo));
    end

    % End ky encoding

    % Gx, Gz spoilers. Gy rewinder gradients
    seq.addBlock(gxSpoil, ...
        mr.scaleGrad(gySpoil, 0), ...
        mr.scaleGrad(gzSpoil, 1));

    % Achieve desired TR
    if TR > minTR
        seq.addBlock(mr.makeDelay(TRdelay));
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

%% Save k-space trajectory for EPI ghost correction
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
kxo = ktraj_adc(1, 1:Nfid);
kxe = ktraj_adc(1, Nfid+1:Nfid*2);

save(sprintf('kxoe%d.mat', Nx),'kxo','kxe','-v7.3');

%% Plot in pulseq
f = seq.plot('timeRange', [0 2*max(minTR, TR)], 'stacked', 1);

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
PNSwt = [0.8 1 0.7]; % PNS channel/direction weights, less for L/R & S/I directions
PNSwt = [0 0 0]; % for phantom
pge2.seq2ge(fn_seq, sysPGE2, pislquant, PNSwt);

%% Check for forbidden gradient frequencies
ceq = seq2ceq(seq);
S = pge2.plot(ceq, sysPGE2, 'blockRange', [1 10], 'rotate', false, 'interpolate', true, 'wt', [0.8 1 0.7]);
check_grad_acoustics(reshape([S.gx.signal S.gy.signal S.gz.signal], [length(S.gx.signal), 1, 3])/100, 'xrm');
return;

%% Detailed sequence report
% Slow but useful for testing during development,
% e.g., for the real TE, TR or for staying within slew rate limits
rep = seq.testReport;
fprintf([rep{[1:9, 11:end]}]); % print report except warnings
return;