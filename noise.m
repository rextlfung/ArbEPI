%% Noise prescan for noise-prewhitening
% Rex Fung

%% Define experiment parameters
run('params.m');
Ncoils = 32;

%% Path and options
seqname = 'noise';

%% Calculate number of samples needed for noise scan
Nsamples_noise = 20 * Ncoils^2;
Nreps = ceil(Nsamples_noise / Nfid);

%% Calculate delay needed to match TR
% manually set to 0 to avoid annoying warnings. 
% Shouldn't be a problem since I don't have back-to-back blocks with adc.
sys.adcDeadTime = 0;

adc_total_dur = sys.adcDeadTime + mr.calcDuration(adc); % total ADC event duration
pad_duration  = mr.calcDuration(gro) - adc_total_dur;

if pad_duration > 1e-9   % only add if non-trivial
    delay_block = mr.makeDelay(pad_duration);
else
    delay_block = [];
end

%% Assemble sequence using parts from anyEPI.m
seq = mr.Sequence(sys);
for rep = 1:Nreps
    if ~isempty(delay_block)
        % Label the first block in each "unique" section with TRID (see Pulseq on GE manual)
        seq.addBlock(adc, mr.makeLabel('SET','TRID',1));
        seq.addBlock(delay_block);
    else
        seq.addBlock(adc, mr.makeLabel('SET','TRID',1));
    end
end

%% Add dummy rf pulse at the end to make scanner happy
seq.addBlock(rf, mr.makeLabel('SET','TRID',2));

%% Validate sequence
fprintf('Validating sequence...\n');
[ok, errorReport] = seq.checkTiming();
if ok
    fprintf('  Timing check PASSED.\n\n');
else
    fprintf('  Timing check WARNINGS/ERRORS:\n');
    fprintf('    %s\n', errorReport{:});
    fprintf('\n');
end

%% Write to .seq file
fn_seq = strcat(seqname, '.seq');
seq.write(fn_seq);
fprintf('Sequence written to: %s\n', fn_seq);

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
% PNSwt = [0 0 0]; % for phantom
pge2.seq2ge(fn_seq, sysPGE2, pislquant, PNSwt);

%% Quick plot
n_plot = min(10, Nreps);
fprintf('Plotting first %d blocks...\n', n_plot);
seq.plot('TimeRange', [0, n_plot * mr.calcDuration(gro)]);