function rfsat = make_fatsat_rf(fatsat, sys, fatOffresFreq)
%% make_fatsat_rf  Create fat-saturation RF pulse object.
%   fatsat:       struct with fields flip, slThick, tbw, dur
%   sys:          Pulseq system struct (mr.opts)
%   fatOffresFreq: fat off-resonance frequency (Hz)

% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ... % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);

% Create pulseq object.
% makeArbitraryRf scales the pulse as: signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi)
% so we back out the correct flip angle to pass in.
flip_ang = fatsat.flip/180*pi;
flipAssumed = abs(sum(rfp));
rfsat = mr.makeArbitraryRf(rfp, ...
    flip_ang*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'system', sys, ...
    'use', 'saturation');
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = -fatOffresFreq; % Hz
end
