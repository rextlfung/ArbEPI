function ceq = write_to_ge(seq, filepath, sysPGE2, PNSwt, pislquant)
%% write_to_ge  Convert Pulseq sequence to GE TOPPE format and write .pge file.
%   seq:       mr.Sequence object
%   filepath:  full output path without extension (e.g. fullfile(outputDir, seqname))
%   sysPGE2:   GE system struct from pge2.opts(...)
%   PNSwt:     PNS direction weights [x y z]
%   pislquant: number of ADC events at scan start for receive gain calibration
%
% Returns ceq (compact sequence representation) for optional downstream use
% (acoustics check, plotting).

ceq = seq2ceq(seq);
figure; S = pge2.plot(ceq, sysPGE2, 'blockRange', [1 10], 'rotate', false, 'interpolate', true, 'wt', PNSwt);

% Check for forbidden gradient frequencies on GE MR750 (xrm)
check_grad_acoustics(reshape([S.gx.signal S.gy.signal S.gz.signal], [length(S.gx.signal), 1, 3])/100, 'xrm', [0, 0]);

pge_params = pge2.check(ceq, sysPGE2, 'wt', PNSwt);
pge2.writeceq(ceq, [filepath '.pge'], 'pislquant', pislquant, 'params', pge_params);
end
