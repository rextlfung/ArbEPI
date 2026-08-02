%% plot_epi.m  Visualization scripts for ArbEPI sequence and sampling masks.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'lib'));
%
% Run this script after ArbEPI() to produce diagnostic plots.
% Requires the following variables in the workspace:
%   seq, omegas, schedules, fov, res, Ny, Nz, Nshots, Nframes

%% Plot sequence timing in Pulseq viewer
seq.plot('timeRange', [0 max(minTR, TR)], 'stacked', 1);
return;

%% Plot only gradients in pge2
figure; plotPGE2grads(ceq, sysPGE2, 'blockRange', [1 10], 'showBlocks', true, 'rotate', false, 'interpolate', true, 'wt', [0.8 1 0.7]);
return;

%% Plot k-space trajectories (slow — strings together all ADC sample points)
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, ~] = seq.calculateKspacePP();

figure; hold on;
tic;
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

deltak = 1./fov;
samps = 1:(length(ktraj_adc)/Nframes);
plot(ktraj_adc(2,samps), ktraj_adc(3,samps), 'r.', 'MarkerSize', 12);

hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom'); uistack(hy, 'bottom');
axis equal;
title(sprintf('3D-EPI trajectory. R = %d', round(R)), 'FontSize', 18);
xlabel('k_y (m^{-1})', 'FontSize', 18); ylabel('k_z (m^{-1})', 'FontSize', 18);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);
toc;
return;

%% Plot sampling mask on k-space grid
[kys, kzs] = meshgrid(-Ny*deltak(2)/2:deltak(2):Ny*deltak(2)/2, ...
                      -Nz*deltak(3)/2:deltak(3):Nz*deltak(3)/2);

figure; hold on;
hx = plot([-Ny*deltak(2)/2, Ny*deltak(2)/2], [0 0], '-k');
hy = plot([0 0], [-Nz*deltak(3)/2, Nz*deltak(3)/2], '-k');
uistack(hx, 'bottom'); uistack(hy, 'bottom');
plot(kys(:), kzs(:), 'k.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 12);
plot(ktraj_adc(2,samps), ktraj_adc(3,samps), 'r.', 'MarkerSize', 12);
axis equal;
title(sprintf('2D sampling mask. R = %d', round(R)), 'FontSize', 18);
xlabel('k_y (m^{-1})', 'FontSize', 18); ylabel('k_z (m^{-1})', 'FontSize', 18);
xlim([-Ny*deltak(2)/2, Ny*deltak(2)/2]); ylim([-Nz*deltak(3)/2, Nz*deltak(3)/2]);
return;

%% Plot point spread function in y-z space
omega = omegas(:,:,1);
psf = ifftshift(ifft2(fftshift(omega)));
[Y, Z] = meshgrid(res(2)*(-Ny/2:Ny/2-1), res(3)*(-Nz/2:Nz/2-1));
figure; surf(Y.', Z.', abs(psf), 'FaceColor', 'interp');
xlim([-Ny*res(2)/2, Ny*res(2)/2]); ylim([-Nz*res(3)/2, Nz*res(3)/2]);
zlim([0, max(abs(psf(:)))]);
xlabel('y (m)', 'FontSize', 18); ylabel('z (m)', 'FontSize', 18);
zlabel('magnitude (a.u.)', 'FontSize', 18);
title('Point spread function in y-z space', 'FontSize', 18);
zticks(linspace(0, max(abs(psf(:))), 5)); zticklabels(linspace(0, 1, 5));
return;

%% Detailed sequence report (slow)
rep = seq.testReport;
fprintf([rep{[1:9, 11:end]}]); % skip warnings
return;
