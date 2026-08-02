%% main.m  Generate all 4 sequences for an ArbEPI scan session.
%
% Edit params.m to configure the experiment, then run this script.
% Sequences are written to the output/ directory.
%
% Order matters: ArbEPI must run first because EPIcal and noise load
% samp_locs.mat that it produces.

clear; close all;
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
addpath(fullfile(projRoot, 'src'));
params;

%% 1. Generate sampling masks and main EPI sequence
omegas = gen_sampling_masks(R);
plot_kt_mask(omegas);
ArbEPI(omegas);

%% 2. Calibration sequence (ghost correction + receiver gain)
EPIcal();

%% 3. Gradient echo reference (sensitivity maps)
GRE();

%% 4. Noise prescan (noise covariance)
noise();