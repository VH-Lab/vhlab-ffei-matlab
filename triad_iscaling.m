% TRIAD_ISCALING: Generates the panels of Figure 3.
%   Examination of inhibitory scaling + noise in retina-LGN 1-to-1
%   circuits.

clear;

% Generate voltage traces (panel A)
disp('Voltage traces under construction...')
% Parameters for each of the shown voltage traces:
alpha = [0 1 1.25 2 5];
Pmax_e = [1.6976e-7*0.47 1.6976e-7*0.15 1.6976e-7*0.42 1.6976e-7*1.4 1.6976e-7*0.95];
tau1i = [0 0.02 0.02 0.02 0.02];

% Generate FC figures (panels (B-D))
disp('Generating comparisons of FFE vs FFEI over alpha...')
alpha_fc(1);