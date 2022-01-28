% TRIAD_VS_INJECTION: Generates the panels of Figure 2.
%   Comparison of retina-LGN 1-to-1 circuits to direct current injections.

% Generate voltage traces due to current injection (panel A)
clear;
disp('Generating voltage trace panels...')
injection_voltagetrace;

% Generate FC figures (panels B-D)
clear;
disp('Generating comparisons of FFE vs FFEI vs current injection...')
injection_fc;

% Generate FC figures for E only cell (panels E-G)
disp('Generating comparisons of FFE over tau_m...')
triad_fc(3);