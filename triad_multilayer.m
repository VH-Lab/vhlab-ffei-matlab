% TRIAD_MULTILAYER: Generates the panels of Figure 5.
%   Comparison of 4-layer retina-LGN 1-to-1 circuits, exhibiting 1) noise 
%   only, 2)feedforward-E transmission only or 3) feedforward-E/I
%   transmission (triad synapse).

clear;
% Generate spike traces (panel A)
disp('Generating spike traces...')
multilayer_spiketrace();

clear;
% Generate FC plots (panel B-E)
disp('Generating FC figures...')
multilayer_fc();



