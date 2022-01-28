% TRIAD_DEPRESSION: Generates the panels of Figure 8.
%   Comparison of retina-LGN 1-to-1 circuits, exhibiting postsynaptic 1)
%   AMPA+NMDA channels with synaptic depression, or 2)
%   AMPA+NMDA+GABAA+GABAB channels with synaptic depression.

clear;
% Generate voltage traces (panels A-B)
disp('Generating voltage traces...')
ang_depression_plot(0);

clear;
% Generate FC plots (panels C-E)
disp('Generating FC figures...')
ang_depression_plot(1);

clear;
% Fit channel depression data to a model of synaptic depression (Supp. Fig.
% 3)
disp("Generating comparisons between channel depression data and fitted depression model...")
triad_depression_fit;
