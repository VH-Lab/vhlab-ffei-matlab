% TRIAD_AMPAGABANMDA: Generates the panels of Figure 7.
%   Comparison of retina-LGN 1-to-1 circuits, exhibiting postsynaptic 1)
%   AMPA channels only, 2) AMPA+NMDA channels only (similar to FFE case),
%   or 3) AMPA+NMDA+GABA channels (similar to FFEI case). 

clear;
% Generate voltage traces (panel A)
disp('Generating voltage traces...')
ang_plot(0,0);

clear;
% Generate FC plots without E background noise (panels B-D)
disp('Generating FC figures without excitatory background noise...')
ang_plot(1,0);

clear;
% Generate FC plots with E background noise (panels E-G)
disp('Generating FC figures with excitatory background noise...')
ang_plot(1,1);