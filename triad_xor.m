% TRIAD_XOR: Generates the panels of Figure 6 and Figure 6 supplement.
%   Comparison of neuromorphic retina-LGN 1-to-1 circuits that perform an
%   XOR computation, exhibiting 1) one-sign transmission or 2) paired
%   two-sign tranmission across layers.

clear;

% Generate spike traces for Figure 6 (panel A)
disp('Generating voltage trace panels (Figure 6)...')
xor_spiketrace(5,pi)
clear;

% Generate FC figures for Figure 6 (panels (B-D))
disp('Generating comparisons of FFE vs FFEI (Figure 6)...')
xor_fc(pi)
clear;


% Generate spike traces for Figure 6supp (panel A)
disp('Generating voltage trace panels (Figure 6supp)...')
xor_spiketrace(5,pi/2)
clear;

% Generate FC figures for Figure 6supp (panels (B-D))
disp('Generating comparisons of FFE vs FFEI (Figure 6supp)...')
xor_fc(pi/2)
clear;