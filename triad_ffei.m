% TRIAD_FFEI: Generates the panels of Figure 1.
%   Comparison of retina-LGN 1-to-1 circuits, exhibiting 1) feedforward-E
%   transmission only or 2) feedforward-E/I transmission (triad synapse).

clear;

% Generate voltage traces (panels A-B)
disp('Generating voltage trace panels...')
triad_voltagetrace(1);
clear;
triad_voltagetrace(2);
clear;

% Generate FC figures (panels (C-F))
disp('Generating comparisons of FFE vs FFEI over tau1i...')
triad_fc(1);
clear;
disp('Generating comparisons of FFE vs FFEI over output power...')
triad_fc(2);
clear;

% Generate voltage traces (Supp. Fig. 1B (panel A is identical to Fig. 1A)
triad_voltagetrace(3);
clear;



