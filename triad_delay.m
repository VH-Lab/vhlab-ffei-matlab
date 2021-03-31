% TRIAD_DELAY: Generates the panels of Figure 4.
%   Examination of inhibitory delay in retina-LGN 1-to-1 FFEI circuits.

clear;

% Generate voltage traces (panel A)
disp('Generating voltage trace panels...')
delay = [0.001 0.005 0.020];
Pmax_e = [1.6976e-7*0.305 1.6976e-7*0.395 1.6976e-7*0.48];
tau1i = [0.02 0.02 0.02];
for i = 1:length(delay)
   delay_voltagetrace(delay(i), Pmax_e(i), tau1i(i));
end

clear;
disp('Generating comparisons of FFE vs FFEI over inhibitory delay...')
delay_fc;

