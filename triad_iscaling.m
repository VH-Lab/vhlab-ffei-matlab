% TRIAD_ISCALING: Generates the panels of Figure 3.
%   Examination of inhibitory scaling + noise in retina-LGN 1-to-1
%   circuits.

clear;

% Generate voltage traces (panel A)
disp('Generating voltage traces...')
% Parameters for each of the shown voltage traces:
alpha = [0 1 1.25 2 5];
Pmax_e = [1.6976e-7*0.47 1.6976e-7*0.15 1.6976e-7*0.42 1.6976e-7*1.4 1.6976e-7*0.95];
tau1i = [0 0.02 0.02 0.02 0.02];

for i = 1:length(alpha)
   alpha_voltagetrace(alpha(i), Pmax_e(i), tau1i(i)); 
end

% Generate FC figures (panels (B-D))
disp('Generating comparisons of FFE vs FFEI over alpha...')
alpha_fc(0, -0.075);

% Note: Panels E-I involve changing the parameters in line 20 and manually
% measuring FC_F/FC_avg at F = 50Hz and F = 100Hz.


