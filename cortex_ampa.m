% CORTEX_AMPA: Generates the panels of Figure 10.
%   Comparison of cortical 100-to-1 circuits, exhibiting 1) FFE projections
%   to cortical cell only, 2) FFE projections to inhibitory interneuron
%   and cortical cell + FFI projection from interneuron to cortical cell,
%   or 3) Same as 2, but with AMPA only inputs to inhibitory interneuron

clear;

% Generate voltage traces (panel B)
F = [5 20 90];
disp('Generating voltage trace panels...')
for i = 1:length(F)
    cortex_voltagetrace(F(i), 1);
end

clear;

% Generate FC figures (panels (C-F))
disp('Generating comparisons of FFE vs FFEI over taum of I cell...')
cortex_fc;
clear;