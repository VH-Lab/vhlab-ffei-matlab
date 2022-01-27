% TRIAD DEPRESSION FIT: fits experimental data on AMPAR, NMDAR, GABAR
% adaptation to a depression model given by depression_model_fit.m (2
% depressing factors, 0 facilitating factors).

% Note that depression models have been pre-generated using this script and
% are saved in this directory under "depression_model_...". Running this
% script as-is will generate identical depression models to the saved
% versions, and will plot their behavior relative to the experimental data.

% Channel adaptation data

% ampa_interval = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 4 8];
% ampa_ratio = [0.238 0.279 0.378 0.418 0.503 0.619 0.722 0.820 0.892 0.934];
% 
% nmda_interval = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 4 8];
% nmda_ratio = [0.510 0.510 0.559 0.607 0.675 0.765 0.775 0.854 0.900 0.962];
% 
% % since depression is mostly presynaptic in the case of GABA channels, we
% % will assume that GABA_A and GABA_B follow the same depression model
% % (measured for GABA_A)
% gaba_interval = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 4];
% gaba_ratio = [0.103 0.131 0.194 0.303 0.340 0.588 0.673 0.838 0.984];


load('AMPA_depression_data')
load('NMDA_depression_data')
load('GABA_depression_data')

% construct spike trains with pairs of spikes separated by 20s
% e.g time 0: 1, time 0.01: 0.510

% SEE generate_depression_exp.m FOR FUNCTION THAT RECONSTRUCTS PAIRED PULSE
% EXPERIMENTS FROM DATA

[ampa_spiketimes, ampa_current] = generate_depression_exp(ampa_interval, ampa_ratio);
[nmda_spiketimes, nmda_current] = generate_depression_exp(nmda_interval, nmda_ratio);
[gaba_spiketimes, gaba_current] = generate_depression_exp(gaba_interval, gaba_ratio);

% Fit depression model to data for each channel type

% fit AMPA data to depression model using depression_model_fit
[a0, f, ftau, d, dtau, err] = depression_model_fit(ampa_spiketimes, ampa_current, 0, 2);

% Export depression model
depression_model_AMPA = var2struct('a0', 'f', 'ftau', 'd', 'dtau');

% check if depression model is a good fit for the data (comparison figure at bottom)
%   apply depression model to spike timings used for raw data, and compute
%   resulting s2/s1 ratio
ampa_fit = depression_model_comp(ampa_spiketimes, a0, f, ftau, d, dtau);
ampa_fit = ampa_fit(2:2:end) ./ ampa_fit(1:2:end);

% the above is repeated for other channel types

% NMDA

[a0, f, ftau, d, dtau, err] = depression_model_fit(nmda_spiketimes, nmda_current, 0, 2);

depression_model_NMDA = var2struct('a0', 'f', 'ftau', 'd', 'dtau');

nmda_fit = depression_model_comp(nmda_spiketimes, a0, f, ftau, d, dtau);
nmda_fit = nmda_fit(2:2:end) ./ nmda_fit(1:2:end);

% GABA

[a0, f, ftau, d, dtau, err] = depression_model_fit(gaba_spiketimes, gaba_current, 0, 2);

depression_model_GABA = var2struct('a0', 'f', 'ftau', 'd', 'dtau');

gaba_fit = depression_model_comp(gaba_spiketimes, a0, f, ftau, d, dtau);
gaba_fit = gaba_fit(2:2:end) ./ gaba_fit(1:2:end);



% save('depression_model_AMPA')
% save('depression_model_NMDA')
% save('depression_model_GABA')

% depression model fit comparisons
figure;

subplot(1,3,1)
plot(ampa_interval, ampa_ratio, 'ko:')
hold on;
plot(ampa_interval, ampa_fit, 'ro-')
title('AMPA depression model')
xlabel('Interspike interval (s)')
ylabel('S2/S1')
legend('Raw data', 'Fit', 'Location', 'southeast')
ylim([0 1])
box off;

subplot(1,3,2)
plot(nmda_interval, nmda_ratio, 'ko:')
hold on;
plot(nmda_interval, nmda_fit, 'ro-')
title('NMDA depression model')
xlabel('Interspike interval (s)')
ylabel('S2/S1')
legend('Raw data', 'Fit', 'Location', 'southeast')
ylim([0 1])
box off;

subplot(1,3,3)
plot(gaba_interval, gaba_ratio, 'ko:')
hold on;
plot(gaba_interval, gaba_fit, 'ro-')
title('GABA depression model')
xlabel('Interspike interval (s)')
ylabel('S2/S1')
legend('Raw data', 'Fit', 'Location', 'southeast')
ylim([0 1])
box off;