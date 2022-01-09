% CORTEX_FC: Generates figures comparing output power at the input
% modulation frequency of FFE vs. FFEI cortex circuit models in
% response to Poisson inputs with rectified sinusoidal rate (Figure 9).
% (e.g. LGN->cortex transmission)
% In FFEI case, the membrane time constant of the inhibitory interneuron is
% varied.

% Time and input ST parameters
F_vec = logspace(log10(5), log10(1000), 50);
n_trials = 10;
Pmax_e = 1e-7; % conductance of E output produced by multiple inputs
Pmax_i = zeros(1,7); % conductance of I output by inhibitory interneuron
Pmax_i(1:6) = Pmax_e;

% tau_m of I cell is varied in the FFEI circuit
tau_m_I = [0.005 0.006 0.01 0.02 0.05 0.1 0.01]; % (default)

FC_vec = zeros(length(Pmax_i), length(F_vec)); % 1: PFFEI, 2: E only
FCavg_vec = zeros(length(Pmax_i), length(F_vec));
FCpct_vec = zeros(length(Pmax_i), length(F_vec));

% Note on order of inputs:
% run_cortex_model takes 1) frequency, 2) tau_m of I cell, 3) any additional
% parameters


for i = 1:length(F_vec)
    
    if mod(i-1,5) == 0
    disp(['Progress: ' num2str(i-1) ' frequencies analyzed'])
    end
    
    % store specific output power (FC), overall output power (FCavg), and
    % normalized output power (FCpct) from each model (FFE, FFEI with
    % varied tau_m_I) in arrays
    FC = zeros(length(Pmax_i), n_trials);
    FCavg = zeros(length(Pmax_i), n_trials);
    FCpct = zeros(length(Pmax_i), n_trials);
    
    for j = 1:n_trials
        for tau = 1:length(tau_m_I)
            output = run_cortex_model(F_vec(i), tau_m_I(tau), 'Pmax_e', Pmax_e, 'Pmax_i', Pmax_i(tau));
            FC(tau,j) = output.FC;
            FCavg(tau,j) = output.FC_avg;
            FCpct(tau,j) = output.FC_pct;
        end
    end
    
    % for each model and over all trials, average the output power for a
    % signal of frequency denoted by Fvec(i)
    for model = 1:length(Pmax_i)
        FC_vec(model,i) = mean(FC(model,:));
        FCavg_vec(model,i) = mean(FCavg(model,:));
        FCpct_vec(model,i) = mean(FCpct(model,:));
    end
end

cmap = parula(length(Pmax_i));
figure;
hold on;
for i = 1:length(Pmax_i)-1
    plot(F_vec,FC_vec(i,:), 'color', cmap(i,:));
end
plot(F_vec, FC_vec(length(Pmax_i),:), 'r')

xlabel('Input modulation frequency (Hz)')
ylabel('FC magnitude of output at the frequency F (Hz) (FC_F)')
legend('0.005s', '0.006s', '0.01s', '0.02s', '0.05s', '0.1s', 'E only') % (C)
title('Avg. output power of FFEI multi-input circuit at different input frequencies')
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
for i = 1:length(Pmax_i)-1
    plot(F_vec,FCavg_vec(i,:), 'color', cmap(i,:));
end
plot(F_vec, FCavg_vec(length(Pmax_i),:), 'r')

xlabel('Input modulation frequency (Hz)')
ylabel({'Avg. FC magnitude over entire spectrum  of  output (Hz) (FC_a_v_g)'})
legend('0.005s', '0.006s', '0.01s', '0.02s', '0.05s', '0.1s', 'E only') % (C)
title('Avg. overall output power of FFEI multi-input circuit')
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
for i = 1:length(Pmax_i)-1
    plot(F_vec,FCpct_vec(i,:), 'color', cmap(i,:));
end
plot(F_vec, FCpct_vec(length(Pmax_i),:), 'r')

plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])
xlabel('Input modulation frequency (Hz)')
ylabel('FC_F/FC_a_v_g')
legend('0.005s', '0.006s', '0.01s', '0.02s', '0.05s', '0.1s', 'E only') % (C)
title('Normalized output power of FFEI multi-input circuit')
set(gca, 'XScale', 'log')
box off;



% For unscaled Pmax plot, manually scale all FC_pct plots to start at 1 Hz
% power
FCpct_scaled_vec = zeros(length(Pmax_i), length(F_vec));
for i = 1:length(Pmax_i)
    FCpct_scaled_vec(i,:) = FCpct_vec(i,:)/FCpct_vec(i,1);
end

figure;
hold on;
for i = 1:length(Pmax_i)-1
    plot(F_vec,FCpct_scaled_vec(i,:), 'color', cmap(i,:));
end
plot(F_vec, FCpct_scaled_vec(length(Pmax_i),:), 'r')

plot([5 1e3],[1/2 1/2],'--', 'Color', [0.5 0.5 0.5])
xlabel('Input modulation frequency (Hz)')
ylabel('FC_F/FC_a_v_g')
legend('0.005s', '0.006s', '0.01s', '0.02s', '0.05s', '0.1s', 'E only') % (C)
title('Normalized output power of FFEI multi-input circuit, scaled to 1 for F=5Hz')
set(gca, 'XScale', 'log')
box off;
