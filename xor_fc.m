function [] = xor_fc(phase_shift)
% XOR_FC: Plots the output power (FC) over input modulation frequency of a
% neuromorphic XOR model circuit with 1) FFE or 2) FFEI projections between
% units, recieving phase-shifted inputs (Figure 6).

dt = 1e-4/5;
Freq = logspace(log10(5), log10(1000), 50);
tmax = 1;

% Conductance parameters
is_inhibition_present = [0 1];
Pmax_base = 1.6976e-7;
if phase_shift == pi
    Pmax_e = Pmax_base * [0.39 0.33]; % 2xF parameters
    FC_freq_factor = 2;
elseif phase_shift == pi/2
    Pmax_e = Pmax_base * [0.4 0.3]; % 1xF parameters
    FC_freq_factor = 1;
else
    error('Error: xor_voltagetrace currently only balances transmission for phase shifts = pi/2 or pi.')
end

% store fourier coefficient values
FC = zeros(length(is_inhibition_present), length(Freq));
FC_avg = zeros(length(is_inhibition_present), length(Freq));
FC_pct = zeros(length(is_inhibition_present), length(Freq));

n_trials = 10; % run test n_trials times

for i = 1:length(is_inhibition_present)
    disp(['ANALYZING case - is inhibition present? = ' num2str(is_inhibition_present(i)) '...'])
    for F = 1:length(Freq)
        if  F == 1 || mod(F,20) == 0
        disp(['Analyzing frequencies = ' num2str(Freq(F)) 'Hz and above...']);
        end
        
        % Set frequency at which FCs are computed
        FC_freq = FC_freq_factor*Freq(F);
        
        trial_FC = zeros(1,n_trials);
        trial_FC_avg = zeros(1,n_trials);
        trial_FC_pct = zeros(1,n_trials);
        
        for n_trial = 1:n_trials

            % Generate phase-shifted inputs
            input1 = generate_input(Freq(F), 1, 0, 'phase_shift', 0, 'dt', dt, 'tmax', tmax);
            input2 = generate_input(Freq(F), 1, 0, 'phase_shift', phase_shift, 'dt', dt, 'tmax', tmax);

            output = run_xor_model(input1, input2, Pmax_e(i), is_inhibition_present(i), FC_freq);
            
            trial_FC(n_trial) = output.FC;
            trial_FC_avg(n_trial) = output.FC_avg;
            trial_FC_pct(n_trial) = output.FC_pct;
            
        end
        
        % take average fourier coefficient for 10 trials
        
        FC(i,F) = mean(trial_FC);
        FC_avg(i,F) = mean(trial_FC_avg);
        FC_pct(i,F) = mean(trial_FC_pct);
        
    end
end

figure;
hold on;
plot(Freq, FC(2,:), 'Color', [0.3 0.5 1])
plot(Freq, FC(1,:), 'Color', [1 0.3 0.3])

xlabel('Input modulation frequency (Hz)')
ylabel('Avg. FC at the frequency (Hz) (FC_F)')
title({'Avg. power of the XOR circuit output at different input modulation' ...
    ['frequencies: phase shift = ' num2str(phase_shift) ' rad, peak rate = 100Hz,'] ...
    'excitatory vs. paired EI transmission'})
legend(['FFEI, ' num2str(FC_freq_factor) 'xF'], ['E only, ' num2str(FC_freq_factor) 'xF'])
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
plot(Freq, FC_avg(2,:), 'Color', [0.3 0.5 1])
plot(Freq, FC_avg(1,:), 'Color', [1 0.3 0.3])

xlabel('Input modulation frequency (Hz)')
ylabel('Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)')
title({'Avg. overall power of the XOR circuit output for inputs with different' ...
    ['modulation frequencies: phase shift = ' num2str(phase_shift) ' rad,'] ...
    'peak rate = 100Hz, excitatory vs. paired EI transmission'})
legend('FFEI', 'E only')
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
plot(Freq, FC_pct(2,:), 'Color', [0.3 0.5 1])
plot(Freq, FC_pct(1,:), 'Color', [1 0.3 0.3])
plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])

xlabel('Input modulation frequency (Hz)')
ylabel('FC_F/FC_a_v_g')
title({'(Avg. power at input modulation frequency)/(Avg. power overall) of then' ...
    ['XOR circuit output: phase shift = ' num2str(phase_shift) ' rad,'] ...
    'peak rate = 100Hz, excitatory vs. paired EI transmission'})
legend(['FFEI, ' num2str(FC_freq_factor) 'xF'], ['E only, ' num2str(FC_freq_factor) 'xF'])
set(gca, 'XScale', 'log')
box off;

end