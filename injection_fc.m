% INJECTION_FC: Generates figures comparing spiking behavior of
% model neurons in response to rectified sinusoidal current injections and
% FFE or FFEI synaptic input.

% Input characteristics
PR = 100;
F = logspace(log10(5), log10(1000), 50);
dt = 1e-4; % time bins
tmax = 5;
num_inputs = 10;

% Conductance parameters
Pmax_e_scaled = 1.6976e-7;
Pmax_e = [1.6976e-7*0.305 1.6976e-7*0.47 0]; % scaled Pmax so that power at 5Hz = 75Hz, alpha = 1.25
tau1i = [0.02 0 0]; % 50ms % 0.02 0.025 0.03 0.05 0.25 1 0
Im_amp = [0 0 8.38e-9];

FC = zeros(length(tau1i), length(F)); % store FCs
FC_avg = zeros(length(tau1i), length(F)); % store avg power for E/I transmission
FC_pct = zeros(length(tau1i), length(F)); % store % power at each freq for E/I transmission

disp('Starting simulations...')

for freq = 1:length(F)
    
    if mod(freq-1,5) == 0
        disp(['Progress: ' num2str(freq-1) ' frequencies analyzed'])
    end
    
    for tau = 1:length(tau1i)
        
        FC_n = zeros(1, num_inputs); % store FC of LIF neuron for each spike train
        FC_a = zeros(1, num_inputs); % store avg power of LIF neuron for each spike train
        FC_p = zeros(1, num_inputs); % store percent power for each spike train
        
        for n = 1:num_inputs
            
            % simulate output power of triad synapse (or cell given current
            % injection) using run_triad_model.m
            [FC1,FC2,FC3] = run_triad_model(1, 'F', F(freq), 'tmax', tmax, 'Pmax_e', Pmax_e(tau), 'tau1i', tau1i(tau), 'Im_amp', Im_amp(tau));
            
            FC_n(n) = FC1;
            FC_a(n) = FC2;
            FC_p(n) = FC3;
            
        end
        
        FC(tau, freq) = mean(FC_n);
        FC_avg(tau,freq) = mean(FC_a);
        FC_pct(tau,freq) = mean(FC_p);
        
    end
end

% Plot FC figures
cmap = parula(5);
figure;
hold on;
plot(F, FC(1, :), 'color', cmap(1,:))
plot(F, FC(2, :), 'r-')
plot(F, FC(3, :), 'm-')

xlabel("Input modulation frequency (Hz)")
ylabel(["FC magnitude at the given input", ...
    "modulation frequency (Hz) (FC_F)"])
legend("Balanced FFEI", "FFE", "8.38nA current injection")
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
plot(F, FC_avg(1, :), 'color', cmap(1,:))
plot(F, FC_avg(2, :), 'r-')
plot(F, FC_avg(3, :), 'm-')

xlabel("Input modulation frequency (Hz)")
ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
legend("Balanced FFEI", "FFE", "8.38nA current injection")
set(gca, 'XScale', 'log')
box off;

figure;
hold on;
plot(F, FC_pct(1, :), 'color', cmap(1,:))
plot(F, FC_pct(2, :), 'r-')
plot(F, FC_pct(3, :), 'm-')

xlabel("Input modulation frequency (Hz)")
ylabel("FC_F/FC_a_v_g")
legend("Balanced FFEI", "FFE", "8.38nA current injection")
set(gca, 'XScale', 'log')
box off;