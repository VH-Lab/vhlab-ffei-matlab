% DELAY_FC: Generates figures comparing output power at the input
% modulation frequency of FFE vs. FFEI triad synapse circuit models with
% different inhibitory delays (Figure 4).s

% Input characteristics
PR = 100;
F = logspace(log10(5), log10(1000), 50);
dt = 1e-4; % time bins
tmax = 5;
num_inputs = 10; % number of trials

% Conductance parameters
Pmax_e_scaled = 1.6976e-7;
Pmax_e = [1.6976e-7*0.305 1.6976e-7*0.335 1.6976e-7*0.395 1.6976e-7*0.45 1.6976e-7*0.48 1.6976e-7*0.47]; % for delay set 1, no noise

tau1e = 0.020; % 20ms
tau2e = 0.001; % 1ms
tau1i = zeros(1, length(Pmax_e));
tau1i(1:end-1) = 0.02;
tau2i = tau2e; % 1ms
delay = [0.001 0.002 0.005 0.010 0.020 0.001];

tvec = 0:dt:tmax;

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
tau_m = 1.0e-2; % time
V_syn_e = 0; % synaptic resting potential
V_syn_i = -0.08; % synaptic resting potential
rmgs = 0.25; % original: 0.05
Nt = length(tvec); % number of time intervals
Im = zeros(Nt, 1); % input current


% BALANCING THE E/I CURRENTS in FFEI case

alpha = 1.25;



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
            
            % 
            output = run_triad_model(1, 'F', F(freq), 'tmax', tmax, 'Pmax_e', Pmax_e(tau), 'tau1i', tau1i(tau), 'delay', delay(tau));
            FC_n(n) = output.FC;
            FC_a(n) = output.FC_avg;
            FC_p(n) = output.FC_pct;
            
        end

        FC(tau, freq) = mean(FC_n);
        FC_avg(tau,freq) = mean(FC_a);
        FC_pct(tau,freq) = mean(FC_p);
        
    end
end


cmap = parula(6);
figure;
hold on;

for i = 1:5
    plot(F,FC(i,:), 'color', cmap(i,:));
end
plot(F, FC(6, :), 'r-')

title(["Avg. power at the input modulation frequency of the postsynaptic", ...
    "firing rate, given different synaptic inhibitory fall time constants,", ...
    ['synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
xlabel("Input modulation frequency (Hz)")
ylabel(["FC magnitude at the given input", ...
    "modulation frequency (Hz) (FC_F)"])
legend("delay = " + num2str(delay(1)), "delay = " + num2str(delay(2)), "delay = " + num2str(delay(3)), "delay = " + num2str(delay(4)), "delay = " + num2str(delay(5)), "excitatory")
set(gca, 'XScale', 'log')
box off;


figure;
hold on;

for i = 1:5
    plot(F,FC_avg(i,:), 'color', cmap(i,:));
end
plot(F, FC_avg(6, :), 'r-')

title(["Avg. overall power of the postsynaptic firing rate for inputs with", ...
    "different modulation frequencies, given different synaptic inhibitory", ...
    ['fall time constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])

xlabel("Input modulation frequency (Hz)")
ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
%ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
legend("delay = " + num2str(delay(1)), "delay = " + num2str(delay(2)), "delay = " + num2str(delay(3)), "delay = " + num2str(delay(4)), "delay = " + num2str(delay(5)), "excitatory")
set(gca, 'XScale', 'log')
box off;



figure;
hold on;

for i = 1:5
    plot(F,FC_pct(i,:), 'color', cmap(i,:));
end
plot(F, FC_pct(6, :), 'r-')

plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])

title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
    "postsynaptic firing rate, given different synaptic inhibitory fall time", ...
    ['constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])

xlabel("Input modulation frequency (Hz)")
ylabel("FC_F/FC_a_v_g")
%ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
legend("delay = " + num2str(delay(1)), "delay = " + num2str(delay(2)), "delay = " + num2str(delay(3)), "delay = " + num2str(delay(4)), "delay = " + num2str(delay(5)), "excitatory")
set(gca, 'XScale', 'log')
box off;
