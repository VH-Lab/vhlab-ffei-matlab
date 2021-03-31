function [] = alpha_fc(noise_level)
%ALPHA_FC: Generates figures comparing output power at the input
% modulation frequency of triad synapse circuit models with different
% inhibitory scaling (Figure 3).


% Input characteristics
PR = 100;
F = logspace(log10(5), log10(1000), 50);
dt = 1e-4/5; % time bins
tmax = 2;
num_inputs = 10;

% Conductance parameters
Pmax_e_scaled = 1.6976e-7;

switch noise_level
    case 0
        Pmax_e = [1.6976e-7*0.15 1.6976e-7*0.42 1.6976e-7*0.8 1.6976e-7*1.2 1.6976e-7*1.4 1.6976e-7*1.12 1.6976e-7*0.95 1.6976e-7*0.72 1.6976e-7*0.53 1.6976e-7*0.47 0]; % for alpha set 1, no noise (last 2 values for E only and noise only)
    case 0.5
        Pmax_e = [1.6976e-7*0.135 1.6976e-7*0.41 1.6976e-7*0.8 1.6976e-7*1.21 1.6976e-7*1.3885 1.6976e-7*1.24 1.6976e-7*0.935 1.6976e-7*0.7 1.6976e-7*0.56 1.6976e-7*0.36 0]; % for alpha set 1, noise level x0.5
    case 1
        Pmax_e = [1.6976e-7*0.12 1.6976e-7*0.4 1.6976e-7*0.8 1.6976e-7*1.2 1.6976e-7*1.37 1.6976e-7*1.26 1.6976e-7*0.94 1.6976e-7*0.73 1.6976e-7*0.53 1.6976e-7*0.26 0]; % for alpha set 1, noise level x1
    case 1.5
        Pmax_e = [1.6976e-7*0.065 1.6976e-7*0.355 1.6976e-7*0.5 1.6976e-7*0.8 1.6976e-7*0.96 1.6976e-7*1.17 1.6976e-7*0.88 1.6976e-7*0.7 1.6976e-7*0.53 1.6976e-7*0.28 0]; % for alpha set 1, noise level x1.5
    case 2
        Pmax_e = [1.6976e-7*0.04 1.6976e-7*0.33 1.6976e-7*0.24 1.6976e-7*0.4 1.6976e-7*0.6 1.6976e-7*1.12 1.6976e-7*0.9 1.6976e-7*0.67 1.6976e-7*0.54 1.6976e-7*0.3 0]; % for alpha set 1, noise level x2
    case 5
        Pmax_e = [1.6976e-7*0.03 1.6976e-7*0.21 1.6976e-7*0.073 1.6976e-7*0.045 1.6976e-7*0.032 1.6976e-7*0.013 1.6976e-7*0.006 1.6976e-7*0.0035 1.6976e-7*0.0025 1.6976e-7*0.3 0]; % for alpha set 1, noise level x5
    otherwise
        error('Error: No set of synaptic strengths to match the given noise level. Valid noise levels are 0, 0.5, 1, 1.5, 2, 5.')
end

tau1i = zeros(1, length(Pmax_e));
tau1i(1:end-2)=0.02;

% BALANCING THE E/I CURRENTS in FFEI case

alpha = [1 1.25 1.5 1.75 2 3 5 7 10 1 1];



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
            
            [FC1, FC2, FC3] = run_triad_model('dt', dt, 'F', F(freq), 'tmax', tmax, 'Pmax_e', Pmax_e(tau), 'tau1i', tau1i(tau), 'noise_level', noise_level, 'alpha', alpha(tau));
            
            FC_n(n) = FC1;
            FC_a(n) = FC2;
            FC_p(n) = FC3;
            
        end

        FC(tau, freq) = mean(FC_n);
        FC_avg(tau,freq) = mean(FC_a);
        FC_pct(tau,freq) = mean(FC_p);
        
    end
end


cmap = parula(7);
figure;
hold on;
plot(F, FC(1,:), 'color', cmap(1,:));
for i = 2:6
    plot(F,FC(i+3,:), 'color', cmap(i,:));
end
plot(F, FC(10, :), 'r-')
plot(F, FC(11, :), 'Color', [0.3 0.3 0.3])

title(["Avg. power at the input modulation frequency of the postsynaptic", ...
    "firing rate, given different synaptic inhibitory fall time constants,", ...
    ['synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
xlabel("Input modulation frequency (Hz)")
ylabel(["FC magnitude at the given input", ...
    "modulation frequency (Hz) (FC_F)"])
legend("alpha = " + num2str(alpha(1)), "alpha = " + num2str(alpha(5)), "alpha = " + num2str(alpha(6)), "alpha = " + num2str(alpha(7)), "alpha = " + num2str(alpha(8)), "alpha = " + num2str(alpha(9)), "excitatory", "noise x" + num2str(noise_level))
set(gca, 'XScale', 'log')
box off;


figure;
hold on;
plot(F, FC_avg(1,:), 'color', cmap(1,:));
for i = 2:6
    plot(F,FC_avg(i+3,:), 'color', cmap(i,:));
end
plot(F, FC_avg(10, :), 'r-')
plot(F, FC_avg(11, :), 'Color', [0.3 0.3 0.3])

title(["Avg. overall power of the postsynaptic firing rate for inputs with", ...
    "different modulation frequencies, given different synaptic inhibitory", ...
    ['fall time constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])

xlabel("Input modulation frequency (Hz)")
ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
%ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
legend("alpha = " + num2str(alpha(1)), "alpha = " + num2str(alpha(5)), "alpha = " + num2str(alpha(6)), "alpha = " + num2str(alpha(7)), "alpha = " + num2str(alpha(8)), "alpha = " + num2str(alpha(9)), "excitatory", "noise x" + num2str(noise_level))
set(gca, 'XScale', 'log')
box off;



figure;
hold on;
plot(F, FC_pct(1,:), 'color', cmap(1,:));
for i = 2:6
    plot(F,FC_pct(i+3,:), 'color', cmap(i,:));
end
plot(F, FC_pct(10, :), 'r-')
plot(F, FC_pct(11, :), 'Color', [0.3 0.3 0.3])

plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])

title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
    "postsynaptic firing rate, given different synaptic inhibitory fall time", ...
    ['constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])

xlabel("Input modulation frequency (Hz)")
ylabel("FC_F/FC_a_v_g")
%ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
legend("alpha = " + num2str(alpha(1)), "alpha = " + num2str(alpha(5)), "alpha = " + num2str(alpha(6)), "alpha = " + num2str(alpha(7)), "alpha = " + num2str(alpha(8)), "alpha = " + num2str(alpha(9)), "excitatory", "noise x" + num2str(noise_level))
set(gca, 'XScale', 'log')
box off;


% additional figure showing normalized FC for alpha 1-2
cmap2 = parula(6);
figure;
hold on;
for i = 1:5
    plot(F,FC_pct(i,:), 'color', cmap2(i,:));
end
plot(F, FC_pct(10, :), 'r-')
plot(F, FC_pct(11, :), 'Color', [0.3 0.3 0.3])

plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])

title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
    "postsynaptic firing rate, given different synaptic inhibitory fall time", ...
    ['constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])

xlabel("Input modulation frequency (Hz)")
ylabel("FC_F/FC_a_v_g")
%ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
legend("alpha = " + num2str(alpha(1)), "alpha = " + num2str(alpha(2)), "alpha = " + num2str(alpha(3)), "alpha = " + num2str(alpha(4)), "alpha = " + num2str(alpha(5)), "excitatory", "noise x" + num2str(noise_level))
set(gca, 'XScale', 'log')
box off;

end