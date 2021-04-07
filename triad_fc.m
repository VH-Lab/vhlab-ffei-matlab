function [] = triad_fc(comparison)
% TRIAD_FC: Generates figures comparing output power at the input
% modulation frequency of FFE vs. FFEI triad synapse circuit models in
% response to Poisson inputs with rectified sinusoidal rate (Figure 1).
    % If comparison = 1, then compare models with different
    % inhibitory fall time constants tau1i.
    % If comparison = 2, then compare FFE and FFEI over range of
    % output powers.
    

% Input characteristics
PR = 100;
F = logspace(log10(5), log10(1000), 50);
dt = 1e-4; % time bins
tmax = 5;
num_inputs = 10;

% Conductance parameters
%Pmax_e = 1.6976e-7;
Pmax_e_scaled = 1.6976e-7;
Pmax_e = [1.6976e-7*0.305 1.6976e-7*0.52 1.6976e-7*0.56 1.6976e-7*0.44 1.6976e-7*0.47]; % scaled Pmax so that power at 5Hz = 75Hz, alpha = 1.25
Pmax_factor_ffei = [0.41 0.75 1 1.2 1.31]; % FFEI, for alpha = 1.25
Pmax_factor_ffe = [0.4 0.675 1 1.5 2]; % FFE
tau1i = [0.02 0.025 0.03 0.05 0]; % 50ms % 0.02 0.025 0.03 0.05 0.25 1 0

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
            
            % generate input
            input = generate_input(F(freq), 1, 0, 'tmax', tmax);
            
            % use when comparing different tau1i
            if comparison == 1
                output = run_triad_model(input, 1, 'Pmax_e', Pmax_e(tau), 'tau1i', tau1i(tau));
            elseif comparison == 2
            % use when comparing balanced E/I over output power
                output = run_triad_model(input, 1, 'Pmax_e', Pmax_e(1)*Pmax_factor_ffei(tau), 'tau1i', tau1i(1));
            else
                error('Error: Comparision type not correctly defined. Set input to 1 if comparing models over tau1i, 2 if comparing over output power.')
            end
            
            FC_n(n) = output.FC;
            FC_a(n) = output.FC_avg;
            FC_p(n) = output.FC_pct;
            
        end
        
        FC(tau, freq) = mean(FC_n);
        FC_avg(tau,freq) = mean(FC_a);
        FC_pct(tau,freq) = mean(FC_p);
        
    end
end


% Generate figures for comparison over tau1i:
if comparison == 1
    cmap = parula(5);
    figure;
    hold on;
    for i = 1:4
        plot(F,FC(i,:), 'color', cmap(i,:));
    end
    plot(F, FC(5, :), 'r-')
    
    title(["Avg. power at the input modulation frequency of the postsynaptic", ...
        "firing rate, given different synaptic inhibitory fall time constants,", ...
        ['synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
    xlabel("Input modulation frequency (Hz)")
    ylabel(["FC magnitude at the given input", ...
        "modulation frequency (Hz) (FC_F)"])
    legend("0.02s", "0.025s", "0.03s", "0.05s", "excitatory")
    set(gca, 'XScale', 'log')
    box off;
    
    
    figure;
    hold on;
    for i = 1:4
        plot(F,FC_avg(i,:), 'color', cmap(i,:));
    end
    plot(F, FC_avg(5, :), 'r-')
    
    title(["Avg. overall power of the postsynaptic firing rate for inputs with", ...
        "different modulation frequencies, given different synaptic inhibitory", ...
        ['fall time constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
    
    xlabel("Input modulation frequency (Hz)")
    ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
    %ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
    legend("0.02s", "0.025s", "0.03s", "0.05s", "excitatory")
    %legend("P_m_a_x*1/4", "P_m_a_x*1/2", "P_m_a_x", "P_m_a_x*2", "P_m_a_x*4")
    set(gca, 'XScale', 'log')
    box off;
    
    
    
    figure;
    hold on;
    for i = 1:4
        plot(F,FC_pct(i,:), 'color', cmap(i,:));
    end
    plot(F, FC_pct(5, :), 'r-')
    plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])
    
    title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
        "postsynaptic firing rate, given different synaptic inhibitory fall time", ...
        ['constants, synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
    
    xlabel("Input modulation frequency (Hz)")
    ylabel("FC_F/FC_a_v_g")
    %ylabel("Avg. postsynaptic firing rate as spikes/duration of trial")
    legend("0.02s", "0.025s", "0.03s", "0.05s", "excitatory")
    %legend("P_m_a_x*1/4", "P_m_a_x*1/2", "P_m_a_x", "P_m_a_x*2", "P_m_a_x*4")
    set(gca, 'XScale', 'log')
    box off;
    
    
elseif comparison == 2
    
    
    cmap = parula(5);
    figure;
    hold on;
    for i = 5:-1:2
        plot(F,FC(i,:), 'color', cmap(6-i,:));
    end
    plot(F, FC(1, :), 'r-')
    
    title(["Avg. power at the input modulation frequency of the postsynaptic", ...
        "firing rate, given different output powers,", ...
        ['synaptic conductance scaled from ' num2str(Pmax_e_scaled) ' S']])
    xlabel("Input modulation frequency (Hz)")
    ylabel(["FC magnitude at the given input", ...
        "modulation frequency (Hz) (FC_F)"])
    set(gca, 'XScale', 'log')
    box off;
    
    % run comparison over output power for FFE case
    FC = zeros(length(tau1i), length(F)); % store FCs
    FC_avg = zeros(length(tau1i), length(F)); % store avg power for E/I transmission
    FC_pct = zeros(length(tau1i), length(F)); % store % power at each freq for E/I transmission
    
    disp('Starting simulations...')
    
    for freq = 1:length(F)
        
        if mod(freq-1,5) == 0
            disp(['Progress: ' num2str(freq-1) ' frequencies analyzed'])
        end
        
        for tau = 1:length(tau1i)
            
            % Get Poisson spike trains using input characteristics
            
            FC_n = zeros(1, num_inputs); % store FC of LIF neuron for each spike train
            FC_a = zeros(1, num_inputs); % store avg power of LIF neuron for each spike train
            FC_p = zeros(1, num_inputs); % store percent power for each spike train
            
            for n = 1:num_inputs
                
                input = generate_input(F(freq), 1, 0, 'tmax', tmax);
                output = run_triad_model(input, 1, 'Pmax_e', Pmax_e(5)*Pmax_factor_ffe(tau), 'tau1i', tau1i(5));
                
                FC_n(n) = output.FC;
                FC_a(n) = output.FC_avg;
                FC_p(n) = output.FC_pct;
            
            end
            
            FC(tau, freq) = mean(FC_n);
            FC_avg(tau,freq) = mean(FC_a);
            FC_pct(tau,freq) = mean(FC_p);
            
        end
    end
    
    for i = 5:-1:2
        plot(F,FC(i,:), '--', 'color', cmap(6-i,:));
    end
    plot(F, FC(1, :), 'r--')
    
end

end 