function [] = ang_plot(plot_type, is_noise_present)
% RUN_TRIAD_MODEL: Performs simulations of retina-LGN 1-to-1 circuit in
% response to Poisson spike train inputs with rectified sinusoidal rate,
% given transmission via 1) AMPARs, 2) AMPARs + NMDARs, or 3) AMPARs +
% NMDARs + GABARs.
% Inputs specify the following:
%   plot_type:          0: plot voltage traces for input with F = 50Hz
%                       1: plot FC at F, overall FC, and normalized FC at F
%                          for a range of frequencies F
%   is_noise_present:   0: no background activity
%                       1: low level excitatory background activity (scaled
%                          as in Figure 5)

if plot_type == 0 % voltage trace plot
    F = 50;
elseif plot_type == 1 % fc plot
    F = logspace(log10(5), log10(1000), 50);
else
    error("Error: Invalid plot type (0 or 1). Please see help text for plot type options.")
end

if is_noise_present == 0
    noise_strength = 0; % no noise
elseif is_noise_present == 1
    noise_strength = 1; % low level excitatory noise
else
    error("Error: Invalid value for is_noise_present (boolean).")
end

% Noisy input parameters
tau1e = 0.02; % 20ms % taufall
tau2e = 0.001; % 1ms % taurise

% Input characteristics
PR = 100;
dt = 1e-4; % time bins
tmax = 0.5;
num_inputs = 10;

% Conductance parameters - AMPA parameters for excitatory, GABA parameters
% for inhibitory
Pmax_e_a = 1.6976e-7;
Pmax_e_n = 0.5 * Pmax_e_a;
Pmax_factor = [4 0.3 0.09 0];
figure_title = ["AMPA", 'AMPA-NMDA', 'AMPA-NMDA/GABA', ['noise only x ' num2str(noise_strength)]];
tau1e_a = 0.00072;
tau2e_a = 0.0007; % AMPA conductance parameters
tau1e_n = 0.100;
tau2e_n = 0.0032; % NMDA conductance parameters
is_inhibition_present = [0 0 1 0];
ampa_only = [1 0 0 0];
delay = 0.001; % 1ms delay
alpha = 1.25; % inhibitory scaling coefficient

tvec = 0:dt:tmax;
tvec_buffer = 0:dt:tmax + delay;

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
tau_m = 1.0e-2; % time
V_syn_e = 0; % synaptic resting potential (excitatory)
V_syn_i = -0.08; % synaptic resting potential (inhibitory)
rmgs = 0.25; % original: 0.05
Nt = length(tvec); % number of time intervals


FC = zeros(length(is_inhibition_present), length(F)); % store FCs
FC_avg = zeros(length(is_inhibition_present), length(F)); % store avg power for E/I transmission
FC_pct = zeros(length(is_inhibition_present), length(F)); % store % power at each freq for E/I transmission

disp('Starting simulations...')

for freq = 1:length(F)
    
    if mod(freq-1,5) == 0
        disp(['Progress: ' num2str(freq-1) ' frequencies analyzed'])
    end
    
    for tau = 1:length(is_inhibition_present)
             
        FC_n = zeros(1, num_inputs); % store FC of LIF neuron for each spike train
        FC_a = zeros(1, num_inputs); % store avg power of LIF neuron for each spike train
        FC_p = zeros(1, num_inputs); % store percent power for each spike train
        
        for n = 1:num_inputs
            
            % Generate input + noise spike trains
            input = generate_input(F(freq), 1, 1, 'tmax', tmax);
            spktrain = input.signal;
            noisy_input = input.noise;
            
            [Ps_E_noise] = ppsc_constantsum(noisy_input, noise_strength*1.6976e-7/75, tau1e, tau2e, 0, 0, 0, dt);
            
            % push-pull conductance triggered by Poisson spike train
            if ampa_only(tau) == 1 % case 1: ampa only conductance
                [Ps_E, Ps_I] = ppsc_constantsum(spktrain, Pmax_e_a*Pmax_factor(tau), tau1e_a, tau2e_a, 0, 0, delay, dt);
            elseif Pmax_factor(tau) == 0 % case 4: noise only
                Ps_E = zeros(1, length(tvec));
                Ps_I = zeros(1, length(tvec));
            else % case 2 and 3: ampa+nmda (w or w/o gaba) conductance
                [Ps_E, Ps_I] = ppsc_constantsum_multcond(spktrain, [Pmax_e_a Pmax_e_n]*Pmax_factor(tau), [tau1e_a tau1e_n], [tau2e_a tau2e_n], is_inhibition_present(tau), delay, dt);
            end
            
            % Simulate voltage behavior of a neuron recieving push-pull conductance
            % input
            
            Vm = zeros(Nt, 1); % store membrane voltage at each time point
            Vm(1) = V_reset; % initial membrane voltage is V_reset
            
            % track spikes from Euler's method
            tspk = []; % count spikes
            post_spktrain = zeros(Nt, 1);
            
            
            for t = 1:Nt-1 % Euler method
                
                if Vm(t) >= V_th
                    Vm(t+1) = V_reset;
                    tspk = [tspk; t*dt];
                    post_spktrain(t) = 1;
                else
                    I_leak = -(Vm(t) - V_e)/Rm;
                    I_syn = (Ps_E(t) + Ps_E_noise(t)) * (V_syn_e - Vm(t)) + alpha * Ps_I(t) * (V_syn_i - Vm(t));
                    dvdt = Rm * (I_leak + I_syn) / tau_m;
                    %dvdt = (-(Vm(t) - V_e) - (Rm * (Ps_E(t) * (Vm(t) - V_syn_e) + alpha * Ps_I(t) * (Vm(t) - V_syn_i)))) / tau_m;
                    Vm(t+1) = Vm(t) + dt * dvdt;
                end
            end
            
            postratevec = post_spktrain/dt;
            
            FC_all = fouriercoeffs(postratevec, dt);
            
            FC_n(n) = abs(fouriercoeffs_tf2(postratevec, F(freq), 1/dt));
            FC_a(n) = mean(abs(FC_all));
            if FC_a(n) == 0
                FC_p(n) = 0;
            else
                FC_p(n) = FC_n(n)/FC_a(n);
            end
            
        end
        
        FC(tau, freq) = mean(FC_n);
        FC_avg(tau,freq) = mean(FC_a);
        FC_pct(tau,freq) = mean(FC_p);
        
        if plot_type == 0 % plot sample traces for F
            signal = max(0,PR*sin(2*pi*F(freq)*tvec));
            
            figure;
            subplot(3,1,1)
            plot(tvec, spktrain*50, '-', 'Color', [0.8 0.8 0])
            hold on;
            plot(tvec, signal)
            ylabel('Firing rate (Hz)')
            xlabel('Time (s)')
            title(figure_title(tau))
            box off;
            
            subplot(3,1,2)
            plot(tvec, Ps_E)
            hold on;
            plot(tvec, Ps_I)
            plot(tvec, Ps_E_noise)
            ylabel('Conductance (S)')
            xlabel('Time (s)')
            legend('Ps_E', 'Ps_I', 'Ps_E_,_n_o_i_s_e')
            box off;
            
            subplot(3,1,3)
            plot(tvec, signal * (0.1/PR) - 0.1, ':', 'Color', [1 0.733 0.318])
            hold on;
            plot(tvec, Vm)
            ylim([-0.1 0.05])
            box off;
            
            if isempty(tspk) == 0
                plot([tspk(:) tspk(:)], [-0.08 0], 'Color', [0.396 0.569 1], 'HandleVisibility', 'off')
                plot(tspk, 0, 'x', 'Color', 'red')
            end
        end
        
    end
end


if plot_type == 1
    
    figure;
    plot(F, FC(1, :), 'k-')
    hold on;
    plot(F, FC(2, :), 'r-')
    plot(F, FC(3, :), 'b-')
    plot(F, FC(4, :), 'g-')
    
    title(["Avg. power at the input modulation frequency of the postsynaptic", ...
        "firing rate, given AMPA transmission alone, AMPA-NMDA transmission alone,", ...
        "or feed-forward AMPA-NMDA/GABA transmission"])
    
    xlabel("Input modulation frequency (Hz)")
    ylabel(["FC magnitude at the given input", ...
        "modulation frequency (Hz) (FC_F)"])
    
    legend("AMPA", "AMPA-NMDA", "AMPA-NMDA/GABA", ['noise only x ' num2str(noise_strength)])
    set(gca, 'XScale', 'log')
    box off;
    
    
    figure;
    plot(F, FC_avg(1, :), 'k-')
    hold on;
    plot(F, FC_avg(2, :), 'r-')
    plot(F, FC_avg(3, :), 'b-')
    plot(F, FC_avg(4, :), 'g-')
    
    title(["Avg. overall power of the postsynaptic firing rate for inputs with", ...
        "different modulation frequencies, given AMPA transmission alone, AMPA-NMDA", ...
        "transmission alone, or feed-forward AMPA-NMDA/GABA transmission"])
    xlabel("Input modulation frequency (Hz)")
    ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
    ylimits = ylim;
    ylim([0 ylimits(2)]);
    legend("AMPA", "AMPA-NMDA", "AMPA-NMDA/GABA", ['noise only x ' num2str(noise_strength)])
    set(gca, 'XScale', 'log')
    box off;
    
    
    
    figure;
    plot(F, FC_pct(1, :), 'k-')
    hold on;
    plot(F, FC_pct(2, :), 'r-')
    plot(F, FC_pct(3, :), 'b-')
    plot(F, FC_pct(4, :), 'g-')
    plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])
    
    title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
        "postsynaptic firing rate, given AMPA transmission alone, AMPA-NMDA", ...
        "transmission alone, or feed-forward AMPA-NMDA/GABA transmission"])
    xlabel("Input modulation frequency (Hz)")
    ylabel("FC_F/FC_a_v_g")
    ylimits = ylim;
    ylim([0 ylimits(2)]);
    legend("AMPA", "AMPA-NMDA", "AMPA-NMDA/GABA", ['noise only x ' num2str(noise_strength)])
    set(gca, 'XScale', 'log')
    box off;
end

end
