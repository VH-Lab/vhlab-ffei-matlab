function [] = ang_depression_plot(plot_type)
% ANG_DEPRESSION_PLOT: Performs simulations of retina-LGN 1-to-1 circuit in
% response to Poisson spike train inputs with rectified sinusoidal rate,
% given transmission via 1) AMPARs + NMDARs with depression, or 2) AMPARs +
% NMDARs + GABARs with depression.
%
% Inputs specify the following:
%   plot_type:          0: plot voltage traces for input with F = 20Hz
%                       1: plot FC at F, overall FC, and normalized FC at F
%                          for a range of frequencies F
%
% Interspike times in seconds and s2/s1 ratio of synaptic currents are
% measured from literature (Chen et al. 2002, Blitz and Regehr 2005)

if plot_type == 0 % voltage trace plot
    F = 20;
    tmax = 2;
elseif plot_type == 1 % fc plot
    F = logspace(log10(5), log10(1000), 50);
    tmax = 5;
else
    error("Error: Invalid plot type (0 or 1). Please see help text for plot type options.")
end

% Input characteristics
PR = 50;
phase_shift = 0;
dt = 1e-4; % time bins
tvec = 0:dt:tmax;
Nt = length(tvec);
num_inputs = 10;

% Conductance parameters - AMPA parameters for excitatory, GABA parameters
% for inhibitory
Pmax_e_base = 3e-6;
Pmax_e_a = 2 * Pmax_e_base;
Pmax_e_n = 0.2 * Pmax_e_base;
Pmax_factor = [2.3 4.7];

figure_title = ["AMPA+NMDA",'AMPA+NMDA, GABA_A + GABA_B'];
tau1e_a = 0.0016;
tau2e_a = 0.0007; % AMPA conductance parameters
tau1e_n = 0.090;
tau2e_n = 0.0032; % NMDA conductance parameters

is_inhibition_present = [0 1];
tau1i_a = 0.006;
tau2i_a = 0.0006; % GABA_A conductance parameters
tau1i_b = 0.150;
tau2i_b = 0.03; % GABA_B conductance parameters (approximate) -> NOTE that GABA_B actually has an additional slow fall time constant...
delay = 0.001; % 1ms delay
alpha = 1; % inhibitory scaling coefficient

% Accounting for GABA_B currents...
Pmax_i = 1.5*generate_balanced_EI([Pmax_e_a Pmax_e_n], [tau1e_a tau1e_n], [tau2e_a tau2e_n], [3 2], [tau1i_a tau1i_b], [tau2i_a tau2i_b]);
Pmax_i_a = Pmax_i(1)*0.8;
Pmax_i_b = Pmax_i(2)*0.45;

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
tau_m = 1.0e-2; % time
V_syn_e = 0; % synaptic resting potential (excitatory)
V_syn_i_a = -0.08; % synaptic resting potential (inhibitory, GABA_A)
V_syn_i_b = -0.1; % synaptic resting potential (inhibitory,  GABA_B)


% Channel adaptation data
load('depression_model_AMPA');
load('depression_model_NMDA');
load('depression_model_GABA');
 
Dampa = depression_model_AMPA;
Dnmda = depression_model_NMDA;
Dgaba = depression_model_GABA;

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
            s1 = spiketrain_sinusoidal(PR, F(freq), phase_shift, 0, 0, tmax, dt); % Generate Poisson spike times
            spktrain = spiketimes2bins(s1, tvec); % Convert to spike train over time bins
            
            % apply depression model to input spike train s1, and compute postsynaptic
            % conductances
            syncurrs_ampa = depression_model_comp(s1, Dampa.a0, Dampa.f, Dampa.ftau, Dampa.d, Dampa.dtau);
            spktrain_ampa = spktrain;
            spktrain_ampa(find(spktrain_ampa)) = syncurrs_ampa;
            P_ampa = synaptic_conductance(spktrain_ampa, Pmax_e_a, tau1e_a, tau2e_a, dt, 0);
            
            % the above is repeated for the other channel types
            
            % NMDA
            
            syncurrs_nmda = depression_model_comp(s1, Dnmda.a0, Dnmda.f, Dnmda.ftau, Dnmda.d, Dnmda.dtau);
            spktrain_nmda = spktrain;
            spktrain_nmda(find(spktrain_nmda)) = syncurrs_nmda;
            P_nmda = synaptic_conductance(spktrain_nmda, Pmax_e_n, tau1e_n, tau2e_n, dt, 0);
            
            % GABA_A, B
            
            syncurrs_gaba = depression_model_comp(s1, Dgaba.a0, Dgaba.f, Dgaba.ftau, Dgaba.d, Dgaba.dtau);
            spktrain_gaba = spktrain;
            spktrain_gaba(find(spktrain_gaba)) = syncurrs_gaba;
            P_gabaA = synaptic_conductance(spktrain_gaba, Pmax_i_a, tau1i_a, tau2i_a, dt, delay);
            P_gabaB = synaptic_conductance(spktrain_gaba, Pmax_i_b, tau1i_b, tau2i_b, dt, delay);
            
            Ps_E = (P_ampa+P_nmda)*Pmax_factor(tau);
            Ps_I = P_gabaA+P_gabaB*Pmax_factor(tau);

            
            % Simulate voltage behavior of a neuron recieving input via
            % adaptive channels
            
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
                    I_syn_e =  Pmax_factor(tau) * (P_ampa(t) + P_nmda(t)) * (V_syn_e - Vm(t));
                    I_syn_i = is_inhibition_present(tau) * Pmax_factor(tau) * alpha * (P_gabaA(t) * (V_syn_i_a - Vm(t)) + P_gabaB(t) * (V_syn_i_b - Vm(t)));
                    I_syn = I_syn_e + I_syn_i;
                    dvdt = Rm * (I_leak + I_syn) / tau_m;
                    Vm(t+1) = Vm(t) + dt * dvdt;
                    
                end
            end
            
            postratevec = post_spktrain(ceil(Nt/tmax):end)/dt;
            
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
            plot(tvec, spktrain*PR/2, '-', 'Color', [0.8 0.8 0])
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
            ylabel('Synaptic conductance (S)')
            xlabel('Time (s)')
            legend('Excitatory (AMPA + NMDA)', 'Inhibitory (GABA)')
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
            ylabel('Output membrane potential (V)')
            
        end
        
    end
end


if plot_type == 1
    
    figure;
    plot(F, FC(1, :), 'r-')
    hold on;
    plot(F, FC(2, :), 'b-')
    
    title(["Avg. power at the input modulation frequency of the postsynaptic", ...
    "firing rate, given adapting feed-forward AMPA-NMDA/GABA transmission"])
    
    xlabel("Input modulation frequency (Hz)")
    ylabel(["FC magnitude at the given input", ...
        "modulation frequency (Hz) (FC_F)"])
    
    legend('AMPA+NMDA','AMPA+NMDA, GABA_A + GABA_B')
    set(gca, 'XScale', 'log')
    box off;
    
    
    figure;
    plot(F, FC_avg(1, :), 'r-')
    hold on;
    plot(F, FC_avg(2, :), 'b-')
    
    title(["Avg. overall power of the postsynaptic firing rate for inputs with", ...
    "different modulation frequencies, given adapting feed-forward AMPA-NMDA/GABA transmission"])
    xlabel("Input modulation frequency (Hz)")
    ylabel("Avg. FC magnitude over entire spectrum (Hz) (FC_a_v_g)")
    ylimits = ylim;
    ylim([0 ylimits(2)]);
    legend('AMPA+NMDA','AMPA+NMDA, GABA_A + GABA_B')
    set(gca, 'XScale', 'log')
    box off;
    
    
    
    figure;
    plot(F, FC_pct(1, :), 'r-')
    hold on;
    plot(F, FC_pct(2, :), 'b-')
    plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])
    
    title(["(Avg. power at input modulation frequency)/(Avg. power overall) of the", ...
    "postsynaptic firing rate, given adapting feed-forward AMPA-NMDA/GABA transmission"])
 
    xlabel("Input modulation frequency (Hz)")
    ylabel("FC_F/FC_a_v_g")
    ylimits = ylim;
    ylim([0 ylimits(2)]);
    legend('AMPA+NMDA','AMPA+NMDA, GABA_A + GABA_B')
    set(gca, 'XScale', 'log')
    box off;
end

end
