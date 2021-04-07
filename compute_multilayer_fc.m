function [FC,FC_avg,FC_pct] = compute_multilayer_fc(PR, F, L, Pmax_factor, PFFEI, transmission_exists)
% COMPUTE_MULTILAYER_FC: Computes output power (FC, FCavg, FCpct) for a single multilayer
% simulation.

% Set up time parameters
dt = 1e-4;
phase_shift = 0;
t_start = 0;
t_end = 1;
tvec = t_start:dt:t_end;



% Conductance parameters
Pmax_e = 1.6976e-7*Pmax_factor;
tau1e = 0.02; % 20ms % taufall
tau2e = 0.001; % 1ms % taurise
if PFFEI
    tau1i = 0.02; % 20ms % set to 0 to simulate only excitatory transmission
else
    tau1i = 0;
end
tau2i = 0.001; % 1ms
delay = 0.001; % 1ms

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
alpha = 1.25;

% Using the previously-determined solution for the firing rate of a neuron
% given a constant input current, we can solve for the background input
% current required to make the neuron fire at a rate of PR/pi
T_background = pi/PR;


inputs = zeros(L+1,Nt);
Vm = zeros(L,Nt); % store membrane potential from each layer
Vm(:,1) = V_reset; % initial membrane voltage is V_reset

s1 = spiketrain_sinusoidal(PR, F, phase_shift, 0, tvec(1), tvec(end), dt);
inputs(1,:) = spiketimes2bins(s1, tvec); % comment out for testing purposes


for i = 1:L
    
    [Ps_E, Ps_I] = ppsc_constantsum(inputs(i,:), Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);
    
    if transmission_exists ~= 1
        Ps_E(:) = 0;
        Ps_I(:) = 0;
    end
    
    % 12/10: An alternative method: provide random noise in the form of
    % spike trains --> excitatory conductance to each input
    noisy_input = poisson_spiketrain(dt, 1/T_background, t_end, 1);
    for j=2:50
        noisy_input = noisy_input + poisson_spiketrain(dt, 1/T_background, t_end, 1);
    end
    
    % Scale the noisy input magnitude by original Pmax_e/75
    [Ps_E_noise] = ppsc_constantsum(noisy_input, 1.6976e-7/75, tau1e, tau2e, 0, tau2i, delay, dt);
    %Ps_E_noise(:) = 0;
    
    for t = 1:Nt-1 % Euler method
        
        if Vm(i,t) >= V_th
            Vm(i,t+1) = V_reset;
            inputs(i+1,t) = 1;
        else
            dvdt = (-(Vm(i,t) - V_e) - (Rm * ((Ps_E(t) + Ps_E_noise(t)) * (Vm(i,t) - V_syn_e) + alpha * Ps_I(t) * (Vm(i,t) - V_syn_i))) + Rm*Im(t)) / tau_m;
            Vm(i,t+1) = Vm(i,t) + dt * dvdt;
        end
    end
    
end

output = inputs(L+1, :)/dt; % convert output to rate
FC = abs(fouriercoeffs_tf2(output(:), F, 1/dt));

FC_all = fouriercoeffs(output(:), dt);
FC_avg = mean(abs(FC_all));
FC_pct = FC/FC_avg;

end

