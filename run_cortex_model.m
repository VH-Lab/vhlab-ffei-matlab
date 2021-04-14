function [output] = run_cortex_model(F, tau_m_I, varargin)
% RUN_CORTEX_MODEL: Simulates the response to multiple paired-EI inputs.
% All spiking inputs are generated from the same sinusoidal signal,
% which has a certain amplitude - denoting peak firing rate of signal (PR)
% - and modulation frequency (F), lasting t_end seconds.
%   - An excitatory cell recieves n_inputs via excitatory transmission
%   (the sum of the inputs is scaled by Ps_E_scale).
%   - An inhibitory cell recieves same n_inputs via excitatory transmission
%   (the sum of the inputs is scaled by Ps_I_scale).
% The excitatory cell recieves the response of the inhibitory cell via
% inhibitory transmission.
% Note: associated with Figure 8.

% Time and input ST parameters
PR = 100;
phase_shift = 0;
dt = 1e-4;
t_start = 0;
t_end = 1;

% Conductance parameters
Pmax_e = 1e-7;
Pmax_i = 1e-7;
tau1e = 0.020; % 20ms % taufall
tau2e = 0.001; % 1ms % taurise
tau1i = 0.020; % 20ms
tau2i = 0.001; % 1ms
alpha = 1.25; % inhibitory scaling coefficient

% Define input, input scaling parameters
n_inputs = 100;
Ps_E_scale = 0.02; % scaling factor for summed input->E
Ps_I_scale = 0.01; % scaling factor for summed input->I

% define LIF parameters - E cell
V_reset_E = -0.080; % -80mV
V_e_E = -0.075; % -75mV
V_th_E = -0.040; % -40mV
Rm_E = 1.0e7; % membrane resistance
tau_m_E = 1.0e-2; % membrane time constant
V_syn_E_e = 0; % synaptic resting potential
V_syn_E_i = -0.08; % synaptic resting potential


% define LIF parameters - I cell
V_reset_I = -0.080;
V_e_I = -0.075;
Vth_factor = 1.2;
Rm_I = 1.0e7;
V_syn_I = 0;

% Assign values to any parameters specified in function call:
assign(varargin{:});
% typically: PR, F, tmax, dt, Pmax_e, tau1i, delay, alpha, Im_F, Im_amp...

% Define any parameters based on unique values included in function call:
tvec = t_start:dt:t_end;
Nt = length(tvec); % number of time intervals
V_th_I = -0.040 * Vth_factor; % threshold voltage for I cell; 20-50% lower than for E cell


% Prepare inputs
inputs = zeros(n_inputs, Nt);
for i = 1:n_inputs
    s1 = spiketrain_sinusoidal(PR, F, phase_shift, 0, t_start, t_end, dt);
    inputs(i,:) = spiketimes2bins(s1, tvec);
end


% Prepare conductance of each input - E
Ps_E = zeros(n_inputs, Nt);
for i = 1:n_inputs
    % when inhibitory parameters are set to 0, ppsc_constantsum computes
    % synaptic conductance as if only one form of transmission is present
    % (no additional balancing for paired transmission)
    Ps_E(i,:) = ppsc_constantsum(inputs(i,:), Pmax_e, tau1e, tau2e, 0, 0, 0, dt);
end
Ps_E_total = sum(Ps_E,1) * Ps_E_scale;

% Conductance of each input, scaled - I
Ps_I_total = sum(Ps_E,1) * 0.6 * Ps_I_scale;


% Initialize vectors to store voltage and spikes over time for both cells
Vm_E = zeros(1,Nt);
Vm_E(1) = V_reset_E;
output_E = zeros(1,Nt);

Vm_I = zeros(1,Nt);
Vm_I(1) = V_reset_I;
output_I = zeros(1,Nt);


% Use Euler method to determine output of I neuron over time
for t = 1:Nt-1 % Euler method
    if Vm_I(t) >= V_th_I
        Vm_I(t+1) = V_reset_I;
        output_I(t) = 1;
    else
        I_leak = -(Vm_I(t) - V_e_I)/Rm_I;
        I_syn = Ps_I_total(t) * (V_syn_I - Vm_I(t));
        dvIdt = Rm_I * (I_leak + I_syn) / tau_m_I;
        %dvIdt = (-(Vm_I(t) - V_e_I) - (Rm_I * (Ps_I_total(t)) * (Vm_I(t) - V_syn_I))) / tau_m_I; 
        Vm_I(t+1) = Vm_I(t) + dt * dvIdt;
    end
end


% Determine conductance of output_I to E neuron
Ps_I_output = ppsc_constantsum(output_I, Pmax_i, tau1i, tau2i, 0, 0, 0, dt);

% Use Euler method to determine output of E neuron over time
for t = 1:Nt-1 % Euler method
    if Vm_E(t) >= V_th_E
        Vm_E(t+1) = V_reset_E;
        output_E(t) = 1;
    else
        I_leak = -(Vm_E(t) - V_e_E)/Rm_E;
        I_syn_e = Ps_E_total(t) * (V_syn_E_e - Vm_E(t));
        I_syn_i = alpha * Ps_I_output(t) * (V_syn_E_i - Vm_E(t));
        I_syn = I_syn_e + I_syn_i;
        dvEdt = Rm_E * (I_leak + I_syn) / tau_m_E;
        %dvEdt = (-(Vm_E(t) - V_e_E) - (Rm_E * ((Ps_E_total(t) * (Vm_E(t) - V_syn_E_e)) +  alpha * Ps_I_output(t) * (Vm_E(t) - V_syn_E_i)))) / tau_m_E;
        Vm_E(t+1) = Vm_E(t) + dt * dvEdt;
    end
end

output_rate = output_E/dt; % convert output of E cell to rate - excitatory output

% Calculate power of output
FC = abs(fouriercoeffs_tf2(output_rate(:), F, 1/dt));

FC_all = fouriercoeffs(output_rate(:), dt);
FC_avg = mean(abs(FC_all));
if FC_avg <= 0
    FC_pct = 0;
else
    FC_pct = FC/FC_avg;
end


output = var2struct('FC', 'FC_avg', 'FC_pct', 'Ps_I_total', 'output_I', 'Ps_I_output', 'Ps_E_total', 'output_E', 'tvec');

end