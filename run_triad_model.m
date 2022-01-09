function [output] = run_triad_model(modelinput, is_spike_present, varargin)
% RUN_TRIAD_MODEL: Performs simulations of retina-LGN 1-to-1 circuit in
% response to an input specified by generate_input.m in a modelinput
% struct. Returns FC at input frequency, FC over all frequencies,
% normalized FC, Poisson spike train output, voltage trace, E/I synaptic
% conductance, input spike train, summed noisy input spikes, synaptic
% current, and leak current (contained in a struct (output)).
% Inputs can vary depending on which parameters are modified in the model
% (see how assign.m in vhlab toolbox is used when calling functions).
% When examining the behavior of the model, the parameter is_spike_present
% determines whether the desired model is spiking (1, Vm resets at the
% threshold) or non-spiking (0, Vm does not reset).

% Conductance parameters
%Pmax_e = 1.6976e-7;
Pmax_e = 1.6976e-7*0.305;
tau1e = 0.020; % 20ms
tau2e = 0.001; % 1ms
tau1i = 0.020; % 20ms
tau2i = 0.001; % 1ms
delay = 0.001; % 2ms

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
Cm = 1.0e-9; % membrane capacitance
V_syn_e = 0; % synaptic resting potential
V_syn_i = -0.08; % synaptic resting potential


% BALANCING THE E/I CURRENTS in FFEI case
alpha = 1.25;

% define noise parameters (if any)
noise_base = 1.6976e-7/75;
noise_level = 0; % scaling factor for base noise strength

% Assign values to any parameters specified in function call:
assign(varargin{:});
% typically: PR, F, tmax, dt, Pmax_e, tau1i, delay, alpha, Im_F, Im_amp...

tau_m = Rm*Cm; % membrane time constant; depends on potential assignment of Rm, Cm in function call

% Set input parameters according to modelinput struct (includes time
% parameters for the trial)
F = modelinput.F;
dt = modelinput.dt;
tvec = modelinput.tvec;
input_spktrain = modelinput.signal;
noisy_input = modelinput.noise;
Im = modelinput.Im;
Nt = length(tvec); % number of time intervals

% Scale the noisy input magnitude by original 1.6976e-7/75 (same as in figure 4)
[Ps_E_noise] = ppsc_constantsum(noisy_input, noise_level*noise_base, tau1e, tau2e, 0, tau2i, delay, dt);

% Calculate synaptic conductance due to spike train using
% ppsc_constantsum.m
[Ps_E, Ps_I] = ppsc_constantsum(input_spktrain, Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);

% Simulate voltage behavior of a neuron recieving push-pull conductance
% input
Vm = zeros(Nt, 1); % store membrane voltage at each time point
Vm(1) = V_reset; % initial membrane voltage is V_reset

tspk = []; % count spikes from Euler's method
post_spktrain = zeros(Nt, 1);

% Track current behavior
I_leak = zeros(Nt, 1);
I_syn = zeros(Nt, 1);

% Track behavior when spiking is removed from the model
Vm_noreset = Vm;
I_leak_noreset = I_leak;
I_syn_noreset = I_syn;


for t = 1:Nt-1 % Euler method
    
    I_syn_noreset_e = (Ps_E(t) + Ps_E_noise(t)) * (V_syn_e - Vm_noreset(t));
    I_syn_noreset_i = Ps_I(t) * (V_syn_i - Vm_noreset(t)) * alpha; % multiply inhibitory current by alpha
    I_syn_noreset(t) = I_syn_noreset_e + I_syn_noreset_i; 
    I_leak_noreset(t) = -(Vm_noreset(t) - V_e)/Rm;
    dvdt_noreset = (Rm * (I_leak_noreset(t) + I_syn_noreset(t) + Im(t))) / tau_m;
    Vm_noreset(t+1) = Vm_noreset(t) + dt * dvdt_noreset;
    
    I_leak(t) = -(Vm(t) - V_e)/Rm;
    I_syn(t) = (Ps_E(t) + Ps_E_noise(t)) * (V_syn_e - Vm(t)) + Ps_I(t) * (V_syn_i - Vm(t)) * alpha;
    
    if Vm(t) >= V_th
        Vm(t+1) = V_reset;
        tspk = [tspk; t*dt];
        post_spktrain(t) = 1;
    else
        %dvdt = (-(Vm(t) - V_e) - (Rm * ((Ps_E(t) + Ps_E_noise(t)) * (Vm(t) - V_syn_e) + alpha * Ps_I(t) * (Vm(t) - V_syn_i))) + Rm*Im(t)) / tau_m;
        dvdt = (Rm * (I_leak(t) + I_syn(t) + Im(t))) / tau_m;
        Vm(t+1) = Vm(t) + dt * dvdt;
    end
end

if is_spike_present == 0 % If the desired model is non-spiking, output non-spiking model behavior
    Vm = Vm_noreset;
    I_syn = I_syn_noreset;
    I_leak = I_leak_noreset;
elseif is_spike_present ~= 1 % If desired model is spiking, behavior of spiking model is outputted by default; this error statement is just here to ensure that is_leak_present was inputted as a boolean.
    error("Error: Invalid input for the boolean is_leak_present (must input 0 or 1).")
end

% Calculate fourier coefficients as the firing rate following a particular
% frequency
postratevec = post_spktrain/dt;
FC_all = fouriercoeffs(postratevec, dt);

FC = abs(fouriercoeffs_tf2(postratevec, F, 1/dt));
FC_avg = mean(abs(FC_all));

if FC_avg == 0
    FC_pct = 0; % Set normalized FC to zero if average overall power = 0 (i.e. no output spikes occur)
else
    FC_pct = FC/FC_avg; % Otherwise, perform normalization by dividing output power at F (FC) by average overall power (FC_avg)
end

output = var2struct('FC', 'FC_avg', 'FC_pct', 'post_spktrain', 'Vm', 'Ps_E', 'Ps_I', 'input_spktrain', 'noisy_input', 'I_syn', 'I_leak');

end