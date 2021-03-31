function [FC, FC_avg, FC_pct, post_spktrain, Vm, Ps_E, Ps_I, spktrain] = run_triad_model(varargin)
% RUN_TRIAD_MODEL: Performs simulations of retina-LGN 1-to-1 circuit in
% response to a Poisson input with a rectified sinusoidal rate. Returns
% Poisson spike train output, voltage trace, FC at input frequency, FC over
% all frequencies, and normalized FC.


% Input characteristics
PR = 100;
F = 50;
dt = 1e-4; % time bins
tmax = 5; % signal time
tbuffer = 0; % "buffer" time after signal

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
tau_m = 1.0e-2; % time
V_syn_e = 0; % synaptic resting potential
V_syn_i = -0.08; % synaptic resting potential


% BALANCING THE E/I CURRENTS in FFEI case
alpha = 1.25;

% define noise parameters (if any)
noise_base = 1.6976e-7/75;
noise_level = 0; % scaling factor for base noise strength
noise_rate = PR/pi;
noise_N = 50; % number of noisy inputs

% define parameters for rectified sinusoidal current injection (if any)
Im_amp = 0;

% Assign values to any parameters specified in function call:
assign(varargin{:});
% typically: PR, F, tmax, dt, Pmax_e, tau1i, delay, alpha, Im_F, Im_amp... 

% Time parameters that are adjusted after manual parameter assignment:
t_total = tmax+tbuffer; % total trial length
tsigvec = 0:dt:tmax; % time vector for signal
tvec = 0:dt:t_total; % time vector for trial
Nt = length(tvec); % number of time intervals

% If sinusoidal direct current injection is present (Im_F and Im_amp are
% assigned in function call), then nonzero Im is generated:
Im = zeros(1,length(tvec));
Im(1:length(tsigvec)) = max(Im_amp * sin(2*pi*F*tsigvec), 0);

    


% Generate a Poisson spike train with input characteristics
spktimes = spiketrain_sinusoidal(PR, F, 0, 0, 0, tmax, dt);
spktrain = spiketimes2bins(spktimes, tvec);

% Generate 50 noisy inputs
noisy_input = poisson_spiketrain(dt, noise_rate, t_total, 1);
for j=2:noise_N
    noisy_input = noisy_input + poisson_spiketrain(dt, noise_rate, t_total, 1);
end
% Scale the noisy input magnitude by original 1.6976e-7/75 (same as in figure 4)
[Ps_E_noise] = ppsc_constantsum(noisy_input, noise_level*noise_base, tau1e, tau2e, 0, tau2i, delay, dt);


% Calculate synaptic conductance due to spike train using
% ppsc_constantsum.m
[Ps_E, Ps_I] = ppsc_constantsum(spktrain, Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);

% Simulate voltage behavior of a neuron recieving push-pull conductance
% input
Vm = zeros(Nt, 1); % store membrane voltage at each time point
Vm(1) = V_reset; % initial membrane voltage is V_reset

tspk = []; % count spikes from Euler's method
post_spktrain = zeros(Nt, 1);


for t = 1:Nt-1 % Euler method
    
    if Vm(t) >= V_th
        Vm(t+1) = V_reset;
        tspk = [tspk; t*dt];
        post_spktrain(t) = 1;
    else
        dvdt = (-(Vm(t) - V_e) - (Rm * ((Ps_E(t) + Ps_E_noise(t)) * (Vm(t) - V_syn_e) + alpha * Ps_I(t) * (Vm(t) - V_syn_i))) + Rm*Im(t)) / tau_m;
        Vm(t+1) = Vm(t) + dt * dvdt;
    end
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

end