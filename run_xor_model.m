function [output] = run_xor_model(input_1, input_2, Pmax_e, is_inhibition_present, FC_freq)
% RUN_XOR_MODEL: Performs simulations of neuromorphic circuit that performs
% an XOR computation, in response to two inputs specified by
% generate_input.m input structs (Figure 6).
% Returns FC at a specified frequency (FC_freq), FC over all frequencies,
% normalized FC, Poisson spike train output of each cell in the circuit,
% input spike trains, and expected output (generated by performing XOR
% computation on input trains in sliding 5ms (or 10ms) bins; XOR 
% computation = 1 if only one input is active during a given bin)
% Additional parameters:
%   Pmax_e: transmission strength throughout the circuit (rescaled in FFEI
%           case)
%   is_inhibition_present: Determines whether FFE transmission (0) or FFEI
%                          tranmission (1) is used in projections within
%                          the circuit

PR = 100;
dt = input_1.dt;
tvec = input_1.tvec;

num_inputs = 2;

% For neurons with different phases
input1 = input_1.signal;
input2 = input_2.signal;
inputs = [input1;input2];

% Conductance parameters
tau1e = 0.02; % 50ms % taufall
tau2e = 0.001; % 1ms % taurise
if is_inhibition_present
    tau1i = 0.02;
else
    tau1i = 0;
end
tau2i = 0.001; % 1ms
delay = 0.001; % 1ms
alpha = 1.25; % inhibitory scaling coefficient

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
Cm = 1.0e-9; % membrane capacitance
tau_m = Rm*Cm; % time
V_syn_e = 0; % synaptic resting potential (excitatory)
V_syn_i = -0.08; % synaptic resting potential (inhibitory)
Nt = length(tvec); % number of time intervals

% RULES: weights for LTU - exclusive or

% Some notes on terminology...
% gate: output cell at each layer; where inputs converge and are integrated

num_gates_1 = 2; % in layer 1, there are 2 gate output cells, both of which receive inputs 1 and 2
% layer 1 weights
w1 = [1 -1; -1 1]; % row 1 is weights for gate 1
                   % row 2 is weights for gate 2

num_gates_2 = 1;
% layer 2 weights
w2 = [1 1]; % row 1 is weights for gate 3


% Synaptic conductance from inputs, given that, following a spike:
%   1) Area under summed excitatory conductance pulse is constant for any
%   set of time constants
%   2) Area under independent inhibitory conductance pulse is equal to area
%   under independent excitatory conductance pulse


Psvec_E = zeros(num_inputs, Nt);
Psvec_I = Psvec_E;

for n = 1:num_inputs
    [Ps_E, Ps_I] = ppsc_constantsum(inputs(n,:), Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);
    
    Psvec_E(n, :) = Ps_E;
    Psvec_I(n, :) = Ps_I;
    
end


% get outputs from first layer
Vm = zeros(num_gates_1, Nt); % store membrane voltage at each time point
Vm(:,1) = V_reset; % initial membrane voltage is V_reset

% count spikes from Euler's method
spktrain1 = zeros(num_gates_1, Nt);

for n = 1:num_gates_1
    for t = 1:Nt-1 % Euler method
        
        if Vm(n,t) >= V_th
            Vm(n,t+1) = V_reset;
            spktrain1(n,t) = 1;
        else
            PsE = 0;
            PsI = 0;
            for input = 1:num_inputs
                if w1(n,input) > 0 % if weight on input is positive, delayed inhibition
                    PsE = PsE + Psvec_E(input,t) * w1(n,input);
                    PsI = PsI + Psvec_I(input,t) * w1(n,input);
                else % otherwise, delayed excitation
                    PsI = PsI + Psvec_E(input,t) * -w1(n,input);
                    PsE = PsE + Psvec_I(input,t) * -w1(n,input);
                end
            end
            dvdt = (-(Vm(n,t) - V_e) - (Rm * (PsE * (Vm(n,t) - V_syn_e) + alpha * PsI * (Vm(n,t) - V_syn_i)))) / tau_m;
            Vm(n,t+1) = Vm(n,t) + dt * dvdt;
        end
    end
end

% synaptic conductance to second layer
Psvec_E = zeros(num_gates_1, Nt);
Psvec_I = Psvec_E;

for n = 1:num_gates_1
    [Ps_E, Ps_I] = ppsc_constantsum(spktrain1(n,:), Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);

    Psvec_E(n,:) = Ps_E;
    Psvec_I(n,:) = Ps_I;
    
end


% get outputs from second layer
Vm = zeros(num_gates_2, Nt); % store membrane voltage at each time point
Vm(:,1) = V_reset; % initial membrane voltage is V_reset

% count spikes from Euler's method
spktrain2 = zeros(num_gates_2, Nt);

for n = 1:num_gates_2
    for t = 1:Nt-1 % Euler method
        
        if Vm(n,t) >= V_th
            Vm(n,t+1) = V_reset;
            spktrain2(n,t) = 1;
        else
            PsE = 0;
            PsI = 0;
            for input = 1:num_inputs
                if w2(n,input) > 0 % if weight on input is positive, delayed inhibition
                    PsE = PsE + Psvec_E(input,t) * w2(n,input);
                    PsI = PsI + Psvec_I(input,t) * w2(n,input);
                else % otherwise, delayed excitation
                    PsI = PsI + Psvec_E(input,t) * -w2(n,input);
                    PsE = PsE + Psvec_I(input,t) * -w2(n,input);
                end
            end
            dvdt = (-(Vm(n,t) - V_e) - (Rm * (PsE * (Vm(n,t) - V_syn_e) + alpha * PsI * (Vm(n,t) - V_syn_i)))) / tau_m;
            Vm(n,t+1) = Vm(n,t) + dt * dvdt;
        end
    end
end

postratevec = spktrain2(1,:)/dt; % Convert spike train to rate over time

% Calculate fourier coefficients for frequency specified by FC_freq
FC = abs(fouriercoeffs_tf2(postratevec(:), FC_freq, 1/dt));
FC_avg = mean(abs(fouriercoeffs(postratevec(:), dt)));

if FC_avg == 0 % set all FCs to zero if FC_avg = 0
    FC_pct = 0;
else
    FC_pct = FC/FC_avg;
end


% calculate expected exclusive-or output given 2 inputs
% a sliding 5ms (or 10ms) window will be used to count number of spikes over time
% for the 2 inputs --> if only one of the inputs is active during a given
% window, then exclusive-or is true; if neither or both inputs are active,
% then exclusive-or is false.

window_5ms = ones(1, ceil(0.005/dt)); % 5ms windows

input1_5ms = min(conv(inputs(1,:), window_5ms, 'same'),1);
input2_5ms = min(conv(inputs(2,:), window_5ms, 'same'),1);

xor_expected_5ms = xor(input1_5ms, input2_5ms);

spk_xor_expected_5ms = xor_expected_5ms; % display expected xor output
spk_xor_expected_5ms(1, 1:2:end) = 0;


window_10ms = ones(1, ceil(0.01/dt)); % 10ms windows

input1_10ms = min(conv(inputs(1,:), window_10ms, 'same'),1);
input2_10ms = min(conv(inputs(2,:), window_10ms, 'same'),1);

xor_expected_10ms = xor(input1_10ms, input2_10ms);

spk_xor_expected_10ms = xor_expected_10ms; % display expected xor output
spk_xor_expected_10ms(1, 1:2:end) = 0;


output = var2struct('FC', 'FC_avg', 'FC_pct', 'spktrain1', 'spktrain2', 'inputs', 'spk_xor_expected_10ms', 'spk_xor_expected_5ms');

end

