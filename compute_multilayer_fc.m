function [FC,FC_avg,FC_pct] = compute_multilayer_fc(F, L, Pmax_factor, PFFEI)
% COMPUTE_MULTILAYER_FC: Computes output power (FC, FCavg, FCpct) for a single multilayer
% simulation.

% Set up time parameters
dt = 1e-4;
t_start = 0;
t_end = 1;
tvec = t_start:dt:t_end;
Nt = length(tvec); % number of time intervals


% Conductance parameters
Pmax_e = 1.6976e-7*Pmax_factor;
if PFFEI
    tau1i = 0.02; % 20ms % set to 0 to simulate only excitatory transmission
else
    tau1i = 0;
end


inputs = zeros(L+1,Nt);
inst_input = generate_input(F, 1, 1, 'tmax', t_end); % generate input to first layer (inst_input will be used to store all input features at current layer in for loop)
inputs(1,:) = inst_input.signal;

for i = 1:L
    
    inst_input.signal = inputs(i,:);
    output_data = run_triad_model(inst_input, 1, 'Pmax_e', Pmax_e, 'tau1i', tau1i, 'noise_level', 1);
    inputs(i+1,:) = output_data.post_spktrain;
    
    inst_input = generate_input(F, 1, 1, 'tmax', t_end); % generate (noisy) input to next layer (during (i+1)th iteration of the loop, we assign the output of layer i to the input signal)
    
end

output = inputs(L+1, :)/dt; % Extract final layer output and convert to rate
FC = abs(fouriercoeffs_tf2(output(:), F, 1/dt));

FC_all = fouriercoeffs(output(:), dt);
FC_avg = mean(abs(FC_all));
FC_pct = FC/FC_avg;

end

