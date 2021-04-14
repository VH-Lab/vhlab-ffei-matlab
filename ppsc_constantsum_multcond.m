function [ Ps_E, Ps_I ] = ppsc_constantsum_multcond(spike_train, Pmax_e_target, tau1e, tau2e, is_inhibition_present, delay, dt )
%PUSH-PULL SYNAPTIC CONDUCTANCE - CONSTANT SUM & MULTIPLE CONDUCTANCES: Determines the push-pull
%postsynaptic conductance given a spike train. Pmax is scaled so that area under the
%push = area under the pull, and so that "sum" of excitation in a single
%conductance spike is constant for a given set of excitatory conductance
%parameters.
% This version can perform the same operation as ppsc_constantsum for
% conductances with multiple E (NOT I) components (e.g. ampa and nmda); provide
% time constants and target Pmax_e values in vectors. In the case where
% inhibition is present, it is assumed that the inhibition is a delayed
% "copy" of the excitation.

% Compute taurise and B factor for all components of the excitatory conductance
taurise_e = (tau1e .* tau2e)./(tau1e-tau2e);
B_e = ((tau2e./tau1e).^(taurise_e./tau1e) - (tau2e./tau1e).^(taurise_e./tau2e)).^(-1); %normalization factor

% if is_inhibition_present is given as 0, then only excitatory transmission
% will be assumed:
if is_inhibition_present == 0
    
    trial_length = length(spike_train);
    tvec = 0:1:trial_length; % time range
    tvec = tvec * dt; % convert to s
    
    Ps_E = zeros(1, length(spike_train));
    for i = 1:length(tau1e) % sum multiple E conductances
        Pvec_e = [zeros(1, length(tvec)) Pmax_e_target(i) * B_e(i) * (exp(-tvec/tau1e(i)) - exp(-tvec/tau2e(i)))];
        Ps_E = Ps_E + conv(spike_train, Pvec_e, "same");
    end
    
    Ps_I = zeros(1, length(Ps_E));
    
else % otherwise, proceed with push-pull inhibition
    
    % Get "target area" --> area of a positive conductance pulse with no
    % inhibition, given a specified Pmax for excitatory transmission and tau1, tau2
    target_area = 0;
    for i = 1:length(tau1e) % sum area from multiple conductances
        target_area = target_area + Pmax_e_target(i) * B_e(i) * (tau1e(i)-tau2e(i)); % integral of Pmax * B * (exp...)
    end
    
    % Generate a sample conductance response - the area under the inhibitory conductance
    % pulse = the area under the excitatory conductance pulse
    Pmax_e_placeholder = Pmax_e_target ./ Pmax_e_target(1);
    spike_sample = zeros(1, ceil(max(tau1e)/dt) * 20);
    spike_sample(1) = 1;
    
    trial_length = length(spike_sample);
    tvec = 0:1:trial_length; % time range
    tvec = tvec * dt; % convert to s
    
    Pe = zeros(1, length(spike_sample));
    for i = 1:length(tau1e) % sum multiple E conductances
        Pvec_e = [zeros(1, length(tvec)) Pmax_e_placeholder(i) * B_e(i) * (exp(-tvec/tau1e(i)) - exp(-tvec/tau2e(i)))];
        Pe = Pe + conv(spike_sample, Pvec_e, "same");
    end
    
    Pi = [zeros(1,ceil(delay/dt)), Pe]; % inhibition is delayed copy of excitation
    Pi = Pi(1:end-ceil(delay/dt));

    P_sample = Pe - Pi;
    
    % Determine the positive area under P_sample
    positive_area = sum(P_sample(P_sample>0))* dt;
    
    % Divide the target area by the positive area to get a multiplier for
    % placeholder Pmax_e
    Pmax_multiplier = target_area/positive_area;
    
    % Multiply placeholder Pmax_e by multiplier to get Pmax_e such that
    % excitation sum is consistent
    Pmax_e = Pmax_multiplier * Pmax_e_placeholder;
    
    
    %Determine postsynaptic conductance from spike train
    
    trial_length = length(spike_train);
    tvec = 0:1:trial_length; % time range
    tvec = tvec * dt; % convert to s
    
    
    Ps_E = zeros(1, length(spike_train));
    for i = 1:length(tau1e) % sum multiple E conductances
        Pvec_e = [zeros(1, length(tvec)) Pmax_e(i) * B_e(i) * (exp(-tvec/tau1e(i)) - exp(-tvec/tau2e(i)))];
        Ps_E = Ps_E + conv(spike_train, Pvec_e, "same");
    end

    Ps_I = [zeros(1,ceil(delay/dt)), Ps_E]; % inhibition is delayed copy of excitation
    
    % modify length of inhibitory conductance vector:
    Ps_I = Ps_I(1:end-ceil(delay/dt));
    
    
end

end

