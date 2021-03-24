function [ Ps_E, Ps_I ] = ppsc_constantsum(spike_train, Pmax_e_target, tau1e, tau2e, tau1i, tau2i, delay, dt )
%PUSH-PULL SYNAPTIC CONDUCTANCE - CONSTANT SUM: Determines the push-pull
%postsynaptic conductance given a spike train; takes time constants for
%excitation and inhibition separately. Pmax is scaled so that area under the
%push = area under the pull, and so that sum of excitation in a single
%conductance spike is constant for a given set of excitatory conductance
%parameters.
%   dt is given in seconds

trial_length = length(spike_train);
tvec = 0:1:trial_length; % time range
tvec = tvec * dt; % convert to s

taurise_e = (tau1e * tau2e)/(tau1e-tau2e);
%tau2e = taurise_e * (1 + taurise_e / tau1e);
B_e = ((tau2e/tau1e)^(taurise_e/tau1e) - (tau2e/tau1e)^(taurise_e/tau2e))^(-1); %normalization factor



% if tau1i is given as 0, then only excitatory transmission will be
% assumed:
if tau1i ==0
    Pvec_e = [zeros(1, length(tvec)) Pmax_e_target * B_e * (exp(-tvec/tau1e) - exp(-tvec/tau2e))];
    Ps_E = conv(spike_train, Pvec_e, "same");
    
    Ps_I = zeros(1, length(Ps_E));
else % otherwise, proceed with push-pull inhibition
    
    
    taurise_i = (tau1i * tau2i)/(tau1i -tau2i);
    %tau2i = taurise_i * (1 + taurise_i / tau1i);
    B_i = ((tau2i/tau1i)^(taurise_i/tau1i) - (tau2i/tau1i)^(taurise_i/tau2i))^(-1); %normalization factor
    
    % Get "target area" --> area of a positive conductance pulse with no
    % inhibition, given a specified Pmax for excitatory transmission and tau1, tau2
    target_area = Pmax_e_target * B_e * (tau1e-tau2e); % integral of Pmax * B * (exp...)
    
    % Generate a sample conductance response using
    % pushpull_synaptic_conductance - the area under the inhibitory conductance
    % pulse = the area under the excitatory conductance pulse for this
    % function
    Pmax_e_placeholder = 1;
    spike_sample = zeros(1, ceil(tau1e/dt) * 100);
    spike_sample(1) = 1;
    [Pe, Pi] = pushpull_synaptic_conductance(spike_sample, Pmax_e_placeholder, tau1e, tau2e, tau1i, tau2i, delay, dt);
    P_sample = Pe - Pi;
    
    % Determine the positive area under P_sample
    positive_area = sum(P_sample(P_sample>0))* dt;
    
    % Divide the target area by the positive area to get a multiplier for
    % placeholder Pmax_e
    Pmax_multiplier = target_area/positive_area;
    
    % Multiply placeholder Pmax_e by multiplier to get Pmax_e such that
    % excitation sum is consistent
    
    Pmax_e = Pmax_multiplier * Pmax_e_placeholder;
    
    
    
    % Scale Pmax_i so that, given tau1i and tau2i different from tau1e and
    % tau2e, the area under the inhibitory conductance pulse = the area under
    % the excitatory conductance pulse
    Pmax_i = (Pmax_e * B_e / B_i) * (tau1e-tau2e)/(tau1i-tau2i);
    
    
    %Determine postsynaptic conductance from spike train
    
    Pvec_e = [zeros(1, length(tvec)) Pmax_e * B_e * (exp(-tvec/tau1e) - exp(-tvec/tau2e))];
    Ps_E = conv(spike_train, Pvec_e, "same");
    
    Pvec_i = [zeros(1, length(tvec)) Pmax_i * B_i * (exp(-tvec/tau1i) - exp(-tvec/tau2i))];
    Ps_I = conv(spike_train, Pvec_i, "same");
    Ps_I = [zeros(1,ceil(delay/dt)), Ps_I]; % inhibitory conductance is negative
    % and activates with a delay.
    
    % modify length of excitatory conductance vector:
    % Ps_E = [Ps_E, zeros(1,ceil(delay/dt))];
    
    % modify length of inhibitory conductance vector:
    Ps_I = Ps_I(1:end-ceil(delay/dt));
    
    %Ps_I(:) = 0; % set inhibitory conductance to 0 so that only excitatory
    %conductance is present in the system (only for excitatory transmission)
    
    % Ps = Ps_ E + Ps_I;
    
end

end

