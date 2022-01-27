function [ Ps ] = synaptic_conductance( spike_train, Pmax, tau1, tau2, dt, delay)
%SYNAPTIC CONDUCTANCE: Determines the postsynaptic conductance given a
%spike train
%   dt is given in seconds

if ~exist('delay', 'var')
    delay = 0;
end

trial_length = length(spike_train);
tvec = 0:1:trial_length; %time range (ms)
tvec = tvec * dt;

taurise = (tau1 * tau2)/(tau1-tau2);
B = ((tau2/tau1)^(taurise/tau1) - (tau2/tau1)^(taurise/tau2))^(-1); %normalization factor

Pvec = [zeros(1, length(tvec)) Pmax * B * (exp(-tvec/tau1) - exp(-tvec/tau2))];
Ps = conv(spike_train, Pvec, "same");

% incorporate delay (if not zero)
Ps = [zeros(1,ceil(delay/dt)), Ps];
Ps = Ps(1:end-ceil(delay/dt));


end


