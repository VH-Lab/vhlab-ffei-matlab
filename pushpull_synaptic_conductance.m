function [ Ps_E, Ps_I ] = pushpull_synaptic_conductance(spike_train, Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt )
%SYNAPTIC CONDUCTANCE: Determines the push-pull postsynaptic conductance given a
%spike train; takes time constants for excitation and inhibition
%separately, and scales Pmax_I so that area under the push = area under the
%pull.
%   dt is given in seconds

trial_length = length(spike_train);
tvec = 0:1:trial_length; % time range
tvec = tvec * dt; % convert to s

taurise_e = (tau1e * tau2e)/(tau1e-tau2e);
%tau2e = taurise_e * (1 + taurise_e / tau1e);
B_e = ((tau2e/tau1e)^(taurise_e/tau1e) - (tau2e/tau1e)^(taurise_e/tau2e))^(-1); %normalization factor


taurise_i = (tau1i * tau2i)/(tau1i -tau2i);
%tau2i = taurise_i * (1 + taurise_i / tau1i);
B_i = ((tau2i/tau1i)^(taurise_i/tau1i) - (tau2i/tau1i)^(taurise_i/tau2i))^(-1); %normalization factor

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

end

