function [Pmax_i_ratio] = generate_balanced_EI(Pmax_e, tau1e, tau2e, Pmax_i_ratio, tau1i, tau2i)
%GENERATE BALANCED E/I: For a triad-type synapse with multiple excitatory
%and inhibitory channels with known temporal parameters (tau), computes
%the inhibitory synapse strength (Pmax_i) such that excitation (with known
%strength(s)) is balanced by inhibition.
%   Assumes 1 inhibitory channel type (GABA), any number of excitatory
%   channel types

% Compute taurise and B factor for all components of the E/I conductance
taurise_e = (tau1e .* tau2e)./(tau1e-tau2e);
B_e = ((tau2e./tau1e).^(taurise_e./tau1e) - (tau2e./tau1e).^(taurise_e./tau2e)).^(-1); %normalization factor

taurise_i = (tau1i .* tau2i)./(tau1i-tau2i);
B_i = ((tau2i./tau1i).^(taurise_i./tau1i) - (tau2i./tau1i).^(taurise_i./tau2i)).^(-1); %normalization factor


% Get "target area": area of a positive conductance pulse with no
% inhibition, given a specified Pmax for E transmission and tau1, tau2
target_area_e = 0;
for i = 1:length(tau1e) % sum area from multiple conductances
    target_area_e = target_area_e + Pmax_e(i) * B_e(i) * (tau1e(i)-tau2e(i)); % integral of Pmax * B * (exp...)
end

target_area_i = 0;
for i = 1:length(tau1i) % sum area from multiple conductances
    target_area_i = target_area_i + Pmax_i_ratio(i) * B_i(i) * (tau1i(i)-tau2i(i)); % integral of Pmax * B * (exp...)
end

% Set Pmax_i such that the area under the inhibitory conductance = target
% area from excitatory conductances

Pmax_i_ratio = target_area_e / target_area_i * Pmax_i_ratio;

end

