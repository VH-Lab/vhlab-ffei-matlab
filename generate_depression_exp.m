function [spiketimes,syncurrs] = generate_depression_exp(spikeintervals,ppratio)
%GENERATE DEPRESSION EXPERIMENT: Generates synaptic currents at each
%spike time for a simulated paired pulse depression experiment, to provide
%as input to depression_model_fit.m.
% spikeintervals: vector of intervals between s1 and s2
% ppratio: vector of paired-pulse ratios (s2/s1)

spiketimes = 1:2*length(spikeintervals);
syncurrs = spiketimes;

for i = 1:length(spikeintervals)
    start_time = 30 * (i-1); % separate spike pairs by 30s
    spiketimes(2*i-1) = start_time;
    spiketimes(2*i) = spikeintervals(i) + start_time;
    syncurrs(2*i-1) = 1;
    syncurrs(2*i) = ppratio(i);
end

end

