% MULTILAYER_SPIKETRACE: Generates figures comparing spiking behavior of 
% 4-layer triad synapse circuit models (Figure 5).
% NOTE: this one is a bit trickier to adapt: run_triad_model.m provides
% randomized inputs (background and feed-forward) during each trial, but in
% the spike traces shown in Figure 5A, all models (3 trials) recieve the
% same inputs. Perhaps provide custom inputs for both background noise and
% signal, and generate all inputs in this script? For now, plotting of the
% outputs from each layer also occurs in this script.

