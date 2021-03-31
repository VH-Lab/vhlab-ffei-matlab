function[spktrain] = alpha_voltagetrace(alpha, Pmax_e, tau1i)
% ALPHA_VOLTAGETRACE: Generates figures comparing spiking behavior of FFEI
% triad synapse circuit models with different alpha (FFE, 1, 1.25, 2, 5)

% General simulation parameters
F = 50; % input modulation frequency (Hz)
% tmax = 0.2; % length of simulation (s)
tmax = 0.5;
dt = 0.0001/5; % time step length (s)
tvec = 0:dt:tmax; % time vector (s)
Nt = length(tvec);

% provide a custom input: 50 Hz uniform spiking signal
spktrain = zeros(1,Nt);
spktrain(round(1/(4*F)/dt):round(1/F/dt):end) = 1;

% Run triad synapse simulations using run_triad_model.m
% spktrain: Poisson spike train input over time
% Ps_E, Ps_I: excitatory (or inhibitory) conductance post-synapse due to
%   spktrain
% post_spktrain: spiking output over time
% Vm: membrane potential of output neuron over time 
[~,~,~,~,Vm_noreset,Ps_E,Ps_I,input, I_syn_noreset, I_leak_noreset] = run_triad_model(0, 'manual_signal', spktrain, 'tmax', tmax, 'dt', dt, 'Pmax_e', Pmax_e, 'tau1i', tau1i, 'alpha', alpha);

% Plot model behavior without Vm resets
figure;

subplot(4,1,1)
plot(tvec, input, 'k')
hold on;
%plot(tvec, rate/PR, 'k:')
ylabel('Spike present')
xlabel('Time (s)')
box off;
title(['alpha = ' num2str(alpha)])

subplot(4,1,2)
plot(tvec, Ps_E, 'r')
hold on;
plot(tvec, Ps_I, 'b')
ylabel('Synaptic conductance (S)')
xlabel('Time (s)')
legend('G_e', 'G_i')
box off;

subplot(4,1,3)
hold on;
%plot(tvec, Im_syn_noreset + Im_leak_noreset, 'k')
plot(tvec, I_syn_noreset, 'k')
plot(tvec, -I_leak_noreset, 'g')
%plot([0 tvec(end)], [0 0], 'k:')
ylabel('Current (A)')
ylim([-2e-8 5e-8])
xlabel('Time (s)')
legend('I_s_y_n', '-I_l_e_a_k')
box off;

subplot(4,1,4)
plot(tvec, Vm_noreset, 'k')
ylabel('Voltage (V)')
ylim([-0.09 -0.02])
xlabel('Time (s)')
box off;

end