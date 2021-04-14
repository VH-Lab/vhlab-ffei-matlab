function [] = cortex_voltagetrace(F)

% CORTEX_VOLTAGETRACE: Generates figures comparing spiking behavior of FFEI
% multi-input circuit models in response to inputs with different
% rectified sinusoidal modulation rate (PR = 100Hz) (Figure 8).

% Time parameters
t_end = 0.6;

% define necessary LIF parameters - I cell
tau_m_I = 1.0e-2;

% Run cortex model using default parameters, given input signal with
% frequency F
output = run_cortex_model(F, tau_m_I, 't_end', t_end);

% Extract model behavior from output struct
Ps_E_total = output.Ps_E_total;
Ps_I_total = output.Ps_I_total;
Ps_I_output = output.Ps_I_output;
output_I = output.output_I;
output_E = output.output_E;
tvec = output.tvec;


% Generate plots of model behavior
figure;
subplot(4,1,1)
plot(tvec, Ps_I_total)
title('Input conductance to I')
box off;
subplot(4,1,2)
plot(tvec, output_I)
title('I output')
box off;
subplot(4,1,3)
plot(tvec, Ps_E_total)
hold on;
plot(tvec, Ps_I_output)
title('Input conductance to E')
legend('Summed, scaled excitatory', 'Inhibitory from I')
box off;
subplot(4,1,4)
plot(tvec, output_E)
title('E output')
box off;
xlabel('Time (s)')

end