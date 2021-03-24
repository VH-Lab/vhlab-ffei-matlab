function [] = triad_voltagetrace(transmissiontype)
% TRIAD_VOLTAGETRACE: Compares spiking behavior of FFE vs. FFEI triad
% synapse circuit models in response to a Poisson input with rectified
% sinusoidal rate (PR = 100Hz, F = 50Hz) (Figure 1).

% Input characteristics
PR = 100;  % peak rate of neuron recieving input
F = 50;   % frequency of input
dt = 1e-4; % time bins
tmax = 0.5;
num_inputs = 1;

tvec = 0:dt:tmax;
rate = sin(2 * pi * tvec * F) * PR;

% Generate Poisson spike train using input characteristics
%spktrain = poisson_spiketrain_wave(dt, tmax, F, PR, num_inputs);
spktimes = spiketrain_sinusoidal(PR, F, 0, 0, tvec(1), tvec(end), dt);
spktrain = spiketimes2bins(spktimes, tvec);

% Conductance parameters
%Pmax_e = 6.7903;
if transmissiontype == 1
    Pmax_e = 1.6976e-7*0.47; % E only
    tau1i = 0; % 50ms
elseif transmissiontype == 2
    Pmax_e = 1.6976e-7*0.305; % FFEI
    tau1i = 0.02; % 50ms
else
    error('Error: Transmission type not correctly defined. Set input to 1 if examining FFE model, 2 if examining FFEI model.')
end
tau1e = 0.02; % 50ms % taufall
tau2e = 0.001; % 1ms % taurise
tau2i = tau2e; % 1ms
delay = 0.001; % 1ms

alpha = 1.25; % inhibitory current scaling constant

% tau1 is fall time constant, tau2 is rise time constant
% delay is the time between push and pull

% push-pull conductance triggered by Poisson spike train
[Ps_E, Ps_I] = ppsc_constantsum(spktrain(1,:), Pmax_e, tau1e, tau2e, tau1i, tau2i, delay, dt);


% Simulate voltage behavior of a neuron recieving push-pull conductance
% input

% define integrate fire variables
V_reset = -0.080; % -80mV
V_e = -0.075; % -75mV
V_th = -0.040; % -40mV
Rm = 1.0e7; % membrane resistance
tau_m = 1.0e-2; % time
V_syn_e = 0; % synaptic resting potential (excitatory)
V_syn_i = -0.08; % synaptic resting potential (inhibitory)
rmgs = 0.25; % original: 0.05
Nt = length(tvec); % number of time intervals

Im = zeros(Nt, 1); % input current

Vm = zeros(Nt, 1); % store membrane voltage at each time point
Vm(1) = V_reset; % initial membrane voltage is V_reset

post_spktrain = zeros(Nt, 1); % count spikes from Euler's method


for t = 1:Nt-1 % Euler method
    
    if Vm(t) >= V_th
        Vm(t+1) = V_reset;
        post_spktrain(t) = 1;
    else
        dvdt = (-(Vm(t) - V_e) - (Rm * (Ps_E(t) * (Vm(t) - V_syn_e) + alpha * Ps_I(t) * (Vm(t) - V_syn_i))) + Rm*Im(t)) / tau_m;
        Vm(t+1) = Vm(t) + dt * dvdt;
    end
end

% sum conductances of multiple input neurons

spiketimes = find(post_spktrain) * dt;

postrate = post_spktrain/dt;
postrate_avg = mean(postrate);


% For poster

figure;
subplot(3,1,1)
plot(tvec, rate, 'Color', [1 0.733 0.318])
hold on;
for i = 1:num_inputs
    plot(tvec, spktrain(i,:) * PR / 2, 'Color', [0.004 0.514 0])
end
ylim([0 PR*2])
set(gca, 'Visible', 'off')
plot([0 tmax], [0 0], 'Color', 'black', 'HandleVisibility','off')
legend("Rate", "Spike present", 'Location', 'northwest')
legend(gca,'show')

scaletime1 = 0.215;
scaletime2 = 0.24;
scalerate1 = 150;
scalerate2 = 175;

plot([scaletime1 scaletime1], [scalerate1 scalerate2], 'Color', 'black', 'HandleVisibility','off') %vertical
plot([scaletime1 scaletime2], [scalerate1 scalerate1], 'Color', 'black', 'HandleVisibility','off') %horizontal

txt_x = [num2str((scaletime2-scaletime1)*1000),'ms'];
txt_y = [num2str(scalerate2-scalerate1),'Hz'];
text(0.218, 135, txt_x);
text(0.1925, 163, txt_y);

% title("Poisson spike train for periodic input")
% xlabel("Time (s)")


subplot(3,1,2)
plot(tvec, max(Ps_E, 6e-9), 'Color', [0.541 0.824 0.325])
hold on;
plot(tvec, max(Ps_I, 6e-9), 'Color', [0.639 0.639 0.639])
set(gca, 'Visible', 'off')
ylim([0 max(Ps_E)*2])
plot([0 tmax], [0 0], 'Color', 'black', 'HandleVisibility','off')
legend("G_s_y_n_,_e", "G_s_y_n_,_i", 'Location', 'northwest')
legend(gca,'show')

units = 10^floor(log10(max(Ps_E))); % conductance units for scale bar
top = round(max(Ps_E),1,'significant');

scalecond1 = top + 1*units;
scalecond2 = top + 3*units;

plot([scaletime1 scaletime1], [scalecond1 scalecond2], 'Color', 'black', 'HandleVisibility','off') %vertical
plot([scaletime1 scaletime2], [scalecond1 scalecond1], 'Color', 'black', 'HandleVisibility','off') %horizontal

txt_x = [num2str((scaletime2-scaletime1)*1000),'ms'];
txt_y = [num2str((units)*10^6 * 2),'\muS'];
text(0.218, (scalecond1 - units), txt_x);
text(0.19, (scalecond2 - units), txt_y);


subplot(3,1,3)
plot(tvec, rate * (0.15/PR) - 0.15, ':', 'Color', [1 0.733 0.318])
hold on;
plot(tvec, Vm, 'Color', [0.396 0.569 1])
xlim([0 max(tvec)])
ylim([-0.15 0.15])
if isempty(spiketimes) == 0
    plot([spiketimes(:) spiketimes(:)], [-0.08 0], 'Color', [0.396 0.569 1], 'HandleVisibility', 'off')
    plot(spiketimes, 0, 'x', 'Color', 'red')
end

plot([0 tmax], [-0.15 -0.15], 'Color', 'black', 'HandleVisibility','off')
set(gca, 'Visible', 'off')
legend("Input rate", "Membrane potential", "Spike present", 'Location', 'northwest')
legend(gca,'show')

scalev1 = 0.075;
scalev2 = 0.125;

plot([scaletime1 scaletime1], [scalev1 scalev2], 'Color', 'black', 'HandleVisibility','off') %vertical
plot([scaletime1 scaletime2], [scalev1 scalev1], 'Color', 'black', 'HandleVisibility','off') %horizontal

txt_x = [num2str((scaletime2-scaletime1)*1000),'ms'];
txt_y = [num2str((scalev2-scalev1)*1000),'mV'];
text(0.218, 0.05, txt_x);
text(0.192, 0.1, txt_y);
end