% MULTILAYER_SPIKETRACE: Generates figures comparing spiking behavior of 
% 4-layer triad synapse circuit models (Figure 5).

L = 4; % number of layers in the circuit

% Time and input ST parameters
dt = 1e-4;
PR = 100;
F = 5;
t_start = -0.5;
t_end = 1;
tvec = t_start:dt:t_end;
Nt = length(tvec);


% Conductance parameters
Pmax_e = 1.6976e-7;
tau1i = [0.020 0 0];


% Set up vectors to store inputs to each layer and voltage at each layer
% 1: no conductance between layers, 2: e only, 3: pffei
inputs_1 = zeros(L+1,Nt);

inputs_2 = zeros(L+1,Nt);

inputs_3 = zeros(L+1,Nt);


% Generate input to first layer using generate_input.m
input = generate_input(F, 1, 1, 'tbuffer_start', 0.5, 'tmax', 1);

% Store input signal for plotting
inputs_1(1,:) = input.signal;
inputs_2(1,:) = inputs_1(1,:); % provide same input to all cells
inputs_3(1,:) = inputs_1(1,:);

% Prepare noisy input parameters + store noisy inputs in noisy_input
noisy_input = zeros(L,Nt);
noisy_input(1,:) = input.noise;

% Theoretical signal
signal = max(0,PR*sin(2*pi*F*tvec));
signal(1:-t_start/dt) = 0;


% Plot the theoretical input separately, alongside Poisson input signal to 
% first layer

fig1 = figure;
set(fig1,'defaultAxesColorOrder',[[0 0 0]; [0.4 0.4 0]]);
xticks([])
yyaxis left
plot([-0.4 -0.15], [0.75 0.75], 'k', 'LineWidth', 2)
text(-0.35, 0.9, '0.25s')
hold on;
plot(tvec, inputs_1(1,:), '-', 'Color', [0.8 0.8 0])
plot([0 0], [0 1], 'k--', 'LineWidth', 1)
yticks([0 1])
ylabel({'Layer 1 input' '(spikes)'})
hold on;
yyaxis right
plot(tvec, signal)
yticks([0 50 100])
ylabel({'Theoretical' 'firing rate (Hz)'})
title('Input to layer 1')
set(gca, 'LineWidth', 1)
box off;



figure;
subplot(L+1,4,(L+1)*4 - 3); % scalebar for noisy inputs
plot([0.75 1], [100 100], 'k', 'LineWidth', 2)
text(0.75, 125, '0.25s')
ylim([-100 100])
xlim([-0.5 1])
box off;
set(gca,'XColor', 'none','YColor','none')
P = get(gca,'position');
P(1,4) = P(1,4) / 3 * 2;
P(1,2) = P(1,2) + P(1,4)/2;
set(gca,'position',P); 

subplot(L+1,4,(L+1)*4 - 2); % this and below are stimulus plots
plot(tvec, signal, 'Color', [0.8 0.8 0])
hold on;
plot([0 0], [0 100], 'k:', 'LineWidth', 1) % onset of frequency stimulus
ylim([-100 100])
box off;
set(gca,'XColor', 'none','YColor','none')
P = get(gca,'position');
P(1,4) = P(1,4) / 3 * 2;
P(1,2) = P(1,2) + P(1,4)/2;
set(gca,'position',P); 

subplot(L+1,4,(L+1)*4 - 1);
plot(tvec, signal, 'Color', [0.8 0.8 0])
hold on;
plot([0 0], [0 100], 'k:', 'LineWidth', 1) % onset of frequency stimulus
ylim([-100 100])
box off;
set(gca,'XColor', 'none','YColor','none')
P = get(gca,'position');
P(1,4) = P(1,4) / 3 * 2;
P(1,2) = P(1,2) + P(1,4)/2;
set(gca,'position',P); 

subplot(L+1,4,(L+1)*4);
plot(tvec, signal, 'Color', [0.8 0.8 0])
hold on;
plot([0 0], [0 100], 'k:', 'LineWidth', 1) % onset of frequency stimulus
ylim([-100 100])
box off;
set(gca,'XColor', 'none','YColor','none')
P = get(gca,'position');
P(1,4) = P(1,4) / 3 * 2;
P(1,2) = P(1,2) + P(1,4)/2;
set(gca,'position',P);
plot([0.75 1], [-100 -100], 'k', 'LineWidth', 2)
text(0.75, -75, '0.25s')



for i = 1:L
    
    % 1/2 strength: FFEI 0.09, FFE 0.095; 1 strength: FFEI 0.18, FFE 0.19
    input.signal = inputs_1(i,:);
    output_1 = run_triad_model(input, 1, 'Pmax_e', 0, 'tau1i', tau1i(3), 'noise_level', 1);
    input.signal = inputs_2(i,:);
    output_2 = run_triad_model(input, 1, 'Pmax_e', Pmax_e * 0.095, 'tau1i', tau1i(2), 'noise_level', 1);
    input.signal = inputs_3(i,:);
    output_3 = run_triad_model(input, 1, 'Pmax_e', Pmax_e * 0.09, 'tau1i', tau1i(1), 'noise_level', 1);
    
    inputs_1(i+1,:) = output_1.post_spktrain';
    inputs_2(i+1,:) = output_2.post_spktrain';
    inputs_3(i+1,:) = output_3.post_spktrain';
    
    % Plot noisy input to layer i
    subplot(L+1,4,4*i-3);
    plot(tvec, noisy_input(i,:), 'Color', [0.5 0.5 0.5])
    hold on;
    plot([0 0], [0 5], 'k--', 'LineWidth', 1) % onset of frequency stimulus
    ylim([0 5])
    if i == 1
        title({'Summed noisy spike train inputs' 'at each layer (n = 50)'})
    end
    set(gca, 'LineWidth', 1)
    xticks([])
    box off;
    
    % Plot output of layer i (i.e. input to layer i+1)
    subplot(L+1,4,4*i-2);
    plot(tvec,inputs_1(i+1,:), 'Color', [0 0.8 0])
    hold on;
    plot([0 0], [0 1], 'k--', 'LineWidth', 1) % onset of frequency stimulus
    ylabel(['A' num2str(i) ' output'])
    yticks([])
    xticks([])
    set(gca, 'LineWidth', 1)
    if i == 1
        title('No connections')
    end
    box off;
    subplot(L+1,4,4*i-1);
    plot(tvec,inputs_2(i+1,:), 'r')
    hold on;
    plot([0 0], [0 1], 'k--', 'LineWidth', 1) % onset of frequency stimulus
    yticks([])
    xticks([])
    set(gca, 'LineWidth', 1)
    if i == 1
        title('Excitation only')
    end
    box off;
    subplot(L+1,4,4*i);
    plot(tvec,inputs_3(i+1,:), 'Color', [0.4 0.4 1])
    hold on;
    plot([0 0], [0 1], 'k--', 'LineWidth', 1) % onset of frequency stimulus
    yticks([])
    xticks([])
    set(gca, 'LineWidth', 1)
    if i == 1
       title('Paired feed-forward EI')
    end
    box off;
    
    input = generate_input(F, 1, 1, 'tbuffer_start', 0.5, 'tmax', 1);
    noisy_input(i+1,:) = input.noise;
    
end