% INJECTION_VOLTAGETRACE: Generates figures comparing spiking behavior of
% model neurons in response to rectified sinusoidal current injections of
% differing frequency.

dt = 1e-4;

% Traces for sinewave current
Imex = 8e-9;
Fex = [8 50 100]; % example current frequencies
siglength = 0.25; % length of current injection
sigvec = 0:dt:siglength;

tmax = 0.4;
tvec = 0:dt:tmax;

sinewave_spktrain_F = zeros(length(Fex), length(tvec)); % store output spike trains with variable input frequency
sinewave_V_F = zeros(length(Fex), length(tvec)); % store voltage trace with variable input frequency
sinewave_Im_F = zeros(length(Fex), length(tvec)); % store current signals for variable input frequency

for i = 1:length(Fex) % get spike trains and voltage trace for variable durations
   input = generate_input(Fex(i), 1, 0, 'PR', 0, 'tmax', siglength, 'tbuffer_end', tmax-siglength, 'Im_amp', Imex);
   sinewave_Im_F(i,:) = input.Im;
   output = run_triad_model(input, 1);
   
   post_spktrain = output.post_spktrain;
   sinewave_spktrain_F(i,:) = post_spktrain(:);
   sinewave_V_F(i,:) = output.Vm;
   
end

% Plot injected current and spiking outputs for each frequency

figure;

subplot(3,1,1);
plot(tvec, (sinewave_Im_F(1,:)*1e7)-0.08, 'k:')
hold on;
plot(tvec, sinewave_V_F(1,:), 'k')
for i = 1:length(sinewave_spktrain_F(1,:))
   if sinewave_spktrain_F(1,i) == 1
   line([i*dt i*dt],[-0.08 0.0], 'Color', 'k', 'LineWidth', 1)
   end
end
plot([0 tmax], [-0.1 -0.1], 'Color', 'black', 'HandleVisibility','off')
set(gca, 'Visible', 'off')
ylim([-0.1 0.0])
ylabel(["Membrane", "potential (V)"])
legend('8Hz current injection', 'Location', 'southeast')
%title("Spike train of a LIF neuron given Input 1")

plot([0.285 0.285], [-0.04 0], 'Color', 'black', 'HandleVisibility','off')
plot([0.285 0.31], [-0.04 -0.04], 'Color', 'black', 'HandleVisibility','off')
plot([0.31 0.31], [-0.04 0], ':', 'Color', 'black', 'HandleVisibility','off')

txt_x = '25ms';
txt_y = '40mV';
txt_y2 = '4nA';
text(0.29, -0.06, txt_x);
text(0.263, -0.02, txt_y);
text(0.31, -0.02, txt_y2);


subplot(3,1,2);
plot(tvec, (sinewave_Im_F(2,:)*1e7)-0.08, 'k:')
hold on;
plot(tvec, sinewave_V_F(2,:), 'k')
hold on;
for i = 1:length(sinewave_spktrain_F(2,:))
   if sinewave_spktrain_F(2,i) == 1
   line([i*dt i*dt],[-0.08 0.0], 'Color', 'k', 'LineWidth', 1)
   end
end
plot([0 tmax], [-0.1 -0.1], 'Color', 'black', 'HandleVisibility','off')
set(gca, 'Visible', 'off')
ylim([-0.1 0.0])
ylabel(["Membrane", "potential (V)"])
legend('50Hz current injection', 'Location', 'southeast')
%title("Spike train of a LIF neuron given Input 2")

plot([0.285 0.285], [-0.04 0], 'Color', 'black', 'HandleVisibility','off')
plot([0.285 0.31], [-0.04 -0.04], 'Color', 'black', 'HandleVisibility','off')
plot([0.31 0.31], [-0.04 0], ':', 'Color', 'black', 'HandleVisibility','off')

txt_x = '25ms';
txt_y = '40mV';
txt_y2 = '4nA';
text(0.29, -0.06, txt_x);
text(0.263, -0.02, txt_y);
text(0.31, -0.02, txt_y2);



subplot(3,1,3);
plot(tvec, (sinewave_Im_F(3,:)*1e7)-0.08, 'k:')
hold on;
plot(tvec, sinewave_V_F(3,:), 'k')
hold on;
for i = 1:length(sinewave_spktrain_F(3,:))
   if sinewave_spktrain_F(3,i) == 1
   line([i*dt i*dt],[-0.08 0.0], 'Color', 'k', 'LineWidth', 1)
   end
end
plot([0 tmax], [-0.1 -0.1], 'Color', 'black', 'HandleVisibility','off')
set(gca, 'Visible', 'off')
ylim([-0.1 0.0])
xlabel("Time (sec)")
ylabel(["Membrane", "potential (V)"])
legend('100Hz current injection', 'Location', 'southeast')
%title("Spike train of a LIF neuron given Input 3")

plot([0.285 0.285], [-0.04 0], 'Color', 'black', 'HandleVisibility','off')
plot([0.285 0.31], [-0.04 -0.04], 'Color', 'black', 'HandleVisibility','off')
plot([0.31 0.31], [-0.04 0], ':', 'Color', 'black', 'HandleVisibility','off')

txt_x = '25ms';
txt_y = '40mV';
txt_y2 = '4nA';
text(0.29, -0.06, txt_x);
text(0.263, -0.02, txt_y);
text(0.31, -0.02, txt_y2);