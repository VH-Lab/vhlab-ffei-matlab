function [] = xor_spiketrace(F, phase_shift)
% XOR_VOLTAGETRACE: Given a certain modulation frequency and phase shift
% between 2 inputs, plots the spiking behavior of a neuromorphic XOR
% model circuit with 1) FFE or 2) FFEI projections between units (Figure
% 6).
% Note that this function currently only balances transmission between the
% 2 models for phase shifts = pi/2 or pi.

% Time parameters
PR = 100;
dt = 1e-4/5;
tmax = 0.6;
tvec = 0:dt:tmax;

% Generate phase-shifted inputs
input1 = generate_input(F, 1, 0, 'phase_shift', 0, 'dt', dt, 'tmax', tmax);
input2 = generate_input(F, 1, 0, 'phase_shift', phase_shift, 'dt', dt, 'tmax', tmax);

% Theoretical signals
signal1 = max(0,PR*sin(2*pi*F*tvec));
signal2 = max(0,PR*sin(2*pi*F*tvec + phase_shift));

% Set conductance parameters
is_inhibition_present = [0 1];
Pmax_base = 1.6976e-7;
if phase_shift == pi
    Pmax_e = Pmax_base * [0.39 0.33]; % 2xF parameters
    FC_freq = 2*F;
elseif phase_shift == pi/2
    Pmax_e = Pmax_base * [0.4 0.3]; % 1xF parameters
    FC_freq = 1*F;
else
    error('Error: xor_voltagetrace currently only balances transmission for phase shifts = pi/2 or pi.')
end

for i = 1:length(is_inhibition_present)
    
    output = run_xor_model(input1, input2, Pmax_e(i), is_inhibition_present(i), FC_freq);
    spktrain1 = output.spktrain1;
    spktrain2 = output.spktrain2;
    inputs = output.inputs;
    spk_xor_expected_10ms = output.spk_xor_expected_10ms;
    spk_xor_expected_5ms = output.spk_xor_expected_5ms;
    
    figure;
    
    subplot(3,1,1)
    if i == 1
        plot(tvec,spktrain1(1,:), 'Color', 'r')
        title("Excitation only")
    else
        plot(tvec,spktrain1(1,:), 'Color', [0.4 0.4 1])
        title("Paired feed-forward EI")
    end
    yticks([])
    ylabel("Output 1, A1")
    set(gca, 'LineWidth', 1)
    box off;

    subplot(3,1,2)
    if i == 1
        plot(tvec,spktrain1(2,:), 'Color', 'r')
    else
        plot(tvec,spktrain1(2,:), 'Color', [0.4 0.4 1])
    end
    yticks([])
    ylabel("Output 2, A1")
    set(gca, 'LineWidth', 1)
    box off;

    subplot(3,1,3)
    if i == 1
        plot(tvec,spktrain2(1,:), 'Color', 'r')
    else
        plot(tvec,spktrain2(1,:), 'Color', [0.4 0.4 1])
    end
    yticks([])
    ylabel("Circuit output")
    set(gca, 'LineWidth', 1)
    box off;
    xlabel('Time (s)')
    
end

% Plot inputs 1 and 2
fig1 = figure;
set(fig1,'defaultAxesColorOrder',[[0 0 0]; [0.4 0.4 0]]);

subplot(2,1,1)
yyaxis left
plot(tvec, inputs(1,:), '-', 'Color', [0.8 0.8 0])
yticks([0 1])
ylabel({'Input 1 (spikes)'})
hold on;
yyaxis right
plot(tvec, signal1)
yticks([0 50 100])
ylabel({'Theoretical' 'firing rate (Hz)'})
title('Input 1')
set(gca, 'LineWidth', 1)
box off;

subplot(2,1,2)
yyaxis left
plot(tvec, inputs(2,:), '-', 'Color', [0.8 0.8 0])
yticks([0 1])
ylabel({'Input 2 (spikes)'})
hold on;
yyaxis right
plot(tvec, signal2)
yticks([0 50 100])
ylabel({'Theoretical' 'firing rate (Hz)'})
xlabel('Time (s)')
title('Input 2')
set(gca, 'LineWidth', 1)
box off;



% Plot expected output
figure;
subplot(2,1,1)
plot(tvec, spk_xor_expected_10ms, 'Color', [0.7 0.4 1])
yticks([0 1])
title("Expected XOR output with 10ms windows")
box off;

subplot(2,1,2)
plot(tvec, spk_xor_expected_5ms, 'Color', [0.7 0.4 1])
yticks([0 1])
title("Expected XOR output with 5ms windows")
xlabel("Time (s)")
box off;

end