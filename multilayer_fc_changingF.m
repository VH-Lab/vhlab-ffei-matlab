% MULTILAYER_FC_CHANGINGF: Compares propagation of a signal that changes
% frequency over time, given different forms of transmission (FFE vs FFEI) 
% within a multilayer circuit
% 

clear;

L = 4; % number of layers in the circuit

% Time and input ST parameters
dt = 1e-4;
PR = 100;
t_start = 0;
t_end = 20;
tvec = t_start:dt:t_end;
Nt = length(tvec);

% Prepare changing frequency input: 
%   0s-2s: 5Hz
%   2s-4s: 10Hz
%   4s-10s: 10Hz->100Hz (stepping)
%   10s-18s: 100Hz->5Hz (stepping)
%   18s-20s: 5Hz
bin_size = 1; % size of FC measurement bins
binvec = 0:dt:1;
N_bins = t_end/bin_size; % number of FC measurement bins
F = zeros(1,length(tvec));

for t = 1:length(tvec)
    freq_mult = ceil(tvec(t)); % during frequency stepping phases, step frequency based on current time, rounded up
    if tvec(t) <= 2
        F(t) = 5;
    elseif tvec(t) <= 4
        F(t) = 10;
    elseif tvec(t) <= 10
        F(t) = 10 + 15 * (freq_mult-4);
    elseif tvec(t) <= 18
        F(t) = 100 - 95/8 * (freq_mult-10);
    else
        F(t) = 5;
    end
end

% Number of trials
n_trials = 50;

% Conductance parameters
Pmax_e = 1.6976e-7;
tau1i = [0.020 0]; % 20ms for ffei, 0ms (placeholder) for ffe
Pmax_factor = [0.18 0.19]; % scaled Pmax so that power at 5Hz = 75Hz


% Prepare sinusoidal rate of changing frequency - this will be used to
% generate input spike trains during each trial
signal = max(0,PR*sin(2*pi*F.*tvec));

% Prepare noisy input parameters + store noisy inputs in noisy_input
noisy_input = zeros(L,Nt);
T_background = pi/PR;

% Prepare vectors to store fourier spectrum
    % Spectrum is calculated in 20 bins (dimension 1), and encompasses
    % length(binvec) frequencies
FC_ffei = zeros(N_bins, length(binvec)); 
FC_ffe = zeros(N_bins, length(binvec)); 


for n = 1:n_trials

    disp(['Starting trial ' num2str(n) '...'])
        
    % Set up vectors to store inputs to each layer
    % 1: ffei, 2: ffe
    inputs_1 = zeros(L+1,Nt); 
    inputs_2 = zeros(L+1,Nt);

    % Provide different initial inputs to each circuit
    inputs_1(1,:) = poisson_spiketrain(dt, signal, t_end, 1);
    inputs_2(1,:) = poisson_spiketrain(dt, signal, t_end, 1);
   
    
    for i = 1:L
        
        inst_input_1 = generate_input(5, 1, 1, 'tmax', t_end); % generate (noisy) input to next layer (we assign the input to layer i to the input signal)
        inst_input_2 = generate_input(5, 1, 1, 'tmax', t_end); % note: F = 5 Hz as a placeholder for the signal frequency
        
        inst_input_1.signal = inputs_1(i,:);
        inst_input_2.signal = inputs_2(i,:);
        
        output_data_1 = run_triad_model(inst_input_1, 1, 'Pmax_e', Pmax_e * Pmax_factor(1), 'tau1i', tau1i(1), 'noise_level', 1);
        output_data_2 = run_triad_model(inst_input_2, 1, 'Pmax_e', Pmax_e * Pmax_factor(2), 'tau1i', tau1i(2), 'noise_level', 1);
        
        inputs_1(i+1,:) = output_data_1.post_spktrain;
        inputs_2(i+1,:) = output_data_2.post_spktrain;
        
    end
    
    % Calculate fourier spectrum in 1 second bins, store in 3d matrices
    % FC_ffei and FC_ffe
    
    output_1 = inputs_1(L+1, :)/dt; % convert output to rate
    output_2 = inputs_2(L+1, :)/dt; % convert output to rate
    
    for bin = bin_size:N_bins
        
        % Calculate FC_pct spectrum for each bin
        FC_F_ffei = abs(fouriercoeffs(output_1((bin-bin_size)/dt+1:bin/dt+1),dt));
        FC_avg_ffei = mean(FC_F_ffei);
        FC_ffei(bin,:,n) = FC_F_ffei/FC_avg_ffei;
        
        FC_F_ffe = abs(fouriercoeffs(output_2((bin-bin_size)/dt+1:bin/dt+1),dt));
        FC_avg_ffe = mean(FC_F_ffe);
        FC_ffe(bin,:,n) = FC_F_ffe/FC_avg_ffe;
    end
    
end


FC_ffei_avg = mean(FC_ffei, 3);
FC_ffe_avg = mean(FC_ffe, 3);

FC_ffei_avg_100 = FC_ffei_avg(:,1:100);
FC_ffe_avg_100 = FC_ffe_avg(:,1:100);

FC_ffei_avg_200 = FC_ffei_avg(:,1:200);
FC_ffe_avg_200 = FC_ffe_avg(:,1:200);

figure;
plot(tvec, F, 'k')
xlabel('Time (s)')
ylabel('Input modulation frequency (F)')
box off;


figure;
s1 = imagesc(1:20,1:200,FC_ffei_avg_200');
colormap(gray(256));
set(gca, 'YDir', 'normal')
title('FC spectrum: FFEI')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar;
caxis([0,10]);
box off;

figure;
s2 = imagesc(1:20,1:200,FC_ffe_avg_200');
colormap(gray(256));
set(gca, 'YDir', 'normal')
title('FC spectrum: FFE')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar;
caxis([0,10]);
box off;