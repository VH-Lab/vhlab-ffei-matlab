% MULTILAYER_FC: Generates figures comparing final output power of 
% 4-layer triad synapse circuit models (Figure 5).

% Define number of layers in circuit
L = 4;
Pmax_e = 1.6976e-9;

% 1) Compare paired-ffei, ex, and no transmission given a fixed Pmax
F_vec = logspace(log10(5), log10(1000), 50);
Pmax_factor = [0.18 0.19 0]; % scaled Pmax so that power at 5Hz = 75Hz - originally 0.275 for pffei
n_trials = 10;

FC_vec = zeros(3, length(F_vec));
FCavg_vec = zeros(3, length(F_vec));
FCpct_vec = zeros(3, length(F_vec));

% NOTE on order of inputs:
% compute_multilayer_fc takes 1)frequency, 2) no. layers,
% 3) factor of Pmax relative to 1.6976e-7 S, and 4) whether or not
% paired-ffei is the mode of transmission

disp("Making comparisons between noise only, FFE, and FFEI circuits...")

for i = 1:length(F_vec)
    
    if mod(i-1,5) == 0
    disp(['Progress: ' num2str(i-1) ' frequencies analyzed'])
    end
    
    FC_ffei = zeros(1, n_trials);
    FC_ffe = zeros(1, n_trials);
    FC_noise = zeros(1, n_trials);
    
    FCavg_ffei = zeros(1, n_trials);
    FCavg_ffe = zeros(1, n_trials);
    FCavg_noise = zeros(1, n_trials);
    
    FCpct_ffei = zeros(1, n_trials);
    FCpct_ffe = zeros(1, n_trials);
    FCpct_noise = zeros(1, n_trials);
    
    for j = 1:n_trials
        [FC_ffei(j), FCavg_ffei(j), FCpct_ffei(j)] = compute_multilayer_fc(F_vec(i), L, Pmax_factor(1), 1);
        [FC_ffe(j), FCavg_ffe(j), FCpct_ffe(j)] = compute_multilayer_fc(F_vec(i), L, Pmax_factor(2), 0);
        [FC_noise(j), FCavg_noise(j), FCpct_noise(j)] = compute_multilayer_fc(F_vec(i), L, Pmax_factor(3), 0);
    end
    
   FC_vec(1,i) = mean(FC_ffei);
   FC_vec(2,i) = mean(FC_ffe);
   FC_vec(3,i) = mean(FC_noise);
   
   FCavg_vec(1,i) = mean(FCavg_ffei);
   FCavg_vec(2,i) = mean(FCavg_ffe);
   FCavg_vec(3,i) = mean(FCavg_noise);
   
   FCpct_vec(1,i) = mean(FCpct_ffei);
   FCpct_vec(2,i) = mean(FCpct_ffe);
   FCpct_vec(3,i) = mean(FCpct_noise);
end


figure;
plot(F_vec, FC_vec(1,:), 'b')
hold on;
plot(F_vec, FC_vec(2,:), 'r')
plot(F_vec, FC_vec(3,:), 'Color', [0 0.8 0])
xlabel('Input modulation frequency (Hz)')
ylabel('FC magnitude of final layer output at the frequency F (Hz) (FC_F)')
legend('PFFEI', 'E only', 'noise only')
title({'Avg. power of 4-layer circuit output at different input modulation' ...
    ['frequencies: synaptic conductance = ' num2str(Pmax_e) ' S']})
set(gca, 'XScale', 'log')
box off;

figure;
plot(F_vec, FCavg_vec(1,:), 'b')
hold on;
plot(F_vec, FCavg_vec(2,:), 'r')
plot(F_vec, FCavg_vec(3,:), 'Color', [0 0.8 0])
xlabel('Input modulation frequency (Hz)')
ylabel({'Avg. FC magnitude over entire spectrum  of final', 'layer output (Hz) (FC_a_v_g)'})
legend('PFFEI', 'E only', 'noise only')
title({'Avg. overall power of 4-layer circuit output for inputs with different' ...
    ['modulation frequencies: synaptic conductance = ' num2str(Pmax_e) ' S']})
set(gca, 'XScale', 'log')
box off;

figure;
plot(F_vec, FCpct_vec(1,:), 'b')
hold on;
plot(F_vec, FCpct_vec(2,:), 'r')
plot(F_vec, FCpct_vec(3,:), 'Color', [0 0.8 0])
plot([5 1e3],[1 1],'--', 'Color', [0.5 0.5 0.5])
xlabel('Input modulation frequency (Hz)')
ylabel('FC_F/FC_a_v_g')
legend('PFFEI', 'E only', 'noise only')
title({'(Avg. power at input modulation frequency)/(Avg. power overall) of' ...
    ['4-layer circuit output: synaptic conductance = ' num2str(Pmax_e) ' S']})
set(gca, 'XScale', 'log')
box off;

% Compare ffe vs ffei over different output power levels

Pmax_factor_2_ffei = [1/2 0.83 1 1.115 1.165]; % FFEI
Pmax_factor_2_ffe = [1/2 3/4 1 1.45 1.7]; % E only

FC_vec_ffei = zeros(length(Pmax_factor_2_ffei), length(F_vec));
FC_vec_ffe = zeros(length(Pmax_factor_2_ffei), length(F_vec));

disp("Making comparisons over output power...")

for i = 1:length(F_vec)
    
    if mod(i-1,5) == 0
    disp(['Progress: ' num2str(i-1) ' frequencies analyzed'])
    end
    
    for p = 1:length(Pmax_factor_2_ffei)
        FC_ffei = zeros(1, n_trials);
        FC_ffe = zeros(1, n_trials);
        
        for j = 1:n_trials
            
            FC_ffei(j) = compute_multilayer_fc(F_vec(i), L, Pmax_factor_2_ffei(p)*Pmax_factor(1), 1);
            FC_ffe(j) = compute_multilayer_fc(F_vec(i), L, Pmax_factor_2_ffe(p)*Pmax_factor(2), 0); 
        end
        
    FC_vec_ffei(p,i) = mean(FC_ffei);
    FC_vec_ffe(p,i) = mean(FC_ffe);
    end
end

cmap = parula(5);

figure;
hold on;
for i = 5:-1:2
    plot(F_vec,FC_vec_ffe(i,:), '--', 'color', cmap(6-i,:));
    plot(F_vec,FC_vec_ffei(i,:), '-', 'color', cmap(6-i,:));
end
plot(F_vec, FC_vec_ffe(1, :), 'r--')
plot(F_vec, FC_vec_ffei(1,:), 'r-')

xlabel('Input frequency F (Hz)')
ylabel('Avg. power of final layer output at the frequency F (Hz)')
title({'Avg. power of 4-layer circuit output at different input modulation' ...
    ['frequencies: synaptic conductance scaled from ' num2str(Pmax_e) ' S']})
set(gca, 'XScale', 'log')
box off;