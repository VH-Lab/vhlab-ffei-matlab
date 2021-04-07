function [modelinput] = generate_input(F, n_input, is_noise_present, varargin)
% Default: no current injection, 100Hz peak firing rate, 5 second trial
% with no start or end buffer

% Input characteristics
PR = 100;
dt = 1e-4; % time bins
tmax = 5; % signal time
tbuffer_start = 0;
tbuffer_end = 0; % "buffer" time after signal

% define noise parameters (if any)
noise_rate = PR/pi;
noise_N = 50; % number of noisy inputs

% define parameters for rectified sinusoidal current injection (if any)
Im_amp = 0;

% define parameters for a specified spiking input e.g. spike train of
% constant uniform rate, custom signal or noise (if any)
manual_signal = 0;

% Assign values to any parameters specified in function call:
assign(varargin{:});
% typically: PR, F, tmax, dt, Pmax_e, tau1i, delay, alpha, Im_F, Im_amp...

% Time parameters that are adjusted after manual parameter assignment:
t_total = tbuffer_start+tmax+tbuffer_end; % total trial length
tsigvec = 0:dt:tmax; % time vector for signal
tvec = -tbuffer_start:dt:t_total-tbuffer_start; % time vector for trial
Nt = length(tvec); % number of time intervals
Nbuffer_start = ceil(tbuffer_start/dt); % number of time intervals occupied by start buffer

% If sinusoidal direct current injection is present (Im_amp is
% assigned in function call), then nonzero Im is generated:
Im = zeros(1,length(tvec));
Im(Nbuffer_start+1:length(tsigvec)+Nbuffer_start) = max(Im_amp * sin(2*pi*F*tsigvec), 0);

if manual_signal == 0 % if no manual spiking signal is provided ->
    
    spktimes = spiketrain_sinusoidal(PR, F, 0, 0, 0, tmax, dt);
    signal = spiketimes2bins(spktimes, tvec);
    
    if n_input>1 % if signal is comprised of n>1 inputs ->
        for i = 2:n_input
            % Generate a Poisson spike train with input characteristics
            spktimes = spiketrain_sinusoidal(PR, F, 0, 0, 0, tmax, dt);
            signal = signal + spiketimes2bins(spktimes, tvec);
        end
    end
else 
    if length(manual_signal) == Nt % check that manually-inputted signal is the correct length
       signal = manual_signal;
    else
        error(['Error: Manually-inputted signal is ' num2str(length(manual_signal)) ' time units in length, while trial is ' num2str(Nt) ' time units in length.'])
    end
end


if is_noise_present == 1
    % Generate 50 noisy inputs
    noise = poisson_spiketrain(dt, noise_rate, t_total, 1);
    for j=2:noise_N
        noise = noise + poisson_spiketrain(dt, noise_rate, t_total, 1);
    end   
    
elseif is_noise_present == 0
    
    noise = zeros(1, Nt);
    
else
    error("Error: Invalid input for is_noise_present (boolean).")
end

modelinput = var2struct('F', 'dt', 'tvec', 'signal', 'noise', 'Im');

end