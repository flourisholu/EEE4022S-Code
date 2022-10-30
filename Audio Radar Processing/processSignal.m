%% Author: Flourish Oluwamakinde
% This function proccesses the signal received from the receiver of the
% audio radar and produces all the Doppler frequencies from the
% received signal. Starts by removing delay, then applies a bandpass and
% notch filter; this is followed by baseband conversion before the signal
% is finally passed through a lowpass fiilter to remove the 2fc signals that appear.

%% Define constants and parameters

% Fs -> Sampling rate [Hz]. So Fs samples are obtained in one second
% Fc_Hz -> Carrier Frequency [Hz]
% Time_Duration_s -> Transmit time of signal [s]
% Tx_Signal -> Transmitted signal
% data -> Received signal
% data_bpf -> Received signal after filtered by BPF
% data_notch -> Received signal after notch filtering
% data_shifted -> Received signal after shifted to baseband
% data_out -> Received signal after lowpass filter; this is input to
% spectrogram

% only output -> data_out and t
function [data_out, t, N] = processSignal(Fs, Fc_Hz, TimeDuration_s, TxSignal, data0)
    %% Step 1: shift the received signal by delay [s]. 
    % This is because there is a delay after transmission, before the radar
    % starts receiving the echoed signal.

    %SpeedSoundWave_ms = 343; % [m/s] -> speed of sound wave c
    dt = 1/Fs; % Sampling period
    delay = 0.07; 
    start = floor(delay/(1/Fs));
    TxSignal = TxSignal(start:end-1);
    data = data0(start:end-1);
    t = delay:dt:TimeDuration_s; % time vector for pulse
    N = length(data); %number of samples

    % Plot transmitted and received signal without delay
    data_fft = fft(data);
    data_fftshift = fftshift(data_fft);
    f0 = (-N/2:N/2-1)*(Fs/N);
    y_received = (abs(data_fftshift).^2)/N;
    
    tsignal_fft = fft(TxSignal);
    tsignal_fftshift = fftshift(tsignal_fft);
    y_transmitted = (abs(tsignal_fftshift).^2)/N;

    subplot(8,1,1)
    plot(t,TxSignal) %transmitted signal in time-domain
    title('Transmitted signal in time-domain')
    xlabel('Time (s)')
    ylabel('Amplitude')

    subplot(8,1,2)
    plot(f0/1e3,y_transmitted)
    title('Fourier transform of transmitted signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')

    subplot(8,1,3)
    plot(t,data) %received signal in time-domain
    title('Received signal in time-domain')
    xlabel('Time (s)')
    ylabel('Amplitude')

    subplot(8,1,4)
    plot(f0/1e3,y_received)
    title('Fourier transform of received signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')
    %% Step 2: pass signal through bandpass filter

    % Remove DC offset and noise around DC with BPF - dont want it shifted when we
    % baseband shift.
    
    Fstop1 = 5000;
    Fstop2 = 15100;

    [b,a] = butter(6, [Fstop1 Fstop2]/(Fs/2),'bandpass');
    data_bpf = filter(b, a, data); %filter the received signal
    data_notch_fft = fft(data_bpf);
    data_notch_fftshift = fftshift(data_notch_fft);
    y_bpf= (abs(data_notch_fftshift).^2)/N;
    %fvtool(b,a); %visualise the bandpass filter

    % Plot the frequency response of the filtered signal
    subplot(8,1,5)
    plot(f0/1e3,y_bpf)
    title('Fourier transform of bandpass-filtered received signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')
    %% Step 3: apply a notch filter
    % Remove stationary clutter located at Fc_Hz
    % only going to sense things that are moving

    Q = 1000; % Q-factor of notch filter
    w0 = Fc_Hz/(Fs/2); % must be 0 < w0 < 1, where 1 = pi rad per sample in frequency range
    bw = w0/Q; % bandwidth at -3 dB point
    [b,a]=iirnotch(w0,bw); % view with fvtool(b,a)

    % Filter the received signal using the notch filter
    data_notch = data_bpf - mean(data_bpf);
    data_notch = filter(b,a, data_notch);
    data_notch_fft = fft(data_notch);
    data_notch_fftshift = fftshift(data_notch_fft);
    y_notch= (abs(data_notch_fftshift).^2)/N;
    %fvtool(b,a) %visualise the notch filter

    % Plot the frequency response of the notch-filtered signal
    subplot(8,1,6)
    plot(f0/1e3,y_notch)
    title('Fourier transform of notch-filtered received signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')
    %% Step 4: baseband conversion
    % Use Baseband function to shift received notch signal to baseband
    
    Ts = size(data_notch,1)/Fs;
    t = 0:dt:Ts-dt;
    N = length(data_notch);

    It = data_notch.*cos(2*pi*Fc_Hz*t)'; % real component
    Qt = data_notch.*sin(2*pi*Fc_Hz*t)'; % imaginary component
    data_baseband = complex(It,-Qt); % complex exponential

    %data_shifted = Baseband(data_notch, Fc_Hz, Fs);
    data_baseband_fft = fft(data_baseband);
    data_baseband_fftshift = fftshift(data_baseband_fft);
    y_baseband= (abs(data_baseband_fftshift).^2)/N;

    % Plot the frequency response of the baseband-converted signal
    subplot(8,1,7)
    plot(f0/1e3,y_baseband)
    title('Fourier transform of baseband-shifted received signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')
    %% Step 5: pass the signal through lowpass filter
    % Pass the baseband-shifted signal through a lowpass filter

    [M_lpf,P_lpf] = butter(6,Fc_Hz/(Fs/2),'low');
    data_out = filter(M_lpf,P_lpf,data_baseband);
    data_out_fft = fft(data_out);
    data_out_fftshift = fftshift(data_out_fft);
    y_lpf = (abs(data_out_fftshift).^2)/N;

    % Plot the frequency response of the lowpass-filtered signal
    subplot(8,1,8)
    plot(f0/1e3,y_lpf)
    title('Fourier transform of lowpass-filtered received signal')
    xlabel('Frequency (kHz)')
    ylabel('Magnitude')
end