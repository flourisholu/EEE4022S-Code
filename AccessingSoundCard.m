% Audio Radar
    

%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters

SpeedSoundWave_ms = 343;             % [m/s]  ->speed of wave
Fc_Hz = 10e3;                        % [Hz]
TimeDuration_s = 3;                  % [s]  
Fs = 44.1e3;                         % [Hz]   ->Sampling rate. So 44 100 samples are obtained in one second
Ts = 1/Fs;                           % Sampling period
t = 0:Ts:(TimeDuration_s);           % time vector for pulse

%% Generate the transmit signal

% Pure sinusoid 
TxSignal = sin(2*pi*Fc_Hz*t); % save this to workspace for easy testing


%% Play out transmit signal through the speakers

soundsc(TxSignal,Fs, 24) % Transmit the signal

%% Record received samples from the microphone

RecLength_samples = length(TxSignal);
RecLength_s = RecLength_samples*1/Fs; 
recObj = audiorecorder(Fs,24,1);
recordblocking(recObj, RecLength_s);  % Records audio for a fixed number of seconds
RX_signal = getaudiodata(recObj);   % Store recorded audio signal in double-precision array

figure; axes('fontsize', 12);
subplot(2,1,1);
plot(t, TxSignal); % plot transmit signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Transmit signal', 'fontsize', 12);
grid on;

subplot(2,1,2);
plot(t, RX_signal); % plot received signal
xlabel('Time (s)', 'fontsize', 12);
ylabel('Amplitude (linear)', 'fontsize', 12);
title('Received signal', 'fontsize', 12);
grid on;


