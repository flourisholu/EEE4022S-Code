%% Author: Flourish Oluwamakinde
% This is a custom algorithm to compute the STFT of the signal received from the ultrasonic radar.

%% References
% Short Time Fourier Transform (STFM) Function written by EEE4022 students in 2020 to compute a spectrogram

% Authors: Nadir Aboobaker and Zahier Parker 
% Reviewed: Dr Abdul Gaffar
% Year: 2020/2021
% Adapted by Ian Edwards 2021
%% Define constants and parameters

% y -> Processed received signal
% Fs  -> Sampling rate [Hz]
% N -> Number of samples in a single frame. Should be a power of 2 
% OverlapFactor -> Overlap between successive frames (decimal value). Example: 0.5 means there is 50% overlap of samples between successive  frames
% S -> Computed spectrogram (matrix, complex numbers)
% f -> Frequency vector axis [Hz]
% t -> Time vector axis [s]

function [S, f, t] = customSTFT(y, Fs, N, OverlapFactor)
    zero_pad = 10*N;
    y = y(:); %convert the signal to a column-vector
    ylength = size(y,1);
    window = hamming(N); %use a window function to get short-time of signal
    overlapSamples = floor(OverlapFactor*N);

    % Get the number of frames to be taken
    hop = N - overlapSamples;
    frames = 1 + floor((ylength-N)/hop);

    % Compute FFT, then apply window to get the STFT
    S = zeros(zero_pad/2+1, frames); % only store positive frequencies in the STFT matrix
    for counter = 0:(frames-1)
        y_window = y(1 + counter*hop : N + counter*hop).*window; %apply window to section of sampled data
        y_fft = fft(y_window, zero_pad);
        y_fft = y_fft(1:zero_pad/2+1)*2;
        S(:,1+counter) = y_fft; %add windowed section to STFT matrix
    end
    
    % Determine frequency and time vectors
    t = (0:hop:(frames-1)*hop)/Fs; 
    f = (0:1:(zero_pad/2))*Fs/zero_pad;
end 