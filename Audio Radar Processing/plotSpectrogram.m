%% Author: Flourish Oluwamakinde
% This script performs the processing needed to compute the spectrogram and
% determine the speed of the moving vehicles. A signal of Fc_Hz is
% transmitted
%% Clear variables and command window

clear all;
close all;
clc;

%% Define constants and parameters
load('/MATLAB Drive/EEE4022S Final Year Project/Flourish/Audio Radar Processing/data/movingHand.mat')
[data_out, t, N] = processSignal(Fs, Fc_Hz, TimeDuration_s, TxSignal, RX_signal);

% Compute the spectrogram 
samplesPerFrame =  2^(12);           % Ensure its a power of 2
overlapFactor = 0.9;                    % Overlap factor of 50% between successive frames

[S, f, t_s] = customSTFT(data_out, Fs, samplesPerFrame, overlapFactor);

datan = abs(S);
datan = datan-min(datan(:));
datan = datan/max(datan(:));
datan = 20*log10(datan);

for i = 2:length(f)-1
    for j = 2:length(t_s)-1
       if datan(i,j) < -60
           datan(i,j) = -100;
       end
    end
end

datan = smoothdata(datan, 'movmedian', 50);

% Calculate velocities (max, avg, instantaneous)
max_frequencies = zeros(length(t_s),1);
for frame = 1:length(t_s)
    [max_S, index] = max(datan(:, frame)); %gets max FFT magnitude in frame
    if (max_S > -60)
        max_frequencies(frame) = f(index); %store corresponding frequency value
    else
        max_frequencies(frame) = 0;
    end
end

v_factor = SpeedSoundWave_ms/(2*Fc_Hz);
max_speed = max(abs(max_frequencies))*v_factor; %maximum speed
disp('Maximum speed: ');
disp(max_speed);
avg_speed = mean(abs(max_frequencies))*v_factor; %average speed
disp('Average speed: ');
disp(avg_speed);
instant_v = max_frequencies*v_factor;

%instant_v = smoothdata(instant_v, 'movmean', 5);

% Plot the speeds
figure
plot(t_s, instant_v);
yline(max_speed, "r");
yline(avg_speed, "b--");
title("Maximum, average, and instantaneous speeds");
legend('Instantaneous speed', 'Maximum speed', 'Average speed')

% Plot the spectrogram 
v = f.*v_factor;
figure; imagesc(t_s, v, datan, [-60, 0]);
ylim([0 10])
title("Spectrogram of received signal")
xlabel('Time (s)');
ylabel('Doppler velocity (m/s)');
grid on;
colorbar;
colormap('jet');
set(gca,'YDir','normal')