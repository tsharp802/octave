clear;
close all;
pkg load signal;

fc = 1e6; %tone frequency
fs = 100e6; %sample rate
num_samples = 100;
num_bits = 16;
T = num_samples/fs;
t = 0:(1/fs):T-(1/fs);

% unsigned / offset binary format
signal_i = cos(2 * pi * fc * t);
signal_q = sin(2 * pi * fc * t);
signal_i_out = round((signal_i + ones(1,length(t)))/2*(2^num_bits-1)).';
signal_q_out = round((signal_q + ones(1,length(t)))/2*(2^num_bits-1)).';
save -v7 signal_i.mat signal_i;
save -v7 signal_q.mat signal_q;
data = [1, 0; 0, 1];
interleaved_data = zeros(length(t), 2);
interleaved_data(:,1) = signal_i;
interleaved_data(:,2) = signal_q;
interleaved_data_2 = zeros(1, 2 * length(t));
interleaved_data_2(1:2:end) = signal_i;
interleaved_data_2(2:2:end) = signal_q;
Ch1_CFrequency_Hz = 0.1000;
Ch1_ChannelName = 'IQValues';
Ch1_Clock_Hz = 1.0000e+08;
Ch1_Data = interleaved_data;
Ch1_Samples = num_samples;
Comment = 'tone';
DataType = 'float32';
Format = 'complex';
Name = 'Tone';
NumberOfChannels = 1;
save -v6 signal.mat Ch1_CFrequency_Hz Ch1_ChannelName Ch1_Clock_Hz Ch1_Data Ch1_Samples Comment DataType Format Name NumberOfChannels;
save -v6 test.mat data;

figure;
% Plot the modulated signal
plot(t, signal_i, 'b', t, signal_q, 'r');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
N = length(t);
f_ss = (0:N/2-1) * (fs/N);
window_correction_factor = 2;
w = hann(length(t)).' * window_correction_factor;
signal_i_windowed = w.*signal_i;
signal_q_windowed = w.*signal_q;
fft_i = fft(signal_i_windowed)/N;
fft_i = fft_i(1:(N/2)); % single-sided spectrum
fft_i(2:end-1) = 2*fft_i(2:end-1);
fft_q = fft(signal_q_windowed)/N;
fft_q = fft_q(1:(N/2)); % single-sided spectrum
fft_q(2:end-1) = 2*fft_q(2:end-1);
subplot(2,1,1);
plot(f_ss, 20*log10(abs(fft_i)));
title("In-phase signal FFT");
ylim([-60 0]);
grid on
subplot(2,1,2);
plot(f_ss, 20*log10(abs(fft_q)));
title("Quadrature signal FFT");
ylim([-60 0]);
grid on
