clear all;
close all;
pkg load signal;

fc = 1e6; %tone frequency
fs = 150e6; % max sample rate of SMBV100B
num_samples = fs/fc; %Capture one cycle of tone
T = num_samples/fs;
t = 0:(1/fs):T-(1/fs);

% unsigned / offset binary format
signal_i = cos(2 * pi * fc * t);
signal_q = sin(2 * pi * fc * t);
signal_out = signal_i + j*signal_q;
save -v7 tone_1MHz.mat signal_out;

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
grid on;
subplot(2,1,2);
plot(f_ss, 20*log10(abs(fft_q)));
title("Quadrature signal FFT");
ylim([-60 0]);
grid on;
