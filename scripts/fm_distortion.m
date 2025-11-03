close all
pkg load signal
pkg load communications

% Define parameters
fs = 100e+6; % Sampling frequency (Hz)
f_signal = 1000; % Frequency of the input signal (Hz)
f_signal_2 = 2000;
f_lo = 1000000; % Frequency of the Local Oscillator (LO) (Hz)
f_dev = 5000;

t = 0:1/fs:5/f_signal; % Time vector (1 second duration)

% Generate input signal
signal_in_i = cos(2 * pi * f_signal * t) + cos(2 * pi * f_signal_2 * t);
signal_in_q = sin(2 * pi * f_signal * t) + sin(2 * pi * f_signal_2 * t);

fm_signal_i = fmmod(signal_in_i, 0, fs, f_dev);
% apply 90deg phase shift to FM signal
%Y = abs(fft(fm_signal_i));
%N_y = length(Y);
%f_axis_y = (-N_y/2:N_y/2-1) * (fs/N_y);
%figure;
%plot(f_axis_y, Y);
%Y_shift = Y.*exp(1i*(pi/2));
%fm_signal_q = ifft(Y_shift);
fm_signal_q = fmmod(signal_in_q, 0, fs, f_dev)

% Generate Local Oscillator (LO) signal
phi_harm_i = (pi/2);
phi_harm_q = (pi/2);
amp_harm = 0.5;
lo_signal_i = sin(2 * pi * f_lo * t);
lo_signal_i_harm = sin(2 * pi * f_lo * t) + amp_harm * sin(2 * pi * 3 * f_lo * t + phi_harm_i);
lo_signal_q = sin(2 * pi * f_lo * t + (pi/2));
lo_signal_q_harm = sin(2 * pi * f_lo * t + (pi/2)) + amp_harm * sin(2 * pi * 3 * f_lo * t + phi_harm_q);

% Perform mixing (multiplication in the time domain)
%mixed_signal_i = fm_signal_i .* lo_signal_i;
%mixed_signal_q = fm_signal_q .* lo_signal_q;
mixed_signal_i = signal_in_i .* lo_signal_i;
mixed_signal_q = signal_in_q .* lo_signal_q;
mixed_signal = mixed_signal_i + mixed_signal_q;
%mixed_signal_i_harm = fm_signal_i .* lo_signal_i_harm;
%mixed_signal_q_harm = fm_signal_q .* lo_signal_q_harm;
mixed_signal_i_harm = signal_in_i .* lo_signal_i_harm;
mixed_signal_q_harm = signal_in_q .* lo_signal_q_harm;
mixed_signal_harm = mixed_signal_i_harm + mixed_signal_q_harm;

% Plot the signals
figure;
subplot(2,1,1);
plot(t, signal_in_i);
title('Input Signal (I)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, signal_in_q);
title('Input Signal (Q)');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
subplot(2,1,1);
plot(t, fm_signal_i);
title('FM Signal (I)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, fm_signal_q);
title('FM Signal (Q)');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
subplot(2,1,1);
plot(t, lo_signal_i);
title('Local Oscillator Signal (I)' );
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, lo_signal_q);
title('Local Oscillator Signal (Q)');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
subplot(2,1,1);
plot(t, lo_signal_i_harm);
title('Local Oscillator Signal (I) + Harmonics' );
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, lo_signal_q_harm);
title('Local Oscillator Signal (Q) + Harmonics');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
subplot(2,1,1);
plot(t, mixed_signal);
title('Mixed Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, mixed_signal_harm);
title('Mixed Signal + LO Harmonics');
xlabel('Time (s)');
ylabel('Amplitude');

% Analyze the frequency content of the mixed signal using FFT
N = length(mixed_signal);
f_axis = (-N/2:N/2-1) * (fs/N); % Frequency axis for FFT
%fft_mixed = fftshift(abs(fft(mixed_signal)));
fft_mixed = fft(mixed_signal);

N_harm = length(mixed_signal_harm);
f_axis_harm = (-N_harm/2:N_harm/2-1) * (fs/N_harm); % Frequency axis for FFT
%fft_mixed_harm = fftshift(abs(fft(mixed_signal_harm)));
fft_mixed_harm = ft(mixed_signal_harm);

figure;
plot(f_axis, fft_mixed);
%xlim([f_lo - 10e+3 f_lo + 10e+3])
title('Frequency Spectrum of Mixed Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

figure;
plot(f_axis, fft_mixed_harm);
%xlim([f_lo - 10e+3 f_lo + 10e+3])
title('Frequency Spectrum of Mixed Signal + LO Harmonics');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
