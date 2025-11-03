close all
pkg load signal
pkg load communications

% === Simulation Parameters ===
fs = 10e6;       % Sampling frequency (1 MHz)
fc1 = 10e3;     % Carrier frequency (100 kHz)
fc2 = 1000e3;
fm = 1e3;    % Message frequency for I component (1 kHz)
fm_q = 1e3;     % Message frequency for Q component (1 kHz)
beta = 5;
fdev = beta * fm;   %FM frequency deviation
T_total = 10/fm; % Total duration of the signal (10 ms)
t = 0:1/fs:T_total-(1/fs); % Time vector

% === 1. Generate Baseband I and Q Signals ===
% The in-phase (I) and quadrature (Q) signals are the messages to be transmitted.
% In a real system, these would be the data streams.
%i_signal = cos(2*pi*fm_i*t);
%q_signal = sin(2*pi*fm_q*t);

s = exp(1j*(fdev/fm)*sin(2*pi*fm*t) + 1j*2*pi*fc1*t);

i_signal = real(s);
q_signal = imag(s);

% set LO harmonic amplitude and phase
phi_harm_i = (pi/2);
phi_harm_q = (0);
amp_harm_i = 1;
amp_harm_q = 0.1;

% === 2. IQ Modulation (Upconversion) ===
% Mix the I and Q signals with the carrier frequency.
% An ideal IQ mixer multiplies the I-component with a cosine and the Q-component with a sine.
carrier_i = cos(2*pi*fc2*t);
carrier_i_harm = sin(2 * pi * fc2 * t) + amp_harm_i * sin(2 * pi * 3 * fc2 * t + phi_harm_i);
carrier_q = sin(2*pi*fc2*t);
carrier_q_harm = sin(2 * pi * fc2 * t + (pi/2)) + amp_harm_q * sin(2 * pi * 3 * fc2 * t + phi_harm_q);

% Mix the signals
rf_i = i_signal .* carrier_i_harm;
rf_q = q_signal .* carrier_q_harm;

% Combine the mixed signals to produce the final RF signal.
rf_signal = rf_i + rf_q;

% === 3. Plot Modulated Signal ===
figure;
subplot(3,1,1);
plot(t, i_signal, 'b', t, q_signal, 'r');
title('Baseband I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I-Signal', 'Q-Signal');

subplot(3,1,2);
plot(t, rf_signal);
title('Modulated RF Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% === 4. IQ Demodulation (Downconversion) ===
% Recover the original I and Q signals from the RF signal.
% This involves another mixing process with the same carrier, followed by low-pass filtering.
demod_i_mixed = rf_signal .* carrier_i;
demod_q_mixed = rf_signal .* carrier_q;

% === 5. Low-Pass Filtering ===
% Low-pass filters are necessary to remove the high-frequency components that result from the mixing process.
[b, a] = butter(5, 2*fc2/fs, 'low'); % Design a 5th order Butterworth low-pass filter
recovered_i_signal = filter(b, a, demod_i_mixed);
recovered_q_signal = filter(b, a, demod_q_mixed);

% === 6. Plot Recovered Signals ===
subplot(3,1,3);
plot(t, recovered_i_signal, 'b', t, recovered_q_signal, 'r');
title('Recovered I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Recovered I', 'Recovered Q');

% Optional: Plot the spectrum to visualize the upconversion and downconversion
figure;
N = length(rf_signal);
f = (-N/2:N/2-1) * (fs/N);

% Spectrum of baseband signals
subplot(2,1,1);
fft_i = fftshift(fft(i_signal));
fft_q = fftshift(fft(q_signal));
plot(f, abs(fft_i), 'b', f, abs(fft_q), 'r');
title('Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc1-10*fm, fc1+10*fm]);

% Spectrum of RF signal
subplot(2,1,2);
fft_rf = fftshift(fft(rf_signal));
plot(f, abs(fft_rf));
title('RF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc2+fc1-10*fm, fc2+fc1+10*fm]);
