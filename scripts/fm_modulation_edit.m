close all
pkg load signal
pkg load communications

% === Simulation Parameters ===
fs = 1e6;       % Sampling frequency (1 MHz)
fc1 = 0;     % Carrier frequency (100 kHz)
fc2 = 100e3;
fm = 1e3;    % Message frequency for I component (1 kHz)
T_total = 10/fm; % Total duration of the signal (10 ms)
t = 0:1/fs:T_total-(1/fs); % Time vector

% === 1. Generate Baseband I and Q Signals ===
% The in-phase (I) and quadrature (Q) signals are the messages to be transmitted.
% In a real system, these would be the data streams.
%i_signal = cos(2*pi*fm_i*t);
%q_signal = sin(2*pi*fm_q*t);

%s = exp(1j*(fdev/fm)*sin(2*pi*fm*t) + 1j*2*pi*fc1*t);
amp_sig_fund = 1;
amp_sig_harm = 0;
fdev = 5000;   %FM frequency deviation
sig_harm_index = 10;
m = amp_sig_fund*(fdev/fm)*sin(2*pi*fm*t);
s = exp(1j*m + 1j*amp_sig_harm*(fdev/(sig_harm_index*fm))*sin(2*pi*sig_harm_index*fm*t)+ 1j*2*pi*fc1*t);

i_signal = real(s);
q_signal = imag(s);

phi_harm_i = (pi/2);
phi_harm_q = (0);
amp_harm_i = 0.5;
amp_harm_q = 0.5;

% === 2. IQ Modulation (Upconversion) ===
% Mix the I and Q signals with the carrier frequency.
% An ideal IQ mixer multiplies the I-component with a cosine and the Q-component with a sine.
carrier_i = cos(2*pi*fc2*t);
carrier_i_harm = cos(2*pi*fc2*t) + amp_harm_i * cos(2*pi*3*fc2*t + phi_harm_i);
carrier_q = sin(2*pi*fc2*t);
carrier_q_harm = sin(2*pi*fc2*t) + amp_harm_q * sin(2*pi*3*fc2*t + phi_harm_q);

% Mix the signals
rf_i = i_signal .* carrier_i;
rf_i_harm = (i_signal + 1*i_signal.*i_signal.*i_signal) .* carrier_i_harm;
rf_q = q_signal .* carrier_q;
rf_q_harm = (q_signal + 1*q_signal.*q_signal.*q_signal) .* carrier_q_harm;

% Combine the mixed signals to produce the final RF signal.
rf_signal = rf_i + rf_q;
rf_signal_harm = rf_i_harm + rf_q_harm;

% === 3. Plot Modulated Signal ===
figure;
subplot(4,1,1);
plot(t, i_signal, 'b', t, q_signal, 'r');
title('Baseband I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I-Signal', 'Q-Signal');

subplot(4,1,2);
plot(t, rf_signal);
title('Modulated RF Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,3);
plot(t, rf_signal_harm);
title('Modulated RF Signal with LO Harmonics');
xlabel('Time (s)');
ylabel('Amplitude');

% === 4. IQ Demodulation (Downconversion) ===
% Recover the original I and Q signals from the RF signal.
% This involves another mixing process with the same carrier, followed by low-pass filtering.
demod_i_mixed = rf_signal_harm .* carrier_i_harm;
demod_q_mixed = rf_signal_harm .* carrier_q_harm;

% === 5. Low-Pass Filtering ===
% Low-pass filters are necessary to remove the high-frequency components that result from the mixing process.
[b, a] = butter(5, 2*fc2/fs, 'low'); % Design a 5th order Butterworth low-pass filter
recovered_i_signal = filter(b, a, demod_i_mixed);
recovered_q_signal = filter(b, a, demod_q_mixed);

% === 6. Plot Recovered Signals ===
subplot(4,1,4);
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
subplot(4,1,1);
fft_i = fftshift(fft(i_signal));
fft_q = fftshift(fft(q_signal));
plot(f, 20*log10(abs(fft_i)), 'b', f, 20*log10(abs(fft_q)), 'r');
title('Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc1-10*fm, fc1+10*fm]);
ylim([0,80]);

% Spectrum of RF signal
subplot(4,1,2);
fft_rf = fftshift(fft(rf_signal));
plot(f, 20*log10(abs(fft_rf)));
title('RF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc2-fc1-10*fm, fc2-fc1+10*fm]);
ylim([0,80]);

% Spectrum of RF signal
subplot(4,1,3);
fft_rf_harm = fftshift(fft(rf_signal_harm));
plot(f, 20*log10(abs(fft_rf_harm)));
title('RF Signal Spectrum with LO Harmonics');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc2-fc1-10*fm, fc2-fc1+10*fm]);
ylim([0,80]);

subplot(4,1,4);
fft_i_recovered = fftshift(fft(recovered_i_signal));
fft_q_recovered = fftshift(fft(recovered_q_signal));
plot(f, 20*log10(abs(fft_i_recovered)), 'b', f, 20*log10(abs(fft_q_recovered)), 'r');
title('Recovered Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fc1-10*fm, fc1+10*fm]);
ylim([0,80]);

interleaved_data = zeros(1, 2 * length(i_signal));
interleaved_data(1:2:end) = i_signal;
interleaved_data(2:2:end) = q_signal;

% Open a binary file for writing
fid = fopen('recovered_signal.bin', 'wb');

% Write the float data to the file
fwrite(fid, interleaved_data, 'float32'); % Write as interleaved single-precision floats

fclose(fid);

fid2 = fopen('demod_signal.bin', 'rb');
fm_demod_signal = fread(fid2, Inf, 'float');

figure;
plot(t, postpad(fm_demod_signal,10000));
title('Demodulated FM signal');
xlabel('Time (s)');
ylabel('Amplitude');
