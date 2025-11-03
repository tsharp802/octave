close all
pkg load signal
pkg load communications
%graphics_toolkit("fltk");

% === Simulation Parameters ===
fs = 10.24e6;       % Sampling frequency
fc_if = 100e3;     % Carrier frequency
fc_rf = 1e6
fm = 1e3;    % Tone frequency
fdev = 5e3;   % FM frequency deviation
num_samples = 10*fs/fm;
T_total = num_samples/fs; % Total duration of the signal
t = 0:1/fs:T_total-(1/fs); % Time vector

% === 1. Generate Baseband I and Q Signals ===
amp_sig_fund = 1;
amp_i = 1;
amp_q = 1;
amp_q_distorted = 0.9;
dc_offset_i = 0.0;
dc_offset_q = 0.0;
phi_offset = 0;

m = amp_sig_fund*(fdev/fm)*sin(2*pi*fm*t);

i_signal = amp_i * cos(m);
q_signal = amp_q * sin(m);
q_signal_distorted = amp_q_distorted * sin(m + phi_offset) + dc_offset_q;

amp_carrier_i = 1;
amp_carrier_q = 1;
amp_carrier_q_distorted = 0.9;
phi_offset_carrier_q = pi/18;

% === 2. IQ Modulation (Upconversion) ===
carrier = exp(j*2*pi*fc_if*t);
carrier_i = amp_carrier_i * real(carrier);
carrier_q = amp_carrier_q * imag(carrier);
carrier_q_distorted = amp_carrier_q_distorted * imag(carrier.*exp(j*phi_offset_carrier_q));

% Mix the signals
if_i = i_signal .* carrier_i;
if_i_alt = i_signal .* carrier_q;
if_q = q_signal .* carrier_q;
if_q_alt = q_signal .* carrier_i;
if_q_distorted = q_signal_distorted .* carrier_q;
if_q_carrier_distorted = q_signal .* carrier_q_distorted;

% Combine the mixed signals to produce the final RF signal.
if_signal_i = if_i - if_q;
if_signal_q = if_i_alt + if_q_alt;
if_out = if_signali - j*if_signal_q;
if_signal_q_distorted = if_i - if_q_distorted;
if_signal_q_carrier_distorted = if_i - if_q_carrier_distorted;

amp_carrier_rf_i = 1;
amp_carrier_rf_q = 1;
amp_carrier_rf_q_distorted = 0.9;
phi_offset_carrier_rf_q = pi/18;

% === 2. IQ Modulation (2nd stage Upconversion) ===
carrier_rf = exp(j*2*pi*fc_rf*t);
carrier_rf_i = amp_carrier_rf_i * real(carrier_rf);
carrier_rf_q = amp_carrier_rf_q * imag(carrier_rf);
carrier_rf_q_distorted = amp_carrier_rf_q_distorted * imag(carrier_rf.*exp(j*phi_offset_carrier_rf_q));

% Mix the signals
rf_i = if_signal_i .* carrier_rf_i;
rf_q = if_signal_q .* carrier_rf_q;
rf_signal = rf_i - rf_q;
%rf_q_distorted = if_q_distorted .* carrier_q;
rf_q_carrier_distorted = if_q .* carrier_rf_q_distorted;

% Export waveform to .mat file for import into R&S ARB Toolbox
save -v6 fm_mod_1kHz_tone_5kHz_deviation_1MHz_carrier.mat if_out;

% Plot the spectrum
figure;
N = length(t);
f = (-N/2:N/2-1) * (fs/N);

% Spectrum of baseband signals
subplot(3,1,1);
fft_i = fftshift(fft(i_signal))/N;
fft_q = fftshift(fft(q_signal))/N;
plot(f, abs(fft_i), 'b', f, abs(fft_q), 'r');
title('Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10*fm,10*fm])

% Spectrum of IF signal
subplot(3,1,2);
fft_if = fftshift(fft(if_signal))/N;
fft_if_distorted = fftshift(fft(if_q_distorted))/N;
plot(f, abs(fft_if), 'b', f, abs(fft_if_distorted), 'r');
title('IF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%xlim([fc_rf-10*fm, fc_rf+10*fm])

% Spectrum of RF signal
subplot(3,1,3);
fft_rf = fftshift(fft(rf_signal))/N;
fft_rf_distorted = fftshift(fft(rf_q_carrier_distorted))/N;
plot(f, 20*log10(abs(fft_rf)), 'b', f, 20*log10(abs(fft_rf_distorted)), 'r');
title('RF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%xlim([fc_rf-10*fm, fc_rf+10*fm])
ylim([-100, 0]);

% === 4. IQ Demodulation (Downconversion) ===
% Recover the original I and Q signals from the RF signal.
% This involves another mixing process with the same carrier, followed by low-pass filtering.
demod_i_mixed = 2 * rf_q_carrier_distorted .* carrier_i;
demod_q_mixed = 2 * rf_q_carrier_distorted .* -carrier_q;

% === 5. Low-Pass Filtering ===
% Low-pass filters are necessary to remove the high-frequency components that result from the mixing process.
[b, a] = butter(5, 2*fc/fs, 'low'); % Design a 5th order Butterworth low-pass filter
recovered_i_signal = filter(b, a, demod_i_mixed);
recovered_q_signal = filter(b, a, demod_q_mixed);

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

subplot(3,1,3);
plot(t, recovered_i_signal, 'b', t, recovered_q_signal, 'r');
title('Recovered I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I-Signal', 'Q-Signal');

interleaved_data = zeros(1, 2 * length(t));
interleaved_data(1:2:end) = recovered_i_signal;
interleaved_data(2:2:end) = recovered_q_signal;

% Open a binary file for writing
fid = fopen('/home/tshar/octave/scripts/signal.bin', 'wb');

% Write the float data to the file
fwrite(fid, interleaved_data, 'float32'); % Write as interleaved single-precision floats
fclose(fid);
