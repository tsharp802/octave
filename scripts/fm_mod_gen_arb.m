close all
pkg load signal
pkg load communications
%graphics_toolkit("fltk");

% === Simulation Parameters ===
fs = 10.24e6;   % Sampling frequency
fc_if = 100e3;  % Carrier frequency
fc_rf = 1e6;    % RF carrier frequency
fm = 1e3;       % Message tone frequency
fdev = 5e3;     % FM frequency deviation
num_samples = 10*fs/fm;  % Capture one cycle of the message signal
T_total = num_samples/fs; % Total duration of the signal
t = 0:1/fs:T_total-(1/fs); % Time vector

% === 1. Generate Baseband I and Q Signals ===
amp_m = 1;
amp_i = 1;
amp_q = 1;
amp_q_imbalanced = 0.9;
dc_offset_i = 0.0;
dc_offset_q = 0.1;
phi_offset = pi/18;

% Baseband message signal
m = amp_m * (fdev/fm) * sin(2*pi*fm*t);

i_signal = amp_i * cos(m);
q_signal = amp_q * sin(m);
q_signal_imbalanced = amp_q_imbalanced * sin(m + phi_offset) + dc_offset_q;
bb_out = 1 * (i_signal + j*q_signal);
bb_out_q_imbalanced = 1 * (i_signal + j*q_signal_imbalanced);

% === 2. IQ Modulation (IF upconversion) ===
amp_if_carrier_i = 1;
amp_if_carrier_i_imbalanced = 0.9;
phi_offset_if_carrier_i = pi/18;
amp_if_carrier_q = 1;
amp_if_carrier_q_imbalanced = 0.9;
phi_offset_if_carrier_q = pi/18;
dc_offset_if_q = 0.1;

if_carrier = exp(j*2*pi*fc_if*t);
if_carrier_i = amp_if_carrier_i * real(if_carrier);
if_carrier_i_imbalanced = amp_if_carrier_i_imbalanced * real(if_carrier.*exp(j*phi_offset_if_carrier_i));
if_carrier_q = amp_if_carrier_q * imag(if_carrier);
if_carrier_q_imbalanced = amp_if_carrier_q_imbalanced * imag(if_carrier.*exp(j*phi_offset_if_carrier_q));

% Mix the signals
if_i = i_signal .* if_carrier_i;
if_i_alt = i_signal .* if_carrier_q;
if_q = q_signal .* if_carrier_q;
if_q_alt = q_signal .* if_carrier_i;
if_q_imbalanced = q_signal_imbalanced .* if_carrier_q;
if_q_carrier_imbalanced = q_signal .* if_carrier_q_imbalanced;
if_q_alt_imbalanced = q_signal_imbalanced .* if_carrier_i;
if_q_alt_carrier_imbalanced = q_signal .* if_carrier_i_imbalanced;

% Combine the mixed signals to produce the final IF signal.
if_signal_i = if_i - if_q;
if_signal_q = if_i_alt + if_q_alt;
if_signal_q_imbalanced = if_i_alt + if_q_alt_imbalanced;
if_signal_q_carrier_imbalanced = if_i_alt + if_q_alt_carrier_imbalanced + dc_offset_if_q;
if_out = 1 * (if_signal_i + j*if_signal_q);
if_out_q_imbalanced = 0.1 * (if_signal_i + j*if_signal_q_imbalanced);
if_out_q_carrier_imbalanced = 0.1 * (if_signal_i + j*if_signal_q_carrier_imbalanced);

% === 2. IQ Modulation (RF Upconversion) ===
amp_carrier_rf_i = 1;
amp_carrier_rf_q = 1;
amp_carrier_rf_q_imbalanced = 0.9;
phi_offset_carrier_rf_q = pi/18;

carrier_rf = exp(j*2*pi*fc_rf*t);
carrier_rf_i = amp_carrier_rf_i * real(carrier_rf);
carrier_rf_q = amp_carrier_rf_q * imag(carrier_rf);
carrier_rf_q_imbalanced = amp_carrier_rf_q_imbalanced * imag(carrier_rf.*exp(j*phi_offset_carrier_rf_q));

% Mix the signals
rf_i = if_signal_i .* carrier_rf_i;
rf_q = if_signal_q .* carrier_rf_q;
rf_signal = rf_i - rf_q;
rf_q_carrier_imbalanced = if_q .* carrier_rf_q_imbalanced;
rf_signal_q_carrier_imbalanced = rf_i - rf_q_carrier_imbalanced;

% Export waveform to .mat file for import into R&S ARB Toolbox
save -v6 fm_mod_1kHz_tone_5kHz_deviation_baseband.mat bb_out;
save -v6 fm_mod_1kHz_tone_5kHz_deviation_baseband_q_imbalanced.mat bb_out_q_imbalanced;
save -v6 fm_mod_1kHz_tone_5kHz_deviation_100kHz_if.mat if_out;
save -v6 fm_mod_1kHz_tone_5kHz_deviation_100kHz_if_q_imbalanced.mat if_out_q_imbalanced;
save -v6 fm_mod_1kHz_tone_5kHz_deviation_100kHz_if_q_if_carrier_imbalanced.mat if_out_q_carrier_imbalanced;

% === 4. IQ Demodulation (IF Downconversion) ===
% Recover the original I and Q signals from the RF signal.
% This involves another mixing process with the same carrier, followed by low-pass filtering.
demod_if_i_mixed = 2 * rf_signal .* carrier_rf_i;
demod_if_q_mixed = 2 * rf_signal .* -carrier_rf_q;

% === 5. Low-Pass Filtering ===
% Low-pass filters are necessary to remove the high-frequency components that result from the mixing process.
[b, a] = butter(5, 1.5*fc_rf/fs, 'low'); % Design a 5th order Butterworth low-pass filter
recovered_if_i_signal = filter(b, a, demod_if_i_mixed);
recovered_if_q_signal = filter(b, a, demod_if_q_mixed);
recovered_if_signal = recovered_if_i_signal - recovered_if_q_signal;

% === 4. IQ Demodulation (BB Downconversion) ===
% Recover the original I and Q signals from the RF signal.
% This involves another mixing process with the same carrier, followed by low-pass filtering.
demod_i_mixed_bb = 2 * recovered_if_i_signal .* if_carrier_i;
demod_q_mixed_bb = 2 * recovered_if_i_signal .* -if_carrier_q;

% === 5. Low-Pass Filtering ===
% Low-pass filters are necessary to remove the high-frequency components that result from the mixing process.
[b, a] = butter(5, 1.5*fc_if/fs, 'low'); % Design a 5th order Butterworth low-pass filter
recovered_i_signal = filter(b, a, demod_i_mixed_bb);
recovered_q_signal = filter(b, a, demod_q_mixed_bb);

% === 3. Plot Modulated Signal ===
figure;
subplot(4,1,1);
plot(t, i_signal, 'b', t, q_signal, 'r', t, q_signal_imbalanced, 'g');
title('Baseband I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I-Signal', 'Q-Signal', 'Q-Signal imbalanced');

subplot(4,1,2);
plot(t, if_signal_i, 'b', t, if_signal_q, 'r');
title('Modulated IF Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,3);
plot(t, rf_signal);
title('Modulated RF Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,1,4);
plot(t, recovered_i_signal, 'b', t, recovered_q_signal, 'r');
title('Recovered I and Q Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I-Signal', 'Q-Signal');

% Plot the spectrum
figure;
N = length(t);
f = (-N/2:N/2-1) * (fs/N);

% Spectrum of baseband signals
subplot(4,1,1);
fft_i = fftshift(fft(i_signal))/N;
fft_q = fftshift(fft(q_signal))/N;
plot(f, 20*log10(abs(fft_i)), 'b', f, 20*log10(abs(fft_q)), 'r');
title('Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10*fm,10*fm])
ylim([-100, 0]);

% Spectrum of IF signal
subplot(4,1,2);
fft_if = fftshift(fft(if_signal_i))/N;
fft_if_imbalanced = fftshift(fft(if_signal_q_imbalanced))/N;
plot(f, 20*log10(abs(fft_if)), 'b', f, 20*log10(abs(fft_if_imbalanced)), 'r');
title('IF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%xlim([fc_rf-10*fm, fc_rf+10*fm])
ylim([-100, 0]);

% Spectrum of RF signal
subplot(4,1,3);
fft_rf = fftshift(fft(rf_signal))/N;
fft_rf_imbalanced = fftshift(fft(rf_q_carrier_imbalanced))/N;
plot(f, 20*log10(abs(fft_rf)), 'b', f, 20*log10(abs(fft_rf_imbalanced)), 'r');
title('RF Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%xlim([fc_rf-10*fm, fc_rf+10*fm])
ylim([-100, 0]);

% Spectrum of recovered baseband signals
subplot(4,1,4);
fft_i_recovered = fftshift(fft(recovered_i_signal))/N;
fft_q_recovered = fftshift(fft(recovered_q_signal))/N;
fft_bb_recovered = fftshift(fft(recovered_i_signal - recovered_q_signal))/N;
plot(f, 20*log10(abs(fft_i_recovered)), 'b', f, 20*log10(abs(fft_q_recovered)), 'r', f, 20*log10(fft_bb_recovered), 'g');
title('Recovered Baseband Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10*fm,10*fm])
ylim([-100, 0]);

interleaved_data_if = zeros(1, 2 * length(t));
interleaved_data_if(1:2:end) = recovered_if_i_signal;
interleaved_data_if(2:2:end) = recovered_if_q_signal;

interleaved_data_bb = zeros(1, 2 * length(t));
interleaved_data_bb(1:2:end) = recovered_i_signal;
interleaved_data_bb(2:2:end) = recovered_q_signal;

% Open a binary file for writing
fid_if = fopen('/home/tsharp/octave/scripts/signal_if.bin', 'wb');
fid_bb = fopen('/home/tsharp/octave/scripts/signal_bb.bin', 'wb');

% Write the float data to the file
fwrite(fid_if, interleaved_data_if, 'float32'); % Write as interleaved single-precision floats
fwrite(fid_bb, interleaved_data_bb, 'float32'); % Write as interleaved single-precision floats
fclose(fid_if);
fclose(fid_bb);
