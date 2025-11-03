close all
%pkg load signal
%pkg load communications

filename = '/mnt/c/Users/sha83592/Documents/octave/ADC3669EVM_Input70MHzCW-CHB_SampleRate256MHz_TwosComplement.csv';
data = dlmread(filename, ',', [0 0 inf 0]); % read only first column

max_adc_val = 32768; %16-bit two's complement input data
data_normalized = data/max_adc_val; %normalize magnitude of input data to ADC full-scale level

% === Simulation Parameters ===
fs = 256e6;       % Sampling freqyency
fc = 70e6;        % Input tone frequency
t = (1/fs)*(0:(length(data)-1)); % Time vector
df = fs/length(data); %FFT frequency resolution


% === Plot Time Domain Signal ===
figure;
plot(t, data);
title('Input Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0,10/fc]);
grid();

window_correction_factor = 2;
w = hann(length(data))*window_correction_factor;
data_normalized_windowed = w.*data_normalized;

% === Low-Pass Filtering (optional) ===
%[b, a] = butter(5, 2*fc/fs, 'low'); % Design a 5th order Butterworth low-pass filter
%data_normalized_windowed = filter(b, a, data_normalized_windowed);

% === Plot FFT ===
N = length(t);
f = (0:N/2-1) * (fs/N);

figure;
fft_input = fft(data_normalized_windowed)/N;
fft_input = fft_input(1:(N/2)+1); % single-sided spectrum
fft_input(2:end-1) = 2*fft_input(2:end-1);
plot(f, 20*log10(abs(fft_input)));
title('Input Signal - FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBFS)');
xlim([0, fs/2]);
ylim([-150,0]);
yticks(-150:10:0)
grid();

% Open a binary file for writing
%fid = fopen('/home/tsharp/octave/scripts/recovered_signal.bin', 'wb');

% Write the float data to the file
%fwrite(fid, interleaved_data, 'float32'); % Write as interleaved single-precision floats

%fclose(fid);

%fid2 = fopen('/home/tsharp/octave/scripts/demod_signal.bin', 'rb');
%fm_demod_signal = fread(fid2, Inf, 'float');

