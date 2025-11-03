clc;
clear;
close all;
pkg load signal;
pkg load image;

% Functions
function h = GMSK_gaussian_filter1(T, sps)
    span = 2;
    a = T*0.1;
    BT = sqrt(log(2)/2)./(a);
    t = linspace(-span*T/2,span*T/2,sps)';
    h = (sqrt(pi)/a)*exp(-(pi*t./a).^2);
    %h = BT * sqrt((2*pi) / log(2)) .* exp(-(((2 * pi^2) * (BT^2)) .* t.^2) / log(2));
    h = (h / sum(h));
end

function h = GMSK_matched_filter(T, sps)
    t = (-1.5*T:T/sps:1.5*T); BT = 0.5;
    h = BT * sqrt((2*pi) / log(2)) .* exp(-(((2 * pi^2) * (BT^2)) .* t.^2) / log(2));
    h = (pi / (2 * sum(h))) * h / sqrt(sum(h));
end

function downsampled_output = GMSK_downsample(start_idx, end_idx, sps, input_signal)
    downsampled_output = input_signal(start_idx:sps:end-end_idx);
end

function quantized_signal = GMSK_ADC(input_signal)
    quantized_signal = sign(input_signal);
end

% MSK waveform generation
function [signal_i, signal_q, timeVec, bitSeq, phaseSeq, symbolIdx] = modulateMSK(bits, carrierFreq, baudRate, sampleFreq)
    % Converts a binary bit sequence into an MSK-modulated signal
    % Inputs:
    %   bits        - Binary input sequence
    %   carrierFreq - Carrier frequency
    %   baudRate    - Symbol rate
    %   sampleFreq  - Sampling frequency
    % Outputs:
    %   signal      - Modulated MSK signal
    %   timeVec     - Corresponding time vector

    % Convert bits to NRZ format (-1, 1)
    diffEncBits = 2 * bits - 1;
    %diffEncBits = bits;
    diffEncBits = [diffEncBits(1), diffEncBits]; % Append initial value
    encBits = [bits(1), bits];

    % Define time parameters
    numBits = length(bits);
    symbDur = 1 / baudRate;
    timeVec = 0:1/sampleFreq:numBits * symbDur - (1/sampleFreq);

    % Compute phase shifts
    phaseShift_i = zeros(1, numBits + 1);
    phaseShift_q = zeros(1, numBits + 1);
    for idx = 2:numBits+1
        %phaseShift_i(idx) = mod(phaseShift_i(idx-1) + ((pi * idx) / 2) * ((diffEncBits(idx-1) - diffEncBits(idx))), 2 * pi);
        %phaseShift_q(idx) = mod(phaseShift_q(idx-1) + ((pi * idx) / 2) * ((diffEncBits(idx-1) - diffEncBits(idx))), 2 * pi);
        %phaseShift_q(idx) = mod(phaseShift_q(idx-1) + ((pi * idx) / 2) * ((diffEncBits(idx-1) - diffEncBits(idx))) + pi * (encBits(idx-1) - encBits(idx)), 2 * pi);
        %phaseShift(idx) = phaseShift(idx-1) + (pi/2) * diffEncBits(idx-1);
    end
    phaseShift_i = phaseShift_i(2:end);
    phaseShift_q = phaseShift_q(2:end);
    diffEncBits = diffEncBits(2:end);
    encBits = encBits(2:end);
    phaseSeq = phaseShift_i;
    bitSeq = encBits;

    % Generate MSK waveform
    symbolIdx = floor(timeVec / symbDur) + 1;
    signal_i = cos(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift_i(symbolIdx));
    signal_q = sin(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift_q(symbolIdx));
    %signal_i = zeros(1,length(timeVec));
    %signal_q = zeros(1,length(timeVec));
    %for i = 1:length(timeVec)
      %signal_i = cos(2 * pi * (carrierFreq + diffEncBits_upsampled(i) / (4 * symbDur)) * timeVec(i) + phaseShift_upsampled(i));
      %signal_q = sin(2 * pi * (carrierFreq + diffEncBits_upsampled(i) / (4 * symbDur)) * timeVec(i) + phaseShift_upsampled(i));
    %endfor
endfunction


%Input parameters
fc_m = 0e6;
fc = 120e6;%10e6;
fm = 16e6;
fs = 256e6;
samples_per_bit = fs/fm;
bit_duration = 1/fm;
record_length = 32768;
num_bytes = record_length / (8 * samples_per_bit);

% Define a bit sequence
bitSeq = repmat([0,0,0,0,1,1,1,1],1,num_bytes);
bitSeq_upsampled = kron(bitSeq, ones(1, samples_per_bit));
% Apply Gaussian filter
gauss_filter = GMSK_gaussian_filter1(bit_duration, samples_per_bit);
bitSeq_filt = conv(gauss_filter, bitSeq_upsampled);
%bitSeq_filt = imgaussfilt(bitSeq_upsampled);

% Perform MSK modulation
[modSignal_i, modSignal_q, t, bitSeq, phaseSeq, symbolIdx] = modulateMSK(bitSeq, fc_m, fm, fs);

%t = 0:(1/fs):(length(bitSeq) * (1/fm) - (1/fs));
%t = 0:(1/fs):(1/fs)*(length(bitseq)-1);



carrier_i = cos(2*pi*fc*t);
carrier_q = sin(2*pi*fc*t);

temp_i = (modSignal_i + ones(1,length(modSignal_i)))/2;
temp_q = (modSignal_q + ones(1,length(modSignal_q)))/2;
signal_i = round(temp_i.'* 65535);
signal_q = round(temp_q.'* 65535);

temp_i_rf = (modSignal_i - (ones(1,length(modSignal_i))/2));
temp_q_rf = (modSignal_q - (ones(1,length(modSignal_q))/2));
rf = modSignal_i.*carrier_i - modSignal_q.*carrier_q;

figure;
subplot(3,1,1);
plot(t, modSignal_i, 'r');
subplot(3,1,2);
plot(t, modSignal_q, 'b');
subplot(3,1,3);
plot(t, rf, 'g');

figure;
N = length(t);
f_ss = (0:N/2-1) * (fs/N);
f = (-N/2:N/2-1) * (fs/N);
window_correction_factor = 2;
w = hann(length(t)).' * window_correction_factor;
modSignal_i_windowed = w.*modSignal_i;
modSignal_q_windowed = w.*modSignal_q;
rf_windowed = w.*rf;

fft_input_i = fftshift(fft(modSignal_i_windowed))/N;
fft_input_q = fftshift(fft(modSignal_q_windowed)/N);
fft_input_rf = fft(rf)/N;
fft_input_rf = fft_input_rf(1:(N/2)); % single-sided spectrum
fft_input_rf(2:end-1) = 2*fft_input_rf(2:end-1);
subplot(3,1,1);
plot(f/1e6, 20*log10(abs(fft_input_i)));
xlim([-fs/2,fs/2])
xticks(linspace(-fs/2,fs/2,fs/2/8)/1e6)
ylim([-60 0]);
grid on
subplot(3,1,2);
plot(f/1e6, 20*log10(abs(fft_input_q)));
xlim([-fs/2,fs/2])
xticks((-fs/2:fs/2/8:fs/2)/1e6)
ylim([-60 0]);
grid on
subplot(3,1,3);
plot(f_ss/1e6, 20*log10(abs(fft_input_rf)));
xlim([0,fs/2])
xticks((0:fs/2/8:fs/2)/1e6)
ylim([-60 0]);
grid on
title('Input Signal - FFT');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dBFS)');
%xlim([0, fs/2]);

figure;
% Plot the modulated signal
subplot(3,1,1);
samples = 1:numel(bitSeq);
plot(t, bitSeq(symbolIdx));
title('Original message signal');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-2,2]);

% Plot the modulated signal
subplot(3,1,2);
plot(t, modSignal_i, 'r', t,modSignal_q,'b');
title('MSK baseband signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the modulated signal
subplot(3,1,3);
plot(t, rf);
title('MSK modulated signal');
xlabel('Time (s)');
ylabel('Amplitude');





