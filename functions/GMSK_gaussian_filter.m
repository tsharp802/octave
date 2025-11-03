clear all;
close all;

% Functions
function h = GMSK_gaussian_filter1(T, sps)
    t = (-1.5*T:T/sps:1.5*T); BT = 0.3;
    h = BT * sqrt((2*pi) / log(2)) .* exp(-(((2 * pi^2) * (BT^2)) .* t.^2) / log(2));
    h = (pi / (2 * sum(h))) * h / sqrt(sum(h));
end

function h = GMSK_matched_filter(T, sps)
    t = (-1.5*T:T/sps:1.5*T); BT = 0.75;
    h = BT * sqrt((2*pi) / log(2)) .* exp(-(((2 * pi^2) * (BT^2)) .* t.^2) / log(2));
    h = (pi / (2 * sum(h))) * h / sqrt(sum(h));
end

function downsampled_output = GMSK_downsample(start_idx, end_idx, sps, input_signal)
    downsampled_output = input_signal(start_idx:sps:end-end_idx);
end

function quantized_signal = GMSK_ADC(input_signal)
    quantized_signal = sign(input_signal);
end

% Parameters
fc_m = 0;
fc = 12e6;
fm = 16e6;
fs = 256e6;
samples_per_bit = fs/fm;
bit_duration = 1/fm;
num_bits = 256;
sample_interval = bit_duration / samples_per_bit;
time_vector = 0:sample_interval:(num_bits * bit_duration);
time_vector(end) = [];

% Generate and modulate binary data
binary_data = repmat([0, 1],1,num_bits/2);
modulated_bits = 2 * binary_data - 1;
upsampled_signal = kron(modulated_bits, ones(1, samples_per_bit));
figure;
plot(time_vector, upsampled_signal);
title('Message Signal');

% Apply Gaussian filter
%filtered_signal = conv(GMSK_gaussian_filter1(bit_duration, samples_per_bit), upsampled_signal);
%filtered_signal = [filtered_signal, filtered_signal(end)];
filtered_signal = upsampled_signal;
figure; plot(filtered_signal);
title('Filtered Signal');

% Integration & GMSK modulation
integrated_signal = cumsum(filtered_signal);
gmsk_signal = exp(1i * integrated_signal);
signal_i = round((real(gmsk_signal) + ones(1,length(gmsk_signal)))/2*(2**16-1)).';
signal_q = round((imag(gmsk_signal) + ones(1,length(gmsk_signal)))/2*(2**16-1)).';

% Plotting the real and imaginary parts of the GMSK signal with labels
figure;
plot(real(gmsk_signal), 'b'); % Plot real part in blue
hold on;
plot(imag(gmsk_signal), 'r'); % Plot imaginary part in red
title('GMSK Modulated Signal');
xlabel('Samples');
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part'); % Adding labels to the legend


% Noiseless demodulation & matched filtering
matched_filter = GMSK_matched_filter(bit_duration, 7);
filt_signal = conv(matched_filter, gmsk_signal);
filt_signal = [filt_signal, filt_signal(end)];

% Extract phase, differentiate & downsample
phase_derivative = [unwrap(angle(filt_signal(1))), diff(unwrap(angle(filt_signal)))];
downsampled_signal = GMSK_downsample(70, 71, samples_per_bit, phase_derivative);
digital_output = GMSK_ADC(downsampled_signal);

% Plot demodulated signal
rect_pulses = repelem(digital_output, samples_per_bit);
time_axis = 0:1/samples_per_bit:length(digital_output);
figure;
plot(time_axis(1:end-1), rect_pulses);
title('Demodulated Signal');





