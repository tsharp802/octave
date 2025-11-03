clc;
clear all;
close all;

% M-sequence generator setup
n = 4;  % Degree of the polynomial (can be adjusted)
c = [1 0 0 1 1];  % Feedback polynomial coefficients (example for n=4)
c = c(2:end);
M = 2^n - 1;  % Length of the M-sequence

% Initialize the shift register with a non-zero state
shift_register = [1 zeros(1, n-1)];
a = zeros(1, M);

% Generate M-sequence
for j = 1:M
    feedback_bit = mod(sum(shift_register .* c), 2);
    a(j) = shift_register(end);
    shift_register = [feedback_bit, shift_register(1:end-1)];
end

bp = 0.001;  % Bit period
disp('Binary information at Transmitter:');
disp(a);

% Representation of transmitting binary information as digital signal
bit = [];
for n = 1:length(a)
    if a(n) == 1
        se = ones(1, 100);
    else
        se = zeros(1, 100);
    end
    bit = [bit se];
end
t1 = bp/100:bp/100:100*length(a)*(bp/100);

% MSK Modulation
A = 1;  % Amplitude of carrier signal
br = 1/bp;  % Bit rate
f1 = br;  % Frequency for binary '1'
f2 = br/2;  % Frequency for binary '0'(ensuring frequency deviation iss br/2)
t2 = bp/99:bp/99:bp;
ss = length(t2);
m = [];
phase = 0;  % Initial phase to maintain phase continuity

for i = 1:length(a)
    if a(i) == 1
        freq = f1;
    else
        freq = f2;
    end
    y = A * cos(2 * pi * freq * t2 + phase);
    m = [m y];
    % Update phase for continuity
    phase = phase + 2 * pi * freq * bp;
end

t3 = bp/99:bp/99:bp*length(a);
figure;
hold on;
plot(t1, bit, 'LineWidth', 2, 'DisplayName', 'Digital Signal');% Plot digital signal
grid on;
% Plot MSK modulated signal
plot(t3, m, 'DisplayName', 'MSK Modulated Signal');
xlabel('Time (sec)');
ylabel('Amplitude (volt)');
title('Digital Signal and MSK Modulated Signal');
legend show;
axis([0 bp*length(a) -1.5 1.5]);
