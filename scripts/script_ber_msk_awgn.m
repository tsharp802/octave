
clear all
N = 5*10^5; % number of bits or symbols

fsHz = 1; % sampling period
T    = 4; % symbol duration

Eb_N0_dB = [0:10]; % multiple Eb/N0 values
ct = cos(pi*[-T:N*T-1]/(2*T));
st = sin(pi*[-T:N*T-1]/(2*T));

for ii = 1:length(Eb_N0_dB)

    % MSK Transmitter
    ipBit = rand(1,N)>0.5; % generating 0,1 with equal probability
    ipMod =  2*ipBit - 1; % BPSK modulation 0 -> -1, 1 -> 1

    ai = kron(ipMod(1:2:end),ones(1,2*T));  % even bits
    aq = kron(ipMod(2:2:end),ones(1,2*T));  % odd bits

    ai = [ai zeros(1,T)  ]; % padding with zero to make the matrix dimension match
    aq = [zeros(1,T) aq ];  % adding delay of T for Q-arm

    % MSK transmit waveform
    xt = 1/sqrt(T)*[ai.*ct + j*aq.*st];

    % Additive White Gaussian Noise
    nt = 1/sqrt(2)*[randn(1,N*T+T) + j*randn(1,N*T+T)]; % white gaussian noise, 0dB variance

    % Noise addition
    yt = xt + 10^(-Eb_N0_dB(ii)/20)*nt; % additive white gaussian noise

    %% MSK receiver
    % multiplying with cosine and sine waveforms
    xE = conv(real(yt).*ct,ones(1,2*T));
    xO = conv(imag(yt).*st,ones(1,2*T));

    bHat = zeros(1,N);
    bHat(1:2:end) = xE(2*T+1:2*T:end-2*T) > 0 ; % even bits
    bHat(2:2:end) = xO(3*T+1:2*T:end-T) > 0 ;  % odd bits

    % counting the errors
    nErr(ii) = size(find([ipBit - bHat]),2);

end

simBer = nErr/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 10 10^-5 0.5])
grid on
legend('theory - bpsk', 'simulation - msk');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for MSK modulation');

