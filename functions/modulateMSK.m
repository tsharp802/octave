% The code is developed by SalimWireless.com
% Function to perform MSK modulation
function [signal_i, signal_q, timeVec] = modulateMSK(bits, carrierFreq, baudRate, sampleFreq)
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
    diffEncBits = [-1, diffEncBits]; % Append initial value

    % Define time parameters
    numBits = length(bits);
    symbDur = 1 / baudRate;
    timeVec = 0:1/sampleFreq:numBits * symbDur - (1/sampleFreq);

    % Compute phase shifts
    phaseShift = zeros(1, numBits + 1);
    for idx = 2:numBits+1
        %phaseShift(idx) = mod(phaseShift(idx-1) + ((pi * idx) / 2) * ((diffEncBits(idx-1) - diffEncBits(idx))), 2 * pi);
        phaseShift(idx) = phaseShift(idx-1) + (pi/2) * diffEncBits(idx-1);
    end
    phaseShift = phaseShift(2:end);
    diffEncBits = diffEncBits(2:end);

    % Generate MSK waveform
    symbolIdx = floor(timeVec / symbDur) + 1;
    signal_i = cos(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift(symbolIdx));
    signal_q = sin(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift(symbolIdx));
end
