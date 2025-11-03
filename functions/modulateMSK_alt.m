% The code is developed by SalimWireless.com
% Function to perform MSK modulation
function [signal_i, signal_q, timeVec] = modulateMSK_alt(bits, carrierFreq, baudRate, sampleFreq)
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

    % Define time parameters
    numBits = length(bits);
    symbDur = 1 / baudRate;
    timeVec = 0:1/sampleFreq:numBits * symbDur - (1/sampleFreq);

    % Compute phase shifts
    phaseShift = zeros(1, numBits);
    for idx = 1:numBits
        %phaseShift(idx) = mod(phaseShift(idx-1) + ((pi * idx) / 2) * ((diffEncBits(idx-1) - diffEncBits(idx))), 2 * pi);
        phaseShift(idx) = (pi/2)*sum(diffEncBits(1:idx));
    end

    % Generate MSK waveform
    symbolIdx = floor(timeVec / symbDur) + 1;
    signal_i = cos(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift(symbolIdx));
    signal_q = sin(2 * pi * (carrierFreq + diffEncBits(symbolIdx) / (4 * symbDur)) .* timeVec + phaseShift(symbolIdx));
end
