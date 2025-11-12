function mel_band_frequencies()
    % Total number of Mel bands
    numBands = 40;
    % Frequency limits
    lowFreq = 0;
    highFreq = 11000;

    % Convert frequency limits to Mel scale
    lowMel = hzToMel(lowFreq);
    highMel = hzToMel(highFreq);

    % Generate Mel points for each band
    melPoints = linspace(lowMel, highMel, numBands + 1);

    % Convert Mel points back to frequency
    freqPoints = melToHz(melPoints);

    % Display the frequency bands
    for i = 1:length(freqPoints)-1
        fprintf('Band %d: %f Hz to %f Hz\n', i, freqPoints(i), freqPoints(i+1));
    end
end

function mel = hzToMel(hz)
    % Convert a frequency in Hertz to Mel scale
    mel = 2595 * log10(1 + hz / 700);
end

function hz = melToHz(mel)
    % Convert a frequency in Mel scale back to Hertz
    hz = 700 * (10.^(mel / 2595) - 1);
end
