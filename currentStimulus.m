function I = currentStimulus(t,pulsewidth)
% Parameters
period = 2.5; % Total period of one cycle (2.5 seconds)
pulseDuration = 0.5; % Duration of the stimulus in each cycle (0.5 seconds)
freq = 50; % Frequency of the stimulus (50 Hz)
amplitude = 1; % Amplitude of the current

% Calculate the number of full periods elapsed
numPeriods = floor(t / period);
timeInCurrentPeriod = t - numPeriods * period;

% Determine if within the pulse duration and calculate the stimulus phase
if timeInCurrentPeriod <= pulseDuration
    %phase = mod(freq * timeInCurrentPeriod, 1);
    if (mod(t,1/freq) < pulsewidth) % phase < 0.5
        I = amplitude;
    else
        I = 0;
    end
else
    I = 0; % No current outside the pulse duration
end
end
