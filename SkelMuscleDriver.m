%% Test different frequencies of stimulus
figure
hold on
tSpan = [0 1e5];
freqVec = 0;%logspace(0,2,10);
maxCa = zeros(size(freqVec));
for i = 1:length(freqVec)
    [t,y,SSParam] = SkelMuscleCa_AK1(tSpan, freqVec(i), false);
    plot(t,y(:,9))
    maxCa(i) = max(y(:,9));
    CaStored{i} = [t, y(:,9)]; 
end


%% Plot maximum ca vs. frequency
figure
semilogx(freqVec, maxCa, '-s')