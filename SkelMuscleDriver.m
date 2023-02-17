%% Test different frequencies of stimulus
figure
hold on
tSpan = [0 2];
freqVec = logspace(-1,2,10);
for i = 1:length(freqVec)
    [t,y] = SkelMuscleCa(tSpan, freqVec(i), false);
    plot(t,y(:,8))
end