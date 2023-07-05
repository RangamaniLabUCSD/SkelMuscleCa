%% start parameter optimization
lb = 0.5*ones(3,1);
ub = 2*ones(3,1);
tSpan = [0 1];
freqVals = 0;
expPeaks = 0.1;
pSol = particleSwarmSolve(tSpan, freqVals, expPeaks, lb, ub);

function pSol = particleSwarmSolve(tSpan,freqVals,expPeaks,lb,ub)

    psOptions = optimoptions('particleswarm','UseParallel',false,'HybridFcn',@fmincon,...
                             'PlotFcn','pswplotbestf');%,'MaxStallIterations',15);
    rng default
    numParam = length(lb);
    pSol = particleswarm(@pToObj_SS,numParam,lb,ub,psOptions);
    pToObj_SS(pSol)
    
    function objVal = pToObj(pVec)
        objVal = 0;
        for i = 1:length(freqVals)
            [~,y] = SkelMuscleCa_AK1(tSpan, freqVals(i), false, pVec);
            CaPeak = max(y(:,9));
            objVal = objVal + (CaPeak - expPeaks(i))^2;
        end
    end

    function objVal = pToObj_SS(pVec)
        y_SS = SkelMuscleCa_SS(false, pVec);
        objVal = (y_SS(9) - expPeaks).^2;
    end
end