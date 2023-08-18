%% Function definition
% Input: 
%       tSpan - Time span
%       freqVals - Frequency
%       expPeaks - Experimental values for ion concentrations
%       lb - Lower bound for estimation
%       ub - upper bound for estimation
%
% Output:
%       pSol - Particle swarm solution 

function pSol = SkelMuscleCa_paramEst(tSpan,freqVals,expPeaks,lb,ub,yinit)

    psOptions = optimoptions('particleswarm','UseParallel',false,'HybridFcn',@fmincon,...
                             'PlotFcn','pswplotbestf');%,'MaxStallIterations',15);
    
    numParam = length(lb);
    pSol = particleswarm(@pToObj_SS,numParam,lb,ub,psOptions);
    pToObj_SS(pSol)
    
    function objVal = pToObj_SS(pVec) % Obj function should be scalar
        lowATP = 0;
        [~,y_SS] = SkelMuscleCa1(tSpan,freqVals, lowATP, yinit, pVec);
       
        objVal = sqrt (mean ((y_SS(9) - expPeaks(1,:)).^2)) / mean(expPeaks(1,:)) + ...
                 sqrt (mean ((y_SS(2) - expPeaks(2,:)).^2)) / mean(expPeaks(2,:)) + ...
                 sqrt (mean ((y_SS(7) - expPeaks(3,:)).^2)) / mean(expPeaks(3,:)) + ...
                 sqrt (mean ((y_SS(14) - expPeaks(4,:)).^2))/ mean(expPeaks(4,:)) + ...
                 sqrt (mean ((y_SS(6) - expPeaks(5,:)).^2)) / mean(expPeaks(5,:)); %nRMSE
    end
end 



%rng default
% objVal = sum(((y_SS(9) - expPeaks(1,:)) / expPeaks(1,:)).^2) + ...
        %          sum(((y_SS(2) - expPeaks(2,:)) / expPeaks(2,:)).^2) + ...
        %          sum(((y_SS(7) - expPeaks(3,:)) / expPeaks(3,:)).^2) + ...
        %          sum(((y_SS(14) - expPeaks(4,:))/ expPeaks(4,:)).^2) + ...
        %          sum(((y_SS(6) - expPeaks(5,:)) / expPeaks(5,:)).^2);