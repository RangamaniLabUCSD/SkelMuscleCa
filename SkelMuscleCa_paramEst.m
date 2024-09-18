%% Function for estimating parameters using Particle swarm estimation 
% Input:
%       yinit - Default state variable values
%       expPeaks - Experimental values for ion concentrations
%       lb - Lower bound for estimation
%       ub - upper bound for estimation
%       p_est - Default input parameters
%       Createplot - Logic command to output plots. 1 - Creates plot. 0 -
%       Omits plots
%
% Output:
%       pSol - Particle swarm solution
%       fval - Objective function value
%       exitflag - Reason for ending estimation (Values between -5 and 1.
%                  Explanation of corresponding exit condition avialabel at
%                  https://www.mathworks.com/help/gads/particleswarm.html)

function [pSol,fval,exitflag] = SkelMuscleCa_paramEst(lb,ub,yinit,p,Createplot)

psOptions = optimoptions('particleswarm','UseParallel',false,'HybridFcn',@fmincon,'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 50, 'SwarmSize', 2); %set swarmsize to 30 for TSCC and 50 to stall iter

% pSol Results
% load PSO_22-Aug-2024.mat pSol
param = importdata('InputParam1.xlsx');
p0 =  param.data;
% To check default parameter behavior
numParam = length(p0);
pVec = ones(1,numParam) ; 
pCur = pVec;
pCur = pCur(:) .* p0;
SkelMuscleObj(pCur);


delete(gcp('nocreate'))
% parpool(3) %% **CHANGE SWARMSIZE!** and save file location 
[pSol,fval,exitflag] = particleswarm(@SkelMuscleObj,length(lb),lb,ub,psOptions);
% [pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pCur = pVec;
highSensIdx = [1,5,15,16,19,22,24,30,32,34,35,43,74,83,91,92];
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
pCur = pCur(:) .* p0;
SkelMuscleObj(pCur)
filename = "PSO_" + date +".mat";
save(fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename'));
% save(fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename));

end

