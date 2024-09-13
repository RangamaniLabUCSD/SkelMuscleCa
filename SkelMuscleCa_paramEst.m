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

psOptions = optimoptions('particleswarm','UseParallel',false,'HybridFcn',@fmincon,'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 20, 'SwarmSize', 2); %set swarmsize to 30 for TSCC and 50 to stall iter

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
% parpool(30) %% **CHANGE SWARMSIZE!** and save file location 
[pSol,fval,exitflag] = particleswarm(@SkelMuscleObj,length(lb),lb,ub,psOptions);
% [pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pCur = pVec;
highSensIdx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95];
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
pCur = pCur(:) .* p0;
SkelMuscleObj(pCur)
filename = "PSO_" + date +".mat";
save(fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename'));
% save(fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename));

end

