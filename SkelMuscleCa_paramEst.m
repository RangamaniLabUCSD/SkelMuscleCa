%% Function for estimating parameters using Particle swarm estimation 
% Input:
%       lb - Lower bound for estimation
%       ub - upper bound for estimation
%
% Output:
%       pSol - Particle swarm solution
%       fval - Objective function value
%       exitflag - Reason for ending estimation (Values between -5 and 1.
%                  Explanation of corresponding exit condition avialabel at
%                  https://www.mathworks.com/help/gads/particleswarm.html)

function [pSol,fval,exitflag] = SkelMuscleCa_paramEst(lb,ub)

 %set swarmsize to 30 for TSCC and 50 to stall iter
psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 10, 'SwarmSize', 15);

% pSol Results
% load PSO_18-Dec-2024.mat pSol
param = importdata('InputParamPsoResults.xlsx');
p0 =  param.data;
% To check default parameter behavior
numParam = length(p0); 
pVec = ones(1,numParam) ; 
pCur = pVec ;
pCur = pCur(:) .* p0;
SkelMuscleObj2(pCur)

delete(gcp('nocreate'))
if psOptions.UseParallel
    parpool(3) %% **CHANGE SWARMSIZE!** and save file location
end
saveProgress = true;
fAnon = @(x)SkelMuscleObj2(x, false, saveProgress); % pass whether to save progress 
[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
pCur = pVec;
% highSensIdx = [12,31,34,39,41,42,43,59:75,89,91,95];
highSensIdx = [1,3,4,5,7,8,9,11,12,13,16,17,19,22,25,26,27,29,30,31,34,36,38,39,41,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,70,71,73,74,75,76,82,84,85,87,88,89,92,93,94,95,96,97,98,99,100,101,102,103,104,105];
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
pCur = pCur(:) .* p0;
SkelMuscleObj2(pCur)
filename = "PSO_" + string(datetime("today")) +".mat";
% save(fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename'));
% save(fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename));
save(filename)

end