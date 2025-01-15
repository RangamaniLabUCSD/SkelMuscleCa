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
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 3, 'SwarmSize', 15);

% pSol Results
load PSO_18-Dec-2024.mat pSol
param = importdata('InputParam1.xlsx');
p0 =  param.data;
highSensIdxResults = [2,6,10,14,15,18,20,21,23,24,28,32,33,35,37,40,42,43,45,69,72,77,78,79,80,81,83,86,90,91];
p0(highSensIdxResults) = p0(highSensIdxResults) .* pSol(:);
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

if saveProgress
    % enter your chosen path here, be sure it exists!
    progressPath = '/tscc/lustre/ddn/scratch/jhamid/SkelMuscleCaData';
    if isfolder(progressPath)
        objTest = inf;
        save(fullfile(progressPath,'objTest.mat'),'objTest');
    else
        progressPath = '';
        objTest = inf;
        save(fullfile(progressPath,'objTest.mat'),'objTest');
    end
end

fAnon = @(x)SkelMuscleObj2(x, false, progressPath); % pass path to save progress
[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
pCur = pVec;
% phosphate = 68:72,74,75 ;  CaEfflux_SOCE = [12,31,34,41:44,89,95]; Crossbridge_Cycle = [55:57,59:67,73,84,93,94];
highSensIdx = [12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
pCur = pCur(:) .* p0;
SkelMuscleObj2(pCur)
filename = "PSO_" + string(datetime("today")) +".mat";
% save(fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename'));
% save(fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename));
save(filename)

end