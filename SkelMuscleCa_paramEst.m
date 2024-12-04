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
psOptions = optimoptions('particleswarm','UseParallel',false,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 10, 'SwarmSize', 15);

% pSol Results
% load PSO_22-Aug-2024.mat pSol
param = importdata('InputParam1.xlsx');
p0 =  param.data;
% To check default parameter behavior
numParam = length(p0);
pVec = ones(1,numParam) ; 
pCur = pVec;
pCur = pCur(:) .* p0;
SkelMuscleObj(pCur)

delete(gcp('nocreate'))
if psOptions.UseParallel
    parpool(3) %% **CHANGE SWARMSIZE!** and save file location
end
saveProgress = true;
if saveProgress
    % enter your chosen path here, be sure it exists!
    progressPath = '/tscc/lustre/ddn/scratch/eafrancis/SkelMuscleCaData';
    if isfolder(progressPath)
        objTest = inf;
        save(fullfile(progressPath,'objTest.mat'),'objTest');
    else
        progressPath = '';
    end
end
fAnon = @(x)SkelMuscleObj(x, false, progressPath); % pass path to save progress
[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
pCur = pVec;
% highSensIdx = [12,31,34,39,41,42,43,59:75,89,91,95];
highSensIdx = 1:105;
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
pCur = pCur(:) .* p0;
SkelMuscleObj(pCur)
filename = "PSO_" + string(datetime("today")) +".mat";
% save(fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename'));
% save(fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename));
save(filename)

end