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

 %set swarmsize to 30 for TSCC and stall iter to 50
psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 50, 'SwarmSize', 30);
param = importdata('InputParam1.xlsx');
p0 =  param.data;
highSensIdx = 1:105;%[12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
% Test the objective function for default values
SkelMuscleObj(ones(length(highSensIdx),1))

delete(gcp('nocreate'))
if psOptions.UseParallel 
    parpool(30) %% **CHANGE SWARMSIZE!** and save file location
end
saveProgress = true;

username = getenv('USER');
if saveProgress
    % enter your chosen path here, be sure it exists!
    progressPath = sprintf('/tscc/lustre/ddn/scratch/%s/SkelMuscleCaData',username);
    if isfolder(progressPath)
        objTest = inf;
        save(fullfile(progressPath,'objTest.mat'),'objTest');
    else
        progressPath = pwd;
        objTest = inf;
        save(fullfile(progressPath,'objTest.mat'),'objTest');
    end
    fAnon = @(x)(SkelMuscleObj(x, false, progressPath) + SkelMuscleObj2(x, false, progressPath)/100); % pass path to save progress
else
    fAnon = @(x)(SkelMuscleObj(x, false) + SkelMuscleObj2(x, false)/100);
end

[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
pCur = ones(size(p0));
% phosphate = 68:72,74,75 ;  CaEfflux_SOCE = [12,31,34,41:44,89,95]; Crossbridge_Cycle = [55:57,59:67,73,84,93,94];
highSensIdx = 1:105;%[12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
pCur(highSensIdx) = pSol(:) .* pCur(highSensIdx)';
SkelMuscleObj(pCur)
filename = "PSO_" + string(datetime("today")) +".mat";
save(fullfile(sprintf('/tscc/lustre/ddn/scratch/%s/',username),filename));
end    