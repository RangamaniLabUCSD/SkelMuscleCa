%% Function for estimating parameters using Particle swarm estimation
% Input:
%       lb - Lower bounds for estimation
%       ub - upper bounds for estimation
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
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 10, 'SwarmSize', 10);
delete(gcp('nocreate'))
if psOptions.UseParallel
    parpool(5) %% **CHANGE SWARMSIZE!** and save file location
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
    % fAnon = @(x)(SkelMuscleObj(x, false, progressPath) + SkelMuscleObj2(x, false, progressPath)/20); % pass path to save progress
    fAnon = @(x)(LinCombObj(x, progressPath)); % pass path to save progress
else
    % fAnon = @(x)(SkelMuscleObj(x, false) + SkelMuscleObj2(x, false)/20);
    fAnon = @(x)(LinCombObj(x, '')); % pass path to save progress
end

[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
filename = "PSO_" + string(datetime("today")) +".mat";
save(fullfile(sprintf('/tscc/lustre/ddn/scratch/%s/',username),filename));

    function objVal = LinCombObj(pVec, progressPath)
        if isfolder(progressPath)
            saveProgress = true;
        end
        objVal = SkelMuscleObj(pVec, false);% + SkelMuscleObj2(pVec, false)/20;
        if saveProgress
            try
                load(fullfile(progressPath,'objTest.mat'),'objTest')
                if objVal < objTest
                    objTest = objVal;
                    save(fullfile(progressPath,'objTest.mat'),'objTest')
                    save(fullfile(progressPath,'pBest.mat'),'pVec')
                end
            catch
                fprintf('file was busy I guess\n')
            end
        end
    end
end