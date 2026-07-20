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

function [pSol,fval,exitflag] = SkelMuscleCa_paramEst(lb,ub,withMito,initVec)

%set swarmsize to 30 for TSCC and stall iter to 20
try
    psOptions = optimoptions('particleswarm','UseParallel',true,...%'HybridFcn',@fmincon,...
        'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 20, 'SwarmSize', 30,...
        'InitialPoints', initVec);
catch
    psOptions = optimoptions('particleswarm','UseParallel',true,...%'HybridFcn',@fmincon,...
        'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 20, 'SwarmSize', 30,...
        'InitialSwarmMatrix', initVec);
end
delete(gcp('nocreate'))
if psOptions.UseParallel
    parpool(30) %% **CHANGE SWARMSIZE!** and save file location
end
saveProgress = true;

username = getenv('USER');
if saveProgress
    if withMito
        if length(ub) < 106
            tagString = 'withMito_Vonly';
        else
            tagString = 'withMito';
        end
    else
        if length(ub) < 106
            tagString = 'noMito_Vonly';
        else
            tagString = 'noMito';
        end
    end
    % enter your chosen path here, be sure it exists!
    progressPath1 = sprintf('/tscc/lustre/ddn/scratch/%s/SkelMuscleCaData',username);
    progressPath2 = sprintf('/expanse/lustre/scratch/%s/temp_project/SkelMuscleCaData',username);
    if isfolder(progressPath1)
        progressPath = progressPath1;
        objTest = inf;
        save(fullfile(progressPath,sprintf('objTest_%s.mat', tagString)),'objTest');
    elseif isfolder(progressPath2)
        progressPath = progressPath2;
        objTest = inf;
        save(fullfile(progressPath,sprintf('objTest_%s.mat', tagString)),'objTest');
    else
        progressPath = pwd;
        objTest = inf;
        save(fullfile(progressPath,sprintf('objTest_%s.mat', tagString)),'objTest');
    end
    % fAnon = @(x)(SkelMuscleObj(x, false, progressPath) + SkelMuscleObj2(x, false, progressPath)/20); % pass path to save progress
    fAnon = @(x)(LinCombObj(x, progressPath)); % pass path to save progress
else
    % fAnon = @(x)(SkelMuscleObj(x, false) + SkelMuscleObj2(x, false)/20);
    fAnon = @(x)(LinCombObj(x, '')); % pass path to save progress
end

[pSol,fval,exitflag] = particleswarm(fAnon,length(lb),lb,ub,psOptions);
% filename = "PSO_" + string(datetime("today")) +"_withMito.mat";
filename = sprintf('PSO_%s_%s.mat', string(datetime("today")), tagString);
save(fullfile(progressPath,filename));

    function objVal = LinCombObj(pVec, progressPath)
        if isfolder(progressPath)
            saveProgress = true;
        end
        objVal = SkelMuscleObj(pVec, false,'',withMito);% + SkelMuscleObj2(pVec, false)/20;
        if saveProgress
            try
                load(fullfile(progressPath,sprintf('objTest_%s.mat', tagString)),'objTest')
                if objVal < objTest
                    objTest = objVal;
                    save(fullfile(progressPath,sprintf('objTest_%s.mat', tagString)),'objTest')
                    save(fullfile(progressPath,sprintf('pBest_%s.mat', tagString)),'pVec')
                end
            catch
                fprintf('file was busy I guess\n')
            end
        end
    end
end