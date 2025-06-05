clc 
clearvars
rng(100,'twister')
uqlab
delete(gcp('nocreate'))
ModelOpts.mFile = 'SkelMuscleCa_SAOutput';
myModel = uq_createModel(ModelOpts);
load Data/p0Struct.mat p0Struct
paramNames = p0Struct.names;
numParam = length(paramNames);
starttimer = tic;
% these specifications designate uniform distributions ranging from about 0.5 to 2
for i = 1:numParam
    InputOpts.Marginals(i).Type = 'Uniform';
    InputOpts.Marginals(i).Name = paramNames{i};
    InputOpts.Marginals(i).Parameters = [0.5, 2.0]; 
end
myInput = uq_createInput(InputOpts);

%% Morris SA
MorrisSensOpts.Type = 'Sensitivity';
MorrisSensOpts.Method = 'Morris';
MorrisSensOpts.Morris.Cost = 1e4;
MorrisAnalysis = uq_createAnalysis(MorrisSensOpts);
username = getenv('USER');
fname = sprintf('/tscc/lustre/ddn/scratch/%s/MorrisResults_1e4_%s.mat',username,datetime('today'));
save(fname,'MorrisAnalysis');
endtimer = toc(starttimer);