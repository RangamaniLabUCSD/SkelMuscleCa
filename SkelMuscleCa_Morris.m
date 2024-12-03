clc 
clearvars
rng(100,'twister')
uqlab
delete(gcp('nocreate'))
ModelOpts.mFile = 'SkelMuscleCa_SAOutput';
myModel = uq_createModel(ModelOpts);
input = importdata('InputParam1.xlsx');
paramNames = input.textdata;
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
MorrisSensOpts.Morris.Cost = 1e5;
MorrisAnalysis = uq_createAnalysis(MorrisSensOpts);
% save('/tscc/lustre/ddn/scratch/jhamid/MorrisResults9-16.mat','MorrisAnalysis');
save('MorrisResults_1e5.mat','MorrisAnalysis');
endtimer = toc(starttimer);
%uq_display(MorrisAnalcysis)