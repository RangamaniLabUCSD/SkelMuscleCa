%% Script for running Morris analysis (uqlab must be in the MATLAB path)
rng(100,'twister')
withMito = true;
uqlab
delete(gcp('nocreate'))
ModelOpts.mFile = 'SkelMuscleCa_SAOutput';
myModel = uq_createModel(ModelOpts);
load Data/p0Struct.mat p0Struct
if withMito
    p0Struct.names = {p0Struct.names{:}, 'volFrac_mitoJ','volFrac_mitoB','vol_SA_mito',...
                      'ADP0_mito', 'ATP0_mito', 'Pi0_mito'};
    p0Struct.data = [p0Struct.data; 1; 1; 1; 1; 1; 1];
end
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
MorrisSensOpts.Morris.Cost = 4e4;
MorrisAnalysis = uq_createAnalysis(MorrisSensOpts);
username = getenv('USER');
fname1 = sprintf('/tscc/lustre/ddn/scratch/%s/SkelMuscleCaData', username);
fname2 = sprintf('/expanse/lustre/scratch/%s/temp_project/SkelMuscleCaData', username);
if isfolder(fname1)
    fname = sprintf('%s/MorrisResults_4e4_%s.mat', fname1, datetime('today'));
elseif isfolder(fname2)
    fname = sprintf('%s/MorrisResults_4e4_%s.mat', fname1, datetime('today'));
else
    fname = sprintf('MorrisResults_4e4_%s.mat', datetime('today'));
end
save(fname,'MorrisAnalysis');
endtimer = toc(starttimer);