%% start parameter optimization
% Importing parameters 
highSensIdx = 1:105;%[12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
% Setting bounds for parameters
lb = 0.5*ones(length(highSensIdx),1); 
ub = 2*ones(length(highSensIdx),1);

% Particle Swarm Optimization
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub);
username = getenv('USER');
fname = sprintf('/tscc/lustre/ddn/scratch/%s/psofig_%s.mat',username,datetime('today'));
saveas(gcf, fname);
