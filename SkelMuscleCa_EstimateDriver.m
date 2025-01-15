%% start parameter optimization
% Importing parameters 
highSensIdx = [12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];

% p(95) = p(95)*5; % increase SOCE baseline
% highSensIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,32,33,34,35,40,42,43,44,45,50,52,74,76,77,78,79,80,81,82,83,89,90,91,92];
% highSensIdx = [1,5,15,16,19,22,24,30,32,34,35,43,74,83,91,92];
% highSensIdx = [2,4,5,6,8,10,14,15,16,24,25,26,28,29,32,33,35,37,40,42,43,60,68,74,77,78,80,81,83,90,91,92];
% highSensIdx = [12,31,34,39,41,42,43,59:75,89,91,95];
% Setting bounds for parameters
lb = 0.5*ones(length(highSensIdx),1); 
ub = 2*ones(length(highSensIdx),1);
% limits for NCX, SERCA, PMCA
% lb([20, 42, 43]) = 0.25;
% ub([20, 42, 43]) = 1.0;
% lb(highSensIdx == 20) = 0.25;
% ub(highSensIdx == 20) = 1.0;
% lb(highSensIdx == 42) = 0.25;
% ub(highSensIdx == 42) = 1.0;
% lb(highSensIdx == 43) = 0.1;
% ub(highSensIdx == 43) = 0.4;
% limits for SR Ca2+ leak
% lb(highSensIdx == 44) = 0.1;
% ub(highSensIdx == 44) = 0.4;
% % limits for sodium leak through SL
% lb(highSensIdx == 45) = 0;
% ub(highSensIdx == 45) = 2.5;
% timer1 = tic;

% Particle Swarm Optimization
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub);
% saveas(gcf,filename_fig)
saveas(gcf,fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename_fig));
% saveas(gcf,fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename_fig));
