%% start parameter optimization
yinit = [
    0.0122; 	% yinit(1) is the initial condition for 'SOCEProb'
    1500.0;		% yinit(2) is the initial condition for 'c_SR'
    0.9983;		% yinit(3) is the initial condition for 'h_K'
    0.9091;		% yinit(4) is the initial condition for 'w_RyR'
    -88.0;		% yinit(5) is the initial condition for 'Voltage_PM'
    14700.0;	% yinit(6) is the initial condition for 'Na_i'
    5830.0;		% yinit(7) is the initial condition for 'Cl_i'
    0.1;		% yinit(8) is the initial condition for 'c_i'
    0.003;		% yinit(9) is the initial condition for 'n'
    0.0128;		% yinit(10) is the initial condition for 'm'
    0.8051;		% yinit(11) is the initial condition for 'h'
    0.8487;		% yinit(12) is the initial condition for 'S'
    154500.0;	% yinit(13) is the initial condition for 'K_i'
    387;        % yinit(14) is the initial condition for 'CaParv'
    1020;       % yinit(15) is the initial condition for 'MgParv'
    0.3632;     % yinit(16) is the inital consition for 'CATP'
    10.004;     % yinit(17) is the initial condition for 'CaTrop'
    0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    0;	    	% yinit(19) is the initial condition for 'D_0'
    0;	    	% yinit(20) is the initial condition for 'D_1'
    0;	    	% yinit(21) is the initial condition for 'D_2'
    0;	    	% yinit(22) is the initial condition for 'Pre_Pow'
    0;	    	% yinit(23) is the initial condition for 'Post_Pow'
    0;	    	% yinit(24) is the initial condition for 'MgATP'
    8000;       % yinit(25) is the initial condition for 'ATP'
    3000        % yinit(26) is the initial condition for 'p_i_SR'
    0           % yinit(27) is the initial condition for 'PiCa'
    3000        % yinit(28) is the initial condition for 'Pi_Myo'
    ];

% Importing parameters 
param = importdata('InputParam1.xlsx');
p =  param.data;
highSensIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,30,31,32,33,35,40,41,42,43,44,45,52,71,72,73,78,79,80,81,82,83,84,85,86,89,92,93]; % a vector listing the indices of all parameters we are still including (higher sensitivity values)
% Setting bounds for parameters
lb = 0.8*ones(length(highSensIdx),1); 
ub = 1.25*ones(length(highSensIdx),1);
% limits for NCX, SERCA, PMCA
% lb([20, 42, 43]) = 0.25;
% ub([20, 42, 43]) = 1.0;
lb(highSensIdx == 43) = 0.25;
ub(highSensIdx == 43) = 1.0;
% limits for SR Ca2+ leak
lb(highSensIdx == 44) = 0.1;
ub(highSensIdx == 44) = 0.4;
% % limits for sodium leak through SL
lb(highSensIdx == 45) = 0;
ub(highSensIdx == 45) = 2.5;
% timer1 = tic;

% Particle Swarm Optimization 
Createplot = 0; %Logical input of 0 or 1. 0 for not plotting any outputs and 1 for generatings plots.
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub, yinit, p,Createplot);


% toc(timer1)
filename_fig = "BestFnva"+ date + ".jpg";
% saveas(gcf,filename_fig)
% saveas(gcf,fullfile('/tscc/lustre/ddn/scratch/jhamid/',filename_fig));
saveas(gcf,fullfile('C:/Users/Juliette/Documents/MATLAB/SkelMuscle/',filename_fig));
