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
            0.3632;    % yinit(16) is the inital consition for 'CATP'
            10.004;    % yinit(17) is the initial condition for 'CaTrop'
            ];
 
% Importing parameters 
param = importdata('InputParam1.xlsx');
p =  param.data;

% Setting bounds for parameters
lb = 0.8*ones(length(p),1); 
ub = 1.25*ones(length(p),1);
% limits for NCX, SERCA, PMCA
lb([20, 42, 43]) = 0.25;
ub([20, 42, 43]) = 1.0;
% limits for SR Ca2+ leak
lb(44) = 0.1;
ub(44) = 0.4;
% limits for sodium leak through SL
lb(45) = 0;
ub(45) = 2.5;

tSpan = [0 1];
Createplot = 1; %Logical input of 0 or 1. 0 for not plotting any outputs and 1 for generatings plots.
timer1 = tic;

% Particle Swarm Optimization 
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst(tSpan,lb, ub, yinit, p,Createplot);

toc(timer1)
filename_fig = "BestFnva"+ date + ".jpg";
saveas(gcf,filename_fig)
