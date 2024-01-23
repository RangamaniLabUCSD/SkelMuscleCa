%% start parameter optimization
clc
tSpan = [0 1];
freqVals = 0;

 yinit = [
            0.0122; 	% yinit(1) is the initial condition for 'SOCEProb'
            1500.0;		% yinit(2) is the initial condition for 'c_SR'
            0.9983;		% yinit(3) is the initial condition for 'h_K'
            0.9091;		% yinit(4) is the initial condition for 'w_RyR'
            0.9091;		% yinit(5) is the initial condition for 'w_DHPR'
            -88.0;		% yinit(6) is the initial condition for 'Voltage_PM'
            14700.0;	% yinit(7) is the initial condition for 'Na_i'
            5830.0;		% yinit(8) is the initial condition for 'Cl_i'
            0.1;		% yinit(9) is the initial condition for 'c_i'
            0.003;		% yinit(10) is the initial condition for 'n'
            0.0128;		% yinit(11) is the initial condition for 'm'
            0.8051;		% yinit(12) is the initial condition for 'h'
            0.8487;		% yinit(13) is the initial condition for 'S'
            154500.0;	% yinit(14) is the initial condition for 'K_i'
            387;        % yinit(15) is the initial condition for 'CaParv'
            1020;       % yinit(16) is the initial condition for 'MgParv'
            0.3632;    % yinit(17) is the inital consition for 'CATP'
            ];
 
%param = importdata('InputParam1.xlsx');
%p = param.data;% .* pSol_est'; % Parameter values
%load LM_PSO_1_22.mat pSol_LM 
p = [4380,0.2,0.2,10,18,4.5,1000,13000,10000,20]';% .* pSol_LM';% 
lb = 0.5*ones(10,1); 
ub = 2*ones(10,1);

%% LM + ode15s
timer1 = tic;
pSol_LM = SkelMuscleCa_paramEst_LM(tSpan,lb, ub, yinit, p);
timerstop1 = toc(timer1)
save('LM_PSO_1_23.mat')
%saveas(gcf,'BestFnva_1_22V.jpg')
exportapp(gcf,'BestFnva_1_23V.jpg');
%% ode15s
% timer2 = tic;
% pSol_ode = SkelMuscleCa_paramEst(tSpan,lb, ub, yinit, p);
% timerstop2 = toc(timer2)
% save('ode_PSO_1_8.mat')
% saveas(gcf,'BestFnva_1_8_ode.jpg')



