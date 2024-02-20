%% start parameter optimization
tSpan = [0 1];

 yinit = [
            0.0122; 	% yinit(1) is the initial condition for 'SOCEProb'
            1500.0;		% yinit(2) is the initial condition for 'c_SR'
            0.9983;		% yinit(3) is the initial condition for 'h_K'
            0.9091;		% yinit(4) is the initial condition for 'w_RyR'            
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
 
% Parameter values
param = importdata('InputParam1.xlsx');
p =  param.data; 

% ClampCurrent = p(1) ;K_S = p(2) ;delta = p(3) ;beta_m0 = p(4) ;K_betam = p(5) ;alpha_m0 = p(6) ;K_alpham = p(7) ;K_RyR = p(8) ;f_RyR = p(9) ;
%p = [-25000, 1000000, 0.4, 1380, 18, 288, 10, 4.5, 0.2]' ;%.* pSol_LM';

lb = 0.5*ones(length(p),1); 
ub = 2*ones(length(p),1);

%% LM + ode15s
timer1 = tic;
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst_LM(tSpan,lb, ub, yinit, p);
toc(timer1)

saveas(gcf,'BestFnva_2_16_1_LM.jpg')


