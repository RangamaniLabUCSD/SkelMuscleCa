%% start parameter optimization
tSpan = [0 1];
freqVals = 0;

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
 
param = importdata('InputParam1.xlsx');
p = param.data; % .* pSol_est'; % Parameter values
% load LM_PSO_2_14_1.mat pSol_LM
% ClampCurrent = p(1) ;K_S = p(2) ;delta = p(3) ;beta_m0 = p(4) ;K_betam = p(5) ;alpha_m0 = p(6) ;K_alpham = p(7) ;K_RyR = p(8) ;f_RyR = p(9) ;
% p = [-25000, 1000000, 0.4, 1380, 18, 288, 10, 4.5, 0.2]' .* pSol_LM';
% p(1) = 1.5 *-25000;
% % p(2) = 1.5 * 1000000;
% p(3) = 0.5 * 0.4;  
% p(4) = 1.5 * 1380;
% p(5) = 1.5 * 18;
lb = 0.5*ones(length(p),1); 
ub = 2*ones(length(p),1);
%% Test - Check how long it takes the system to reach SS 

%freq = [100, 100, 67, 67,67, 60, 60, 60,60];
% p(1) = 0.9*p(1);
% p(2) = 0.9*p(2);
% p(3) = 0.9*p(3);
% p(4) = 0.9*p(4);
% p(5) = 1.2 *p(5);
% p(6) = 0.985 * p(6);
% p(7) = 1.0 * p(7);
% p(8) = 0.5 * p(8);
% p(9) = 1.0 * p(9);
% 
% for n = 1 %:5 % length(yinit)
%     initialGuess = yinit';
%     options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIter', 1000);
%     lb1 = 0.5 * initialGuess;
%     ub1 = 2 * initialGuess;
%     lb1(5) = 2 * initialGuess(5);
%     ub1(5) = 0.5 * initialGuess(5);
% 
%     [y_LM] = lsqnonlin(@(initialGuess) objectivefn(initialGuess,n,p),initialGuess, lb1, ub1, options);
% 
%     [tSS,ySS] = SkelMuscleCa1_SS([0 1000],0, 0, y_LM, p,tic,1);
%     yinf = ySS(end,:);
% 
%     % semilogy(tSS,ySS(:,n),'LineWidth',2);
%     % hold on
%     % xlabel('Time (s)');
%     %t = [0 0.1];
%     %p(6) = 1.5*p(6);
%     %p = p;
%     t = [0.0554782608695652	0.235507246400000	0.119565217391304	0.0499324324324324	0.0978927203065133	0.00561386138600000	0.0296590909090909	0.0369811320754717	0.00600000000000000];
%     [Time,y] = SkelMuscleCa1_SS([0 t(n)],freq(n), 0, yinf, p,tic,n);
% 
%     figure
%     yyaxis left
%     plot(Time,y(:,8),'b','LineWidth',2);
%     ylabel('[Ca^{2+}]_i (uM)')
%     yyaxis right
% 
%     plot(Time,y(:,5),'r','LineWidth',2);
%     ylabel('Voltage_{PM} (mV)')
%     xlabel('Time (s)')
% end

%% LM + ode15s
timer1 = tic;
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst_LM(tSpan,lb, ub, yinit, p);
toc(timer1)
save('LM_PSO_2_14_4.mat')
saveas(gcf,'BestFnva_2_14_4_LM.jpg')

%% Edits to fix RyR and DHPR
%Step 1: Remove w-DHPR and reduce yinit/dydt to 16. w_DHPR = w_RyR for
%J_DHPR. Nf reduced to 16
%Step 2: equate f, K, and L for DHPR with RyR. 
% Remove f_DHPR from param list

%% Error
% Warning: Derivative finite-differencing step was artificially reduced to be within bound constraints. This may adversely affect convergence. Increasing distance between bound constraints, in dimension 1, to be at least 3.9056e-09 may improve results.
% > In fwdFinDiffInsideBnds
% In finitedifferences
% In computeFinDiffGradAndJac
% In levenbergMarquardt (line 107)
% In lsqncommon (line 187)
% In lsqnonlin (line 260)
% In SkelMuscleCa_paramEst_LM/pToObj (line 62)
% In parallel_function>make_general_channel/channel_general (line 837)
% In remoteParallelFunction (line 67)

