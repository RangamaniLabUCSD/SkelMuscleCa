%% make plot following parameter estimation
load PSO_06-Feb-2025.mat pSol highSensIdx
param = importdata('InputParam1.xlsx'); % load default parameters
p0 =  param.data;
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);
% load pBest.mat pVec
% pVec(45) = 0;
[objVal, simSaved] = SkelMuscleObj(pPSO, true);

%% make plot of freq-dependent SOCE behavior (further testing below)
load PSO_19-Dec-2024.mat pSol
[objVal, simSaved] = SkelMuscleObj2(pSol, true);

%% Test SOCE vs no SOCE steady states after adjusting parameters
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
    0;%387;        % yinit(14) is the initial condition for 'CaParv'
    0;%1020;       % yinit(15) is the initial condition for 'MgParv'
    0.3632;     % yinit(16) is the initial consition for 'CATP'
    0;%10.004;     % yinit(17) is the initial condition for 'CaTrop'
    0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    0;	    	% yinit(19) is the initial condition for 'D_2'
    0;	    	% yinit(20) is the initial condition for 'Pre_Pow'
    0;	    	% yinit(21) is the initial condition for 'Post_Pow'
    0;	    	% yinit(22) is the initial condition for 'MgATP'
    8000;       % yinit(23) is the initial condition for 'ATP'
    3000;        % yinit(24) is the initial condition for 'p_i_SR'
    0;           % yinit(25) is the initial condition for 'PiCa'
    3000;        % yinit(26) is the initial condition for 'Pi_Myo'
    1300;        % c_o (µM)
    147000.0;   % Na_o (µM)
    4000.0;      % K_o (µM)
    128000.0;   % Cl_o (µM)
    15000; % free CSQ
    ];
ECVals = [1300;        % c_o (µM)
        147000.0;   % Na_o (µM)
        4000.0;      % K_o (µM)
        128000.0  ];
juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
yinit = [yinit(juncLocLogic); yinit(bulkLocLogic)];

param = importdata('InputParam1.xlsx'); % load default parameters
p0 =  param.data;

% load PSO_22-Nov-2024_20size.mat pSol
load PSO_18-Dec-2024.mat pSol highSensIdx
% highSensIdx = 1:105;%[1,5,15,16,19,22,24,30,32,34,35,43,74,83,91,92];
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);
pPSO(45) = 0.002*.5676;
% load pBest.mat pVec
% pVec(45) = 0;
% pPSO = pVec;

% if desired, perturb parameters to test effects
% pPSO(45) = 0*pPSO(45); % Na leak
% pPSO(19) = 1*pPSO(19); % Na channel
% number 69 is the rate of phosphate transport into SR
% pPSO(69) = pPSO(69)*100
% pPSO(21) = 0.1*pPSO(21); % NK pump

% try changing PMCA and SERCA systematically
% pPSO(34) = 1*pPSO(34); % K_SERCA
% pPSO(42) = 1*pPSO(42); % nu_SERCA
% pPSO(31) = 1*pPSO(31); % K_PMCA
% pPSO(43) = 4*pPSO(43); % g_PMCA
% pPSO(20) = 0.5*pPSO(20); % g_NCX
% pPSO(44) = 2*pPSO(44); % SR calcium leak
% pPSO(100) = 1.0*pPSO(100); % D_ion
% 
% pPSO(91) = 0.0*pPSO(91); % g0 DHPR
% pPSO(89) = 0.1*pPSO(89); % SL leak calcium
% pPSO(92) = 10*pPSO(92); % j0_RyR
% 
% pPSO(52) = 0.5*pPSO(52); %Parv conc
% pPSO(41) = 0.01*pPSO(41); % tau SOCE
% pPSO(96) = 0.2*pPSO(96); % total SR buffer
% yinit(sum(juncLocLogic(1:31))) = pPSO(96)*0.5;
% yinit(sum(juncLocLogic) + sum(bulkLocLogic(1:31))) = pPSO(96)*0.5;
% 
% pPSO(65) = 1*pPSO(65); % h0

% compute steady state solution without SOCE or phosphate accumulation
pPSO0 = pPSO;
pPSO0(95) = 0; % no SOCE
[TimeSS_noSOCE,ySS_noSOCE, ~, fluxes, currents] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, false);
% pPSO(12) = 0.25;
pPSO(12) = pPSO(12)*ySS_noSOCE(end,2); % set c_ref according to SS c_SR 
% pPSO(95) = pPSO(95)*10; % gSOCE
[TimeSS_withSOCE,ySS_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, false);

% compute max cross bridge engagement - expt 10 is only crossbridge testing
% (all other dydt = 0)
yinit_maxCa = ySS_withSOCE(end,:);
yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
[Time_maxCa, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa, pPSO, tic, 10, false);
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
maxCrossBridge = Y_maxCa(end,crossbridgeIdx);

vol_Fiber = pi * (20 ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003 ;
vol_myo = 0.95 * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio ;
diffusion_length = pPSO(99);%0.05;
vol_myoJ = SA_TT*diffusion_length;
JFrac = vol_myoJ / vol_myo;
BFrac = 1 - JFrac;

yinf_noSOCE = ySS_noSOCE(end,:);
yinf_withSOCE = ySS_withSOCE(end,:);

%% test single frequency and plot currents and fluxes
tSol = [0, 0.5]; % 420 is max time for HIIT
expt = 1;
phosphateAccum = true;
[Time_withSOCE,Y_withSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, 50, yinf_withSOCE, pPSO, tic, 1, phosphateAccum); % compute time-dependent solution
[Time_noSOCE,Y_noSOCE,~,fluxes2,currents2] = SkelMuscleCa_dydt(tSol, 50, yinf_noSOCE, pPSO, tic, 2, phosphateAccum); % compute time-dependent solution
% fluxes in µM/s = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
%                   LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];% 
% current in pA/s = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                %    I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% plot sodium, KDR, and KIR currents
figure
subplot(3,1,1)
plot(Time_withSOCE, currents(:,10)) % in TT
hold on
plot(Time_withSOCE, currents(:,13+10)) % in bulk (SL)
plot(Time_withSOCE, currents(:,10)+currents(:,13+10))
legend('TT', 'SL', 'Total')
title('Na')
subplot(3,1,2)
plot(Time_withSOCE, currents(:,4)) % in TT
hold on
plot(Time_withSOCE, currents(:,13+4)) % in bulk (SL)
plot(Time_withSOCE, currents(:,4)+currents(:,13+4))
legend('TT', 'SL', 'Total')
title('KDR')
subplot(3,1,3)
plot(Time_withSOCE, currents(:,5)) % in TT
hold on
plot(Time_withSOCE, currents(:,13+5)) % in bulk (SL)
plot(Time_withSOCE, currents(:,5)+currents(:,13+5))
legend('TT', 'SL', 'Total')
title('KIR')
% plot RyR, SERCA, and PMCA fluxes
% multiply all fluxes by JFrac or BFrac to represent in terms of total
% concentration change
figure
subplot(3,1,1)
plot(Time_withSOCE, fluxes(:,6)*JFrac) % from SRJ to MJ
hold on
plot(Time_withSOCE, fluxes(:,8+6)*BFrac) % from bulk SR to bulk myo
plot(Time_withSOCE, fluxes(:,6)*JFrac+fluxes(:,8+6)*BFrac)
legend('SRJ', 'SRB', 'Total')
title('RyR')
subplot(3,1,2)
plot(Time_withSOCE, fluxes(:,7)*JFrac) % from SRJ to MJ
hold on
plot(Time_withSOCE, fluxes(:,8+7)*BFrac) % from bulk SR to bulk myo
plot(Time_withSOCE, fluxes(:,7)*JFrac+fluxes(:,8+7)*BFrac)
legend('SRJ', 'SRB', 'Total')
title('SERCA')
subplot(3,1,3)
plot(Time_withSOCE, fluxes(:,5)*JFrac) % from TT to MJ
hold on
plot(Time_withSOCE, fluxes(:,8+5)*BFrac) % from SL to bulk myo
plot(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,8+5)*BFrac)
legend('TT', 'SL', 'Total')
title('PMCA')

%%
figure
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
hold on
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)

%% test range of frequencies for SOCE vs no SOCE
tSol = [0, 0.5]; % 420 is max time for HIIT
freq = [1, 25, 50, 75, 100, 125, 150, 175, 200, 250];
forceRatio = zeros(size(freq));
peakForce = zeros(2, length(freq));
endPeakRatio = zeros(2, length(freq));
figure
for i = 1:length(freq)
    % expt 1: HIIT stim, with SOCE, expt 2: HIIT stim, no SOCE
    expt = 1;
    phosphateAccum = true;
    [Time_withSOCE,Y_withSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq(i), yinf_withSOCE, pPSO, tic, 1, phosphateAccum); % compute time-dependent solution
    [Time_noSOCE,Y_noSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq(i), yinf_noSOCE, pPSO, tic, 2, phosphateAccum); % compute time-dependent solution
    
    % create plots
    colors = colororder;
    subplot(3,1,1)
    % final entry in color spec gives the opacity
    plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac, 'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac, 'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('Cytosolic calcium conc. (µM)')
    title('Cytosolic calcium')
    legend('With SOCE', 'No SOCE')
    subplot(3,1,2)
    plot(Time_withSOCE, Y_withSOCE(:,2), 'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,2), 'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('SR calcium conc. (µM)')
    title('SR calcium')
    subplot(3,1,3)
    plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge, 'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge, 'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('Relative force')
    title('Force')
    drawnow
    forceRatio(i) = max(Y_withSOCE(:,crossbridgeIdx))/max(Y_noSOCE(:,crossbridgeIdx));
    peakForce(1,i) = max(Y_withSOCE(:,crossbridgeIdx));
    peakForce(2,i) = max(Y_noSOCE(:,crossbridgeIdx));
    endPeakRatio(1,i) = Y_withSOCE(end,crossbridgeIdx)/peakForce(1,i);
    endPeakRatio(2,i) = Y_noSOCE(end,crossbridgeIdx)/peakForce(2,i);
end

peakForce = peakForce / max(peakForce(:));

figure
plot(freq, peakForce)

%% test agreement and make plot for adjusted parameter set
pPSO(12) = pSol(12);
[objVal, simSaved] = SkelMuscleObj2(pPSO, true);