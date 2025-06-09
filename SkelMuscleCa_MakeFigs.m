%% SkelMuscleCa_MakeFigs script
% This file contains all the necessary code to generate figures for the paper 
% "Systems modeling reveals that store-operated calcium entry modulates
% force and fatigue during exercise"

%% Plot results of estimation
load Data/pSol_allParam.mat pSol
SkelMuscleObj(pSol, true)

%% Morris plots for Fig 2 and Fig S2
try uqlab
catch
    error('must first add UQLab to path before starting this section of code')
end
load Data/MorrisResults.mat MorrisAnalysis
graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion",...
    "SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force",...
    "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force",...
    "Avg Voltage", "AP Width"};
expt_names = repmat(graph_names,1,4)';
expt_names = [expt_names; "Objective Value 1"];
muVec = zeros(size(MorrisAnalysis.Results.MuStar));

ActionPotential = [1:6, 8:11, 13,14,16:30,33,36:40,45,76:82, 85:88];
CaInflux_SR_Release = [7, 15,32,35,83,90:92];
CaBuffering = [46:54,58,68:72,74,75,96,97,98];
CaEfflux_SOCE = [12,31,34,41:44,89,95];
Crossbridge_Cycle = [55:57,59:67,73,84,93,94,106];
Diffusion = 98:105; % parameters for diffusion between junctional and bulk compartments
qoiList = [42+9,42+13, 8, 11]; % voltage qois from expt 8, calcium qois from expt 3
sensIdx = cell(length(qoiList),1);

for k = qoiList
    sensIdx{qoiList==k} = [];
    figure
    scatter(MorrisAnalysis.Results.MuStar(ActionPotential,k),...
        MorrisAnalysis.Results.Std(ActionPotential,k),'filled',...
        'MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(CaInflux_SR_Release,k),...
        MorrisAnalysis.Results.Std(CaInflux_SR_Release,k),'filled',...
        'MarkerFaceColor',[0.4660 0.6740 0.1880])
    scatter(MorrisAnalysis.Results.MuStar(CaBuffering,k),...
        MorrisAnalysis.Results.Std(CaBuffering,k),'filled',...
        'MarkerFaceColor',[0 0.4470 0.7410])
    scatter(MorrisAnalysis.Results.MuStar(CaEfflux_SOCE,k),...
        MorrisAnalysis.Results.Std(CaEfflux_SOCE,k),'filled',...
        'MarkerFaceColor',	[1 0 0])
    scatter(MorrisAnalysis.Results.MuStar(Crossbridge_Cycle,k),...
        MorrisAnalysis.Results.Std(Crossbridge_Cycle,k),'filled',...
        'MarkerFaceColor',[0.4940 0.1840 0.5560])
    scatter(MorrisAnalysis.Results.MuStar(Diffusion,k),...
        MorrisAnalysis.Results.Std(Diffusion,k),'filled',...
        'MarkerFaceColor','g')
    xlabel('mu*')
    ylabel('sigma')
    title(expt_names{k})
    hold on
    ten_percent = max(MorrisAnalysis.Results.MuStar(:,k))*0.1;
    xline(ten_percent , '--b' )
    for i = 1:length(MorrisAnalysis.Results.MuStar(:,k))
        mu_star = MorrisAnalysis.Results.MuStar(i,k);
        if mu_star >= ten_percent
            muVec(i,k) = true;     
            sensIdx{qoiList==k} = [sensIdx{qoiList==k},i];
        else
            muVec(i,k)=false;
        end
    end
    for j= 1:length(muVec(:,k))
        shiftVal = 0.01*max(MorrisAnalysis.Results.MuStar(:,k));
        if muVec(j,k) ==1
            text(MorrisAnalysis.Results.MuStar(j,k)+shiftVal, ...
                MorrisAnalysis.Results.Std(j,k)-shiftVal, string(j));
        end
    end
    legend('Action Potential', 'Calcium Influx and SR Release',...
        'Calcium Buffering','Calcium Efflux and SOCE','Cross-Bridge Cycle',...
        'Diffusion','10% Threshold')
    prettyGraph
end

%% Test baseline model behavior - run baseline simulation for plots in Figs 3-4
load Data/p0Struct.mat p0Struct % load in default parameter values
p0 = p0Struct.data;
load Data/pSol_allParam.mat pSol
pSol(12) = 1; % Keep cref at default value
% VOnlyIdx is the list of parameters to which SL voltage is sensitive
VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
if length(pSol) == 28 % then only parameters used in Voltage fit
    highSensIdx = VOnlyIdx;
    VOnly = true;
elseif length(pSol) == 106 % then considering fit to both calcium and voltage data
    highSensIdx = 1:106;
else
    error('This length of pSol is not recognized')
end
% pSol is normalized to reference parameter values, so full parameter
% vector is given after multiplying by p0
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);
phosphateAccum = true; % whether phosphate is allowed to accumulate in the myoplasm and SR

% juncLocLogic and bulkLocLogic dictate whether a given state variable
% appears in junctional and/or bulk regions (see SkelMuscleCa_dydt)
juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions

% define geometric quantities
R_fiber = 20;
vol_SA_SR = 1/10;
vol_Fiber = pi * (R_fiber ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003;
volFraction_myo = 0.95;
volFraction_SR = 0.05;
vol_myo = volFraction_myo * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
diffusion_length = pPSO(99);
vol_myoJ = SA_TT*diffusion_length;
JFrac = vol_myoJ / vol_myo;
BFrac = 1 - JFrac;
SRJ_occupancy = 0.5;
geomParam = [vol_SA_ratio, volFraction_myo, volFraction_SR,...
             volFraction_TT, R_fiber, vol_SA_SR, SRJ_occupancy];

% run simulation for 0.5 s at 100 Hz
tSol = [0, 0.5];
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
SSPhosAccum = false; % whether phosphate accumulation is allowed in SS estimation (always false in tests for this paper)
[Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
    computeSol(pPSO, [], tSol, 1, 100, phosphateAccum, geomParam, SSPhosAccum);

%% Fig 3 plots
figure
subplot(4,1,1)
plot(Time_withSOCE, currents(:,13+13)/1000)
xlim([0 0.2])
ylabel('Applied current (nA)')
prettyGraph
subplot(4,1,2)
plot(Time_withSOCE, Y_withSOCE(:,5))
ylabel('SL voltage (mV)')
xlim([0 0.2])
prettyGraph
subplot(4,1,3)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
xlim([0 0.2])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(4,1,4)
plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
xlim([0 0.2])
ylabel('Rel force')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% Dissect currents for Fig 4A
% current in pA/s = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                %    I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% plot sodium, KDR, and KIR currents
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,sum(juncLocLogic) + sum(bulkLocLogic(1:5))))
xlim([0 .01])
ylabel('SL voltage (mV)')
prettyGraph
% sodium currents
subplot(3,1,2)
plot(Time_withSOCE, (currents(:,10)+currents(:,13+10))/1000)
hold on
plot(Time_withSOCE, (currents(:,7)+currents(:,13+7))/1000)
plot(Time_withSOCE, (currents(:,9)+currents(:,13+9))/1000)
legend('Na', 'NCX', 'NaK')
ylabel('Sodium current (nA)')
prettyGraph
xlim([0 .01])
% Potassium currents
subplot(3,1,3)
plot(Time_withSOCE, (currents(:,4)+currents(:,13+4))/1000)
hold on
plot(Time_withSOCE, (currents(:,5)+currents(:,13+5))/1000)
plot(Time_withSOCE, (currents(:,8)+currents(:,13+8))/1000)
legend('KDR', 'KIR', 'NaK')
xlim([0 .01])
ylabel('Potassium current (nA)')
xlabel('Time (s)')
prettyGraph

%% Dissect fluxes for Fig 4B
% note that plots are split between positive and negative fluxes to allow
% log y axis for each case (subplots 3 and 5 have flipped y axes)
% fluxes in µM/s = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
%                   LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];% 
figure
subplot(5,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
xlim([0 .01])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(5,1,2)
semilogy(Time_withSOCE, fluxes(:,6)*JFrac+fluxes(:,10+6)*BFrac)
hold on
semilogy(Time_withSOCE, fluxes(:,8)*JFrac+fluxes(:,10+8)*BFrac)
legend('RyR', 'Leak')
ylabel('SR fluxes (μM/s)')
prettyGraph
ylim([1e0 5e5])
xlim([0 .01])
subplot(5,1,3)
semilogy(Time_withSOCE, fluxes(:,7)*JFrac+fluxes(:,10+7)*BFrac)
xlim([0 .01])
ylim([1e0 5e5])
set(gca, 'ydir', 'reverse')
legend('SERCA')
ylabel('SR fluxes (μM/s)')
prettyGraph
subplot(5,1,4)
semilogy(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,10+1)*BFrac)
hold on
semilogy(Time_withSOCE, fluxes(:,2)*JFrac+fluxes(:,10+2)*BFrac)
legend('SOCE', 'SL leak')
xlim([0 .01])
ylim([.5e-2 .5e2])
ylabel('PM fluxes (μM/s)')
prettyGraph
subplot(5,1,5)
semilogy(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,10+5)*BFrac)
xlim([0 .01])
ylim([.5e-2 .5e2])
ylabel('PM fluxes (μM/s)')
set(gca, 'ydir', 'reverse')
legend('PMCA')
xlabel('Time (s)')
prettyGraph

%% Plot calcium in different compartments (Fig S3)
figure
subplot(1,3,1)
plot(Time_withSOCE, Y_withSOCE(:,8))
hold on
plot(Time_withSOCE, Y_withSOCE(:,32))
xlim([0 0.2])
legend('Junc','Bulk')
ylabel('Myo calcium (μM)')
prettyGraph
subplot(1,3,2)
plot(Time_withSOCE, Y_withSOCE(:,2))
hold on
plot(Time_withSOCE, Y_withSOCE(:,27))
ylabel('SR calcium (μM)')
xlabel('Time (s)')
xlim([0 0.2])
prettyGraph
subplot(1,3,3)
plot(Time_withSOCE, Y_withSOCE(:,22))
hold on
plot(Time_withSOCE, 1300*ones(size(Time_withSOCE)))
ylabel('EC calcium (μM)')
xlabel('Time (s)')
xlim([0 0.2])
prettyGraph

%% Plot crossbridge cycling for Fig S1
figure
subplot(4,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
hold on
ylabel('Myo calcium (μM)')
prettyGraph
subplot(4,1,2)
CaCaTrop = Y_withSOCE(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:18)));
CaTrop = Y_withSOCE(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:17)));
Apre = Y_withSOCE(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:20)));
Apost = Y_withSOCE(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:21)));
D = Y_withSOCE(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:19)));
Trop = pSol(73)*p0(73) - CaCaTrop - CaTrop - Apre - Apost - D;
plot(Time_withSOCE, CaTrop)
hold on
ylabel('CaTrop (μM)')
xlim([0 0.25])
prettyGraph
subplot(4,1,3)
plot(Time_withSOCE, CaCaTrop)
hold on
ylabel('CaCaTrop (μM)')
xlim([0 0.25])
prettyGraph
xlim([0 0.25])
subplot(4,1,4)
plot(Time_withSOCE, Apost)
hold on
ylabel('Apost (μM)')
prettyGraph
xlim([0 0.25])
xlabel('Time (s)')

%% SOCE vs no SOCE for Fig 5 TG test (Fig 5A)
% code allows for testing a range of SOCE conductance (soceFactor) and
% STIM1 calcium sensitivity (crefFactor), currently just keeps crefFactor
% at 1 for all tests
soceFactor = [0.3,1,3,10];
crefFactor = 1; % could be vector to test range
[soceFactor, crefFactor] = meshgrid(soceFactor, crefFactor);
ciSS = zeros(size(soceFactor));
cSRSS = zeros(size(soceFactor));
soceRef = pSol(95)*p0(95);
crefRef = pSol(12)*p0(12);
tSol = [0, 30*60];

% load in initial condition starting estimate
load Data/yinit0.mat yinit0
yinit0([24,26]) = pPSO(74:75);
ECVals = [1300;  % c_o (µM)
    147000.0;   % Na_o (µM)
    4000.0;      % K_o (µM)
    128000.0];  % Cl_o (µM)
yinit0(27:30) = ECVals;
juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
yinit = [yinit0(juncLocLogic); yinit0(bulkLocLogic)];

% define geometric quantities
vol_Fiber = pi * (20 ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003 ;
vol_myo = 0.95 * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
diffusion_length = pPSO(99);
vol_myoJ = SA_TT*diffusion_length;
JFrac = vol_myoJ / vol_myo;
BFrac = 1 - JFrac;
vol_SR = 0.05*vol_Fiber;
SRJ_occupancy = 0.5;
SA_SRJ = SA_TT * SRJ_occupancy;
vol_SRJ = SA_SRJ*diffusion_length;
JSRFrac = vol_SRJ / vol_SR;
BSRFrac = 1 - JSRFrac;

% compute steady state solution without SOCE
pPSO0 = pPSO;
pPSO0(95) = 0; % no SOCE
[TimeSS_noSOCE,ySS_noSOCE,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, false);%phosphateAccum);
yinf_noSOCE = ySSFinal_noSOCE;
[Time_noSOCE,Y_noSOCE,~,fluxes2,currents2] = SkelMuscleCa_dydt(tSol, 0, yinf_noSOCE, pPSO, tic, 8, phosphateAccum); % compute time-dependent solution
figure
subplot(2,1,1)
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
hold on
prettyGraph
ylabel('Myo calcium (μM)')
subplot(2,1,2)
plot(Time_noSOCE, Y_noSOCE(:,2)*JSRFrac + Y_noSOCE(:,27)*BSRFrac)
hold on
prettyGraph
ylabel('SR calcium (μM)')
xlabel('Time (s)')
finalCa = zeros(size(soceFactor));
finalCaSR = zeros(size(soceFactor));

for i = 1:size(soceFactor,1)
    for j = 1:size(soceFactor,2)
        pPSO(12) = crefFactor(i,j)*crefRef*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        pPSO(95) = soceFactor(i,j)*soceRef; % gSOCE
        pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
        [TimeSS_withSOCE,ySS_withSOCE,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, false);%phosphateAccum);
        yinf_withSOCE = ySSFinal_withSOCE;

        % test single frequency and plot calcium
        [Time_withSOCE,Y_withSOCE] = SkelMuscleCa_dydt(tSol, 0, yinf_withSOCE, pPSO, tic, 7, phosphateAccum); % compute time-dependent solution
        subplot(2,1,1)
        plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
        subplot(2,1,2)
        plot(Time_withSOCE, Y_withSOCE(:,2)*JSRFrac + Y_withSOCE(:,27)*BSRFrac)
        drawnow
        finalCa(i,j) = Y_withSOCE(end,8)*JFrac + Y_withSOCE(end,32)*BFrac;
        finalCaSR(i,j) = Y_withSOCE(end,2)*JSRFrac + Y_withSOCE(end,27)*BSRFrac;
        fprintf('Done with soceFactor=%.2f, crefFactor=%.2f: final ca = %.2f and final SR Ca = %.2f\r\n',...
            soceFactor(i,j), crefFactor(i,j), finalCa(i,j), finalCaSR(i,j))
    end
end

%% plot SOCE vs no SOCE over time for baseline conditions (Fig 5B-C)
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 20*pSol(95)*p0(95); % SOCE conductance
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
phosphateAccum = true; 

% test 1 s of stimulus here
tSol = [0, 1];
SSPhosAccum = false;
[Time_withSOCE,Y_withSOCE,~,~,maxCrossBridge, yinits] =...
    computeSol(pPSO, [], tSol, 1, 100, phosphateAccum, geomParam, SSPhosAccum);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[Time_noSOCE,Y_noSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100,...
    phosphateAccum, geomParam, SSPhosAccum);

figure
subplot(3,1,1)
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
hold on
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
% xlim([0 2])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(3,1,2)
plot(Time_noSOCE, Y_noSOCE(:,2)*JSRFrac + Y_noSOCE(:,27)*BSRFrac)
hold on
plot(Time_withSOCE, Y_withSOCE(:,2)*JSRFrac + Y_withSOCE(:,27)*BSRFrac)
% xlim([0 2])
ylabel('SR calcium (μM)')
prettyGraph
subplot(3,1,3)
plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge)
hold on
plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
% xlim([0 2])
ylabel('Rel force')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% plot phosphate accumulation and force w/o phophate accum (Fig S4)
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 1*pSol(95)*p0(95); % gSOCE
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
tSol = [0, 2];
[tWithPhos,yWithPhos,fluxes,currents,maxCrossBridge,yinits] = computeSol(pPSO, [], tSol, 1, 100, true);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[tNoPhos,yNoPhos] = computeSol(pPSO, yinits{2}, tSol, 1, 100, false);
figure
subplot(4,1,1)
plot(tWithPhos, yWithPhos(:,crossbridgeIdx)/maxCrossBridge)
hold on
plot(tNoPhos, yNoPhos(:,crossbridgeIdx)/maxCrossBridge)
prettyGraph
subplot(4,1,2)
plot(tWithPhos, yWithPhos(:,21)*JFrac + yWithPhos(:,50)*BFrac)
hold on
plot(tNoPhos, yNoPhos(:,21)*JFrac + yNoPhos(:,50)*BFrac)
prettyGraph
subplot(4,1,3)
plot(tWithPhos, yWithPhos(:,19)*JFrac + yWithPhos(:,48)*BFrac)
hold on
plot(tNoPhos, yNoPhos(:,19)*JFrac + yNoPhos(:,48)*BFrac)
prettyGraph
subplot(4,1,4)
plot(tWithPhos, yWithPhos(:,20)*JFrac + yWithPhos(:,49)*BFrac)
hold on
plot(tNoPhos, yNoPhos(:,20)*JFrac + yNoPhos(:,49)*BFrac)
prettyGraph

%% plot hydrolysis rates (Fig S4)
F = 96485.3321;
Jhyd = fluxes(:,9)*JFrac + fluxes(:,19)*BFrac;
J_NKX_ATP = 1e9*((currents(:,8) + currents(:,13+8)))/(2*F*vol_myo);  % one ATP per two potassium pumped, uM/s
J_SERCA_ATP = (fluxes(:,7)*JFrac+fluxes(:,10+7)*BFrac)/2; % uM/s, 1 ATP transports 2 calcium ions
J_PMCA_ATP = fluxes(:,5)*JFrac+fluxes(:,10+5)*BFrac; % uM/s
J_XB_ATP = pPSO(63)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:19)))*BFrac - ...
    pPSO(64)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:20)))*BFrac;
% J_XB_ATP = pPSO(67)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:21)))*BFrac;
ATPjunc = yWithPhos(:,sum(juncLocLogic(1:23)));
ATPbulk = yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:23)));
J_baselinehyd = (pPSO(55)*(ATPjunc./(pPSO(56)+ATPjunc))*JFrac +...
    pPSO(55)*(ATPbulk./(pPSO(56)+ATPbulk))*BFrac);

figure
plot(tWithPhos, Jhyd)
hold on
plot(tWithPhos, J_NKX_ATP)
plot(tWithPhos, J_SERCA_ATP)
plot(tWithPhos, J_PMCA_ATP)
plot(tWithPhos, J_XB_ATP)
plot(tWithPhos, J_baselinehyd)
xlim([0 0.5])
prettyGraph
xlabel('Time (s)')
ylabel('Hydrolysis rate (μM/s)')
legend('tot','NKX','SERCA','PMCA','XB','baseline')

%% plot calcium with vs without SOCE for different initial conditions
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 20*pSol(95)*p0(95); % gSOCE
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
tSol = [0, 1];
[tWithSOCE,yWithSOCE,~,~,maxCrossBridge,yinits] = computeSol(pPSO, [], tSol, 1, 100, true);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[tNoSOCE,yNoSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100, true);
[tNoSOCE_altInit,yNoSOCE_altInit] = computeSol(pPSO, yinits{2}, tSol, 2, 100, true);
figure
subplot(2,1,1)
[t1,y1] = getMaxes(tNoSOCE, yNoSOCE(:,8)*JFrac + yNoSOCE(:,32)*BFrac, 100);
% plot(t1,y1)
plot(tNoSOCE, yNoSOCE(:,8)*JFrac + yNoSOCE(:,32)*BFrac)
hold on
[t2,y2] = getMaxes(tWithSOCE, yWithSOCE(:,8)*JFrac + yWithSOCE(:,32)*BFrac, 100);
% plot(t2,y2)
plot(tWithSOCE, yWithSOCE(:,8)*JFrac + yWithSOCE(:,32)*BFrac)
[t3,y3] = getMaxes(tNoSOCE_altInit, yNoSOCE_altInit(:,8)*JFrac + yNoSOCE_altInit(:,32)*BFrac, 100);
% plot(t3, y3)
plot(tNoSOCE_altInit, yNoSOCE_altInit(:,8)*JFrac + yNoSOCE_altInit(:,32)*BFrac)
ylabel('Myo calcium (μM)')
prettyGraph
subplot(2,1,2)
[t4,y4] = getMaxes(tNoSOCE, yNoSOCE(:,2)*JSRFrac + yNoSOCE(:,27)*BSRFrac, 100);
% plot(t4, y4)
plot(tNoSOCE, yNoSOCE(:,2)*JSRFrac + yNoSOCE(:,27)*BSRFrac)
hold on
[t5,y5] = getMaxes(tWithSOCE, yWithSOCE(:,2)*JSRFrac + yWithSOCE(:,27)*BSRFrac, 100);
% plot(t5,y5)
plot(tWithSOCE, yWithSOCE(:,2)*JSRFrac + yWithSOCE(:,27)*BSRFrac)
[t6,y6] = getMaxes(tNoSOCE_altInit, yNoSOCE_altInit(:,2)*JSRFrac + yNoSOCE_altInit(:,27)*BSRFrac, 100);
% plot(t6,y6)
plot(tNoSOCE_altInit, yNoSOCE_altInit(:,2)*JSRFrac + yNoSOCE_altInit(:,27)*BSRFrac)
ylabel('SR calcium (μM)')
prettyGraph

%% test range of frequencies for SOCE vs no SOCE - Fig 6
tSol = [0, 0.5];
freq = 5:5:200;
inclAltPhos = true; % whether to include additional test with 150% resting phosphate in SOCE KOs
peakForce = zeros(2+inclAltPhos, length(freq));
endPeakRatio = zeros(2+inclAltPhos, length(freq));
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 20*pSol(95)*p0(95); % gSOCE
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
if inclAltPhos
    pPSOAlt = pPSO;
    pPSOAlt(75) = 1.5*pSol(75)*p0(75); % myo phosphate
end
phosphateAccum = true;
SSPhosAccum = false;
for i = 1:length(freq)
    % expt 1: HIIT stim, with SOCE, expt 2: HIIT stim, no SOCE
    if i == 1
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
            computeSol(pPSO, [], tSol, 1, freq(i), phosphateAccum, geomParam, SSPhosAccum);
        pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
        [Time_noSOCE,Y_noSOCE,~,~,~,yinitNoSOCE] =...
            computeSol(pPSO, [], tSol, 2, freq(i), phosphateAccum, geomParam, SSPhosAccum);
        if inclAltPhos
            [~,Y_noSOCEAlt,~,~,~,yinitNoSOCEAlt] =...
                computeSol(pPSOAlt, [], tSol, 2, freq(i), phosphateAccum, geomParam, SSPhosAccum);
            pPSOAlt(12) = pPSOAlt(12)*yinitNoSOCEAlt{1}(2); % define cref in terms of resting SR calcium without SOCE
        end
    else
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
            computeSol(pPSO, yinits{2}, tSol, 1, freq(i), phosphateAccum, geomParam, SSPhosAccum);
        [Time_noSOCE,Y_noSOCE] =...
            computeSol(pPSO, yinitNoSOCE{1}, tSol, 2, freq(i), phosphateAccum, geomParam, SSPhosAccum);
        if inclAltPhos
            [~,Y_noSOCEAlt] =...
                computeSol(pPSOAlt, yinitNoSOCEAlt{1}, tSol, 2, freq(i), phosphateAccum, geomParam, SSPhosAccum);
        end
    end
   
    % create plots
    if any(freq(i)==[20, 80])
        figure
        colors = colororder;
        subplot(3,1,1)
        % final entry in color spec gives the opacity
        plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
        ylabel('Myo calcium (µM)')
        % title('Cytosolic calcium')
        legend('No SOCE', '20x SOCE')
        prettyGraph
        subplot(3,1,2)
        plot(Time_noSOCE, Y_noSOCE(:,2))
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,2))
        ylabel('SR calcium (µM)')
        prettyGraph
        % title('SR calcium')
        subplot(3,1,3)
        plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge)
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
        xlabel('Time (s)')
        ylabel('Rel force')
        prettyGraph
        drawnow
    end
    peakForce(1,i) = max(Y_withSOCE(:,crossbridgeIdx));
    peakForce(2,i) = max(Y_noSOCE(:,crossbridgeIdx));
    endPeakRatio(1,i) = Y_withSOCE(end,crossbridgeIdx)/peakForce(1,i);
    endPeakRatio(2,i) = Y_noSOCE(end,crossbridgeIdx)/peakForce(2,i);
    if inclAltPhos
        peakForce(3,i) = max(Y_noSOCEAlt(:,crossbridgeIdx));
        endPeakRatio(3,i) = Y_noSOCEAlt(end,crossbridgeIdx)/peakForce(3,i);
    end
    fprintf('Done with freq = %.2f Hz\n', freq(i))
end

%% plot and compare freq-dependent behavior with exp
freqExp = [1;25;50;75;100;125;150;175;200;250];
maxForceExpMean_withSOCE = [0.239024373; 0.321951196; 0.595121908; 0.824390184; 0.951219443;...
                            0.999999927; 0.999999927; 0.975609685; 0.946341394; 0.907317007];
maxForceExpSEM_withSOCE = [0.009756097; 0.014634145; 0.019512194; 0.024390242; 0.024390242;...
                           0.024390242; 0.024390242; 0.024390242; 0.029268291; 0.029268291];
forceRatioExpMean_withSOCE = [0.972839506; 0.992592593; 1.002469136; 1.002469136;...
                               0.997530864; 0.982716049; 0.95308642; 0.92345679];
forceRatioExpSEM_withSOCE = [0.009876543; 0.009876543; 0.004938272; 0.004938272;...
                             0.004938272; 0.004938272; 0.009876543; 0.019753086];
maxForceExpMean_noSOCE =  [0.229268276; 0.297560954; 0.546341423; 0.751219457; 0.824390184;...
                           0.848780426; 0.843902377; 0.829268232; 0.80487799; 0.760975554];
maxForceExpSEM_noSOCE = [0.009756097; 0.014634145; 0.024390242; 0.024390242; 0.019512194;...
                         0.019512194; 0.024390242; 0.019512194; 0.019512194; 0.019512194];
forceRatioExpMean_noSOCE = [0.972839506; 0.967901235; 0.958024691; 0.962962963;...
                            0.958024691; 0.928395062; 0.854320988; 0.725925926];
forceRatioExpSEM_noSOCE = [0.004938272; 0.004938272; 0.009876543; 0.009876543;...
                           0.009876543; 0.009876543; 0.014814815; 0.014814815];
% force = 1 in exp corresponds to force of 169.4 mN/mm^2
figure
subplot(1,2,1)
plot(freq, peakForce/maxCrossBridge)
prettyGraph
ylabel('Rel force')
xlabel('Freq (Hz)')
subplot(1,2,2)
errorbar(freqExp, maxForceExpMean_withSOCE*169.4, maxForceExpSEM_withSOCE*169.4)
hold on
errorbar(freqExp, maxForceExpMean_noSOCE*169.4, maxForceExpSEM_noSOCE*169.4)
ylim([25 200])
prettyGraph
ylabel('Specific force (mN/mm2)')
xlabel('Freq (Hz)')

%% plot SOCE vs no SOCE over time for fatiguing conditions (Fig 7)
exercise = 'HIIT'; % set this to either 'HIIT' or 'resistance'
pPSO(12) = 0.5; % set cref ratio
soceVals = 10.^[-.5,.5,1.5]*pSol(95)*p0(95); % gSOCE
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
phosphateAccum = true;
pPSO(59) = 1*pSol(59)*p0(59); % kon1 (trop binding)
SSPhosAccum = false;
tCell = cell(length(soceVals),1);
yCell = cell(length(soceVals),1);
switch exercise
    case 'HIIT'
        exerciseParam = [0.65, 0.25*0.65]; % sprinting
        freq = 100;
        tSol = [0, 20];
    case 'resistance'
        exerciseParam = [6, 3]; % resistance
        freq = 40;
        tSol = [0, 60];
    otherwise
        error('exercise "%s" not recognized, must be either "HIIT" or "resistance"',exercise)
end
figure
for i = 1:length(tCell)
    pPSO(95) = soceVals(i);
    [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
        computeSol(pPSO, [], tSol, [3, exerciseParam], freq, phosphateAccum, geomParam, SSPhosAccum);
    tCell{i} = Time_withSOCE;
    yCell{i} = Y_withSOCE;
    subplot(5,1,1)
    plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
    hold on
    ylabel('Myo calcium (μM)')
    prettyGraph
    subplot(5,1,2)
    plot(Time_withSOCE, Y_withSOCE(:,2)*JSRFrac + Y_withSOCE(:,27)*BSRFrac)
    hold on
    ylabel('SR calcium (μM)')
    prettyGraph
    subplot(5,1,3)
    plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
    hold on
    ylabel('Rel force')
    prettyGraph
    subplot(5,1,4)
    plot(Time_withSOCE, Y_withSOCE(:,21)*JFrac + Y_withSOCE(:,50)*BFrac)
    hold on
    ylabel('Myo phosphate (μM)')
    prettyGraph
    subplot(5,1,5)
    plot(Time_withSOCE, Y_withSOCE(:,20)*JSRFrac + Y_withSOCE(:,49)*BSRFrac)
    hold on
    ylabel('SR calcium phosphate (μM)')
    xlabel('Time (s)')
    prettyGraph
    drawnow
end
set(gcf,'Renderer','painters')

%% plot fatiguing dynamics over time for Fig 7
tLims = 0:exerciseParam(1):tSol(end);
maxForces = zeros(length(tCell),length(tLims));
maxCaCaTrop = zeros(length(tCell),length(tLims));
figure
for i = 1:length(tCell)
    CaCaTropVec = yCell{i}(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:18)));
    subplot(3,1,1)
    plot(tCell{i},CaCaTropVec)
    hold on
    % xlim([0 20])
    prettyGraph
    ylabel('CaCaTrop (μM)')
    phosVec = yCell{i}(:,21)*JFrac + yCell{i}(:,50)*BFrac;
    subplot(3,1,2)
    plot(tCell{i},phosVec)
    hold on
    % xlim([0 20])
    prettyGraph
    ylabel('Myo phosphate (μM)')
    forceVec = yCell{i}(:,crossbridgeIdx);
    subplot(3,1,3)
    plot(tCell{i},forceVec/maxCrossBridge)
    hold on
    % xlim([0 20])
    prettyGraph
    ylabel('Rel force')
    for j = 2:length(tLims)
        curLogic = tCell{i} > tLims(j-1) & tCell{i} <= tLims(j);
        maxForces(i,j-1) = max(forceVec(curLogic)/maxCrossBridge);
        maxCaCaTrop(i,j-1) = max(CaCaTropVec(curLogic));
    end
    maxForces(i,end) = maxForces(i,end-1);
    maxCaCaTrop(i,end) = maxCaCaTrop(i,end-1);
end

%% Phase diagrams for Fig 7
exerciseNames = {'RUNNING', 'RESISTANCE'};
freqVals = [100, 40];
for k = 1:length(exerciseNames)
    exercise = exerciseNames{k};
    figure
    turboMap = colormap('turbo');
    colormap(turboMap(10:240,:))
    curFile = sprintf('Data/sweepSOCE_1.0phosdeg_1.0Trop_%s.mat',exercise);
    load(curFile,'results','freq','soceFactor','crefFactor','tLims')
    maxCrossBridge = results{freq==freqVals(k)}{1};
    dynCell = results{freq==freqVals(k)}{end};
    finalForceVals = zeros(size(dynCell));
    if contains(exercise,'RUNNING')
        testIdx = find(tLims<20,1,'last');
    else
        testIdx = length(dynCell{1});
    end
    for idx = 1:length(finalForceVals(:))
        finalForceVals(idx) = mean(dynCell{idx}(testIdx));
    end
    surf(log10(soceFactor), log10(crefFactor), finalForceVals./maxCrossBridge, 'FaceColor', 'interp')
    xlabel('SOCE flux')
    ylabel('cref')
    colorbar
    view([0 90])
    title(exercise)
    set(gcf,'Renderer','painters')
end

%% Plot array of phase diagrams for Figs S5-S6
exerciseNames = {'RUNNING', 'RESISTANCE'};
freqVals = [100, 40];
konFactor = [0.1, 0.1, 1, 1];
phosDegFactor = [1, 5, 1, 5];
for k = 1:length(exerciseNames)
    exercise = exerciseNames{k};
    figure
    turboMap = colormap('turbo');
    colormap(turboMap(10:240,:))
    for i = 1:length(konFactor)
        subplot(2,2,i)
        curFile = sprintf('Data/sweepSOCE_%.1fphosdeg_%.1fTrop_%s.mat',...
            phosDegFactor(i),konFactor(i),exercise);
        load(curFile,'results','freq','soceFactor','crefFactor','tLims')
        dynCell = results{freq==freqVals(k)}{end};
        maxCrossBridge = results{freq==freqVals(k)}{1};
        finalForceVals = zeros(size(dynCell));
        if contains(exercise,'RUNNING')
            testIdx = find(tLims<20,1,'last');
        else
            testIdx = length(dynCell{1});
        end
        for idx = 1:length(finalForceVals(:))
            finalForceVals(idx) = mean(dynCell{idx}(testIdx));
        end
        surf(log10(soceFactor), log10(crefFactor), finalForceVals./maxCrossBridge, 'FaceColor', 'interp')
        xlabel('SOCE flux')
        ylabel('cref')
        if contains(exercise,'RUNNING') && phosDegFactor(i)==1
            clim([0.115 0.17])
        elseif contains(exercise,'RUNNING') && phosDegFactor(i)==5
            clim([0.26 0.38])
        elseif contains(exercise,'RESISTANCE') && phosDegFactor(i)==1
            clim([0.015 0.105])
        elseif contains(exercise,'RESISTANCE') && phosDegFactor(i)==5
            clim([0.07 0.36])
        end
        colorbar
        view([0 90])
        title(sprintf('phosDeg=%.1f, kon=%.1f, %s', phosDegFactor(i),konFactor(i),exercise))
    end
    set(gcf,'Renderer','painters')
end

function [t,y,fluxes,currents,maxCrossBridge,yinits] = ...
    computeSol(pPSO, yinit, tSol, expt, freq, phosphateAccum, varargin)
    % Utility function to compute solution for different conditions and
    % parameter sets    
    % Inputs:
    %   - pPSO: vector of parameters (not normalized)
    %   - yinit: vector of initial values for all 51 variables (if
    %   previously solved for), otherwise empty and solved for here
    %   - tSol: vector with start time and end time for simulation
    %   - expt: Scalar with the expt number or: 
    %           [expt_n, cycle_time, stim_time], where expt_n is 3 or 4,
    %           cycle_time is the total time per repitition/stride and stim_time
    %           is the time of activation for each repitition/stride.
    %   - freq: Freqency of stimulus in Hz
    %   - phosphateAccum: logical variable determining if phosphate
    %       accumulation is accounted for
    % Optional inputs (stored in varargin):
    %   - geomParam (varargin{1}): vector with stored geometric parameters.
    %       If not provided, geomParam is assigned default values (see code
    %       below)
    %   - SSPhosAccum (varargin{2}): whether phosphate is allowed to
    %       accumulate at steady state (default is false)
    % Outputs:
    %   - t is the vector of times
    %   - y is the vector of state variables
    %   - fluxes: calcium fluxes at each time point, each row consists:
    %        [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, 
    %         LumpedJ_SERCA, J_CaLeak_SR, Jhydtot, total calcium flux]
    %        in units of µM/s for junctional myoplasm, then the same 
    %        entities for the bulk myoplasm
    %   - currents: Ionic and total current at each time point, each row consists of:
    %        [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N,
    %         I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL] in units of pA
    %         for the T-tubules and then the same entities for the SL
    %   - maxCrossBridge: post power stroke cross bridge density at
    %     saturating calcium concentration for these parameters
    %   - yinits: Cell vector with each entry a vector containing the initial values
    %     for all 51 state variables in the model.
    %     Cell length 1 if only SOCE-free condition is tested, or if
    %     initial condition is provided.
    %     Otherwise first cell entry is for case without SOCE and second
    %     entry is for case with SOCE.

    if isscalar(varargin) % length 1, that is
        geomParam = varargin{1};
        SSPhosAccum = false;
    elseif length(varargin)==2
        geomParam = varargin{1};
        SSPhosAccum = varargin{2};
    elseif isempty(varargin)
        geomParam = [0.01, 0.95, 0.05, 0.003, 20, 1/10, 0.5];
        SSPhosAccum = false;
    else
        error('Unknown argument(s)')
    end

    juncLocLogic = true(1,31);
    juncLocLogic(17:21) = false; % cross bridges
    bulkLocLogic = true(1,31);
    bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
    if isempty(yinit)
        % load in initial condition starting estimate
        load Data/yinit0.mat yinit0
        yinit0([24,26]) = pPSO(74:75);
        if yinit0(31) > pPSO(96)
            yinit0(31) = pPSO(96);
        end
        yinit = [yinit0(juncLocLogic); yinit0(bulkLocLogic)];
        
        % compute steady state solution without SOCE
        pPSO0 = pPSO;
        pPSO0(95) = 0; % no SOCE
        [~,~,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit,...
            pPSO0, tic, 2, SSPhosAccum, geomParam);
        pPSO(12) = pPSO(12)*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        if any(expt(1) == [2,4]) || pPSO(95) == 0
            yinf = ySSFinal_noSOCE;
            yinits = {yinf};
        else
            [~,~,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit,...
                pPSO, tic, 1, SSPhosAccum, geomParam);
            yinf = ySSFinal_withSOCE;
            yinits = {ySSFinal_noSOCE, ySSFinal_withSOCE};
        end
    else
        yinf = yinit;
        yinits = {yinit};
        if pPSO(12) < yinit(2)/10
            warning('cref may not have been initialized correctly')
        end
    end
    
    % compute max cross bridge engagement - expt 10 is only crossbridge testing
    % (all other dydt = 0)
    yinit_maxCa = yinf;
    yinit_maxCa([21,50]) = pPSO(75); % assuming no phos accum
    yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
    yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
    [~, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa,...
        pPSO, tic, 10, phosphateAccum, geomParam);
    crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
    maxCrossBridge = Y_maxCa(end,crossbridgeIdx);
    
    % compute time-dependent solution
    [t,y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, yinf, pPSO,...
        tic, expt, phosphateAccum, geomParam);
end

function [tMax,yMax] = getMaxes(t, y, freq)
    % pull out maxes over each time interval specified by 1/freq
    % Inputs:
    %   - t: time vector
    %   - y: data vector
    %   - freq: frequency of characteristic events (time intervals sliced
    %   into chunks of width 1/freq)
    %
    % Outputs:
    %   - tMax: vector of t values at edges of each interval. Final
    %   interval is defined as [tMax(end), t(end)]
    %   - yMax: maximum y values over each interval, where yMax(i) is
    %   assessed over the range of times [tMax(i), tMax(i+1)]

    tMax = min(t):1/freq:max(t);
    yMax = zeros(size(tMax));
    if max(t) > tMax(end)+0.1/freq
        tMax = [tMax, max(t)];
    else
        yMax = yMax(1:end-1);
    end
    for i = 1:length(yMax)
        curLogic = t > tMax(i) & t <= tMax(i+1);
        yMax(i) = max(y(curLogic));
    end
    tMax = tMax(1:end-1);
end