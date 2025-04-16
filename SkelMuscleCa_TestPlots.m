%% Test SOCE vs no SOCE steady states after adjusting parameters
load p0Struct.mat p0Struct
p0 = p0Struct.data;
load PSOpenaltyEst3_07-Apr-2025.mat pSol
VOnlyIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,...
            28,30,32,33,40,76,77,78,79,80,81,82,83,91];
if length(pSol) == 35 % then VOnly
    highSensIdx = VOnlyIdx;
    VOnly = true;
else % then fitting to both calcium and V using ca sens indices
    highSensIdx = 1:105;
    VOnlyStruct = load('pVec_VOnly.mat', 'pVec');
    pVecVOnly = VOnlyStruct.pVec;
    pVecVOnly = pVecVOnly(:); % be sure it is a column vector
    [VOnlyOnly,onlyIdx] = setdiff(VOnlyIdx, highSensIdx); % non overlapping indices
    p0(VOnlyOnly) = pVecVOnly(onlyIdx).*p0(VOnlyOnly); % set p0 according to previous estimation
end
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);

phosphateAccum = true;
% load in initial condition starting estimate
load yinit0.mat yinit0
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

% compute steady state solution without SOCE
pPSO0 = pPSO;
pPSO0(95) = 0; % no SOCE
[TimeSS_noSOCE,ySS_noSOCE,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, phosphateAccum);
pPSO(12) = pPSO(12)*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
[TimeSS_withSOCE,ySS_withSOCE,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, phosphateAccum);
yinf_noSOCE = ySSFinal_noSOCE;
yinf_withSOCE = ySSFinal_withSOCE;

% compute max cross bridge engagement - expt 10 is only crossbridge testing
% (all other dydt = 0)
yinit_maxCa = ySS_withSOCE(end,:);
yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
[~, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa, pPSO, tic, 10, phosphateAccum);
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
maxCrossBridge = Y_maxCa(end,crossbridgeIdx);

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

% test single frequency dynamics
tSol = [0, 0.5];
[Time_withSOCE,Y_withSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, 100, yinf_withSOCE, pPSO, tic, 1, phosphateAccum); % compute time-dependent solution

%% Fig 3A plots
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

%% Dissect currents for Fig 3B
% current in pA/s = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                %    I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% plot sodium, KDR, and KIR currents
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,5))
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

%% Dissect fluxes for Fig 3C
% fluxes in µM/s = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
%                   LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];% 
figure
subplot(5,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
xlim([0 .01])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(5,1,2)
semilogy(Time_withSOCE, fluxes(:,6)*JFrac+fluxes(:,8+6)*BFrac)
hold on
semilogy(Time_withSOCE, fluxes(:,8)*JFrac+fluxes(:,8+8)*BFrac)
legend('RyR', 'Leak')
ylabel('SR fluxes (μM/s)')
prettyGraph
ylim([1e0 5e5])
xlim([0 .01])
subplot(5,1,3)
semilogy(Time_withSOCE, fluxes(:,7)*JFrac+fluxes(:,8+7)*BFrac)
xlim([0 .01])
ylim([1e0 5e5])
set(gca, 'ydir', 'reverse')
legend('SERCA')
ylabel('SR fluxes (μM/s)')
prettyGraph
subplot(5,1,4)
semilogy(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,8+1)*BFrac)
hold on
semilogy(Time_withSOCE, fluxes(:,2)*JFrac+fluxes(:,8+2)*BFrac)
legend('SOCE', 'SL leak')
xlim([0 .01])
ylim([.5e-2 .5e2])
ylabel('PM fluxes (μM/s)')
prettyGraph
subplot(5,1,5)
semilogy(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,8+5)*BFrac)
xlim([0 .01])
ylim([.5e-2 .5e2])
ylabel('PM fluxes (μM/s)')
set(gca, 'ydir', 'reverse')
legend('PMCA')
xlabel('Time (s)')
prettyGraph

%% Dissect SOCE in baseline case
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
ylabel('Myo calcium (μM)')
prettyGraph
subplot(3,1,2)
plot(Time_withSOCE, Y_withSOCE(:,2)*JSRFrac + Y_withSOCE(:,27)*BSRFrac)
ylabel('SR calcium (μM)')
yyaxis right
plot(Time_withSOCE, Y_withSOCE(:,1))
ylabel('SOCE open probability')
prettyGraph
subplot(3,1,3)
plot(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,8+1)*BFrac)
ylabel('SOCE flux (μM/s)')
xlabel('Time (s)')
prettyGraph

%% Plot calcium in different compartments
figure
subplot(2,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8))
hold on
plot(Time_withSOCE, Y_withSOCE(:,32))
xlim([0 0.2])
legend('Junc','Bulk')
ylabel('Myo calcium (μM)')
prettyGraph
subplot(2,1,2)
plot(Time_withSOCE, Y_withSOCE(:,2))
hold on
plot(Time_withSOCE, Y_withSOCE(:,27))
legend('Junc','Bulk')
ylabel('SR calcium (μM)')
xlabel('Time (s)')
xlim([0 0.2])
prettyGraph

%% SOCE vs no SOCE for figure 4A thaps test
soceFactor = [0.5,1,2];%logspace(-1, 1, 6);
crefFactor = 1;%logspace(-1, 1, 6);
[soceFactor, crefFactor] = meshgrid(soceFactor, crefFactor);
ciSS = zeros(size(soceFactor));
cSRSS = zeros(size(soceFactor));
soceRef = pSol(95)*p0(95);
crefRef = pSol(12)*p0(12);
tSol = [0, 30*60];

% pPSO(41) = 0.01*pSol(41)*p0(41); % tau SOCE
% pPSO(96) = 0.25*pSol(96)*p0(96); % total SR buffer
% pPSO(89) = 0*pSol(89)*p0(89); % gSL leak Ca

% load in initial condition starting estimate
load yinit0.mat yinit0
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

%% plot SOCE vs no SOCE over time for baseline conditions (Fig 4)
% CaSensIdx = [2,3,4,6,8,14,15,20,22,23,26,28,29,32,35,38,40,42,43,69,72,77,78,79,80,81,82,83,86,90,91];
% pPSO = params{86};
% pPSO(42) = 1.0*pSol(42)*p0(42); % SERCA
% pPSO(43) = 1*pSol(43)*p0(43); % PMCA
% pPSO(96) = 1.0*pSol(96)*p0(96); % total SR buffer
pPSO(12) = 1*pSol(12)*p0(12); % set cref ratio
pPSO(95) = 2*pSol(95)*p0(95); % gSOCE
pPSO(92) = 1*pSol(92)*p0(92); % j0 RyR
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
% 59,60,93,94
pPSO(59) = 0.2*pSol(59)*p0(59);
% pPSO(60) =
% pPSO(93) = 0.5*pSol(93)*p0(93);
% pPSO(94) =
% pPSO(41) = 1.0*pSol(41)*p0(41); % tau SOCE
% pPSO(74) = pSol(74)*p0(74); %5000; % SR phosphate
% pPSO(75) = 20e3; % myo phosphate
% pPSO(89) = 0*pSol(89)*p0(89); % gSL leak Ca
tSol = [0, 2];
[Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
    computeSol(pPSO, [], tSol, 1, 100, phosphateAccum);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
% pPSO(75) = 0; % myo phosphate
[Time_noSOCE,Y_noSOCE,fluxes] = computeSol(pPSO, yinits{1}, tSol, 2, 100, phosphateAccum);
figure
subplot(3,1,1)
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
hold on
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
% xlim([0 2])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(3,1,2)
plot(Time_noSOCE, Y_noSOCE(:,2))
hold on
plot(Time_withSOCE, Y_withSOCE(:,2))
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

%% test range of frequencies for SOCE vs no SOCE - Fig 5
tSol = [0, 0.5]; % 420 is max time for HIIT
freq = logspace(1,2.2,50);%[1, 25, 50, 75, 100, 125];%, 150];%, 175, 200, 250];
forceRatio = zeros(size(freq));
peakForce = zeros(2, length(freq));
endPeakRatio = zeros(2, length(freq));
% pPSO(96) = 0.25*pSol(96)*p0(96); % total SR buffer
pPSO(12) = 1*pSol(12)*p0(12); % set cref ratio
pPSO(95) = 10*pSol(95)*p0(95); % gSOCE
pPSO(41) = 1*pSol(41)*p0(41); % tau SOCE
phosphateAccum = true;
figure
for i = 1:length(freq)
    % figure
    % expt 1: HIIT stim, with SOCE, expt 2: HIIT stim, no SOCE
    if i == 1
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
            computeSol(pPSO, [], tSol, 1, freq(i), phosphateAccum);
        pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
        [Time_noSOCE,Y_noSOCE,~,~,~,yinitNoSOCE] =...
            computeSol(pPSO, yinits{1}, tSol, 2, freq(i), phosphateAccum);
    else
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
            computeSol(pPSO, yinits{2}, tSol, 1, freq(i), phosphateAccum);
        [Time_noSOCE,Y_noSOCE,~,~,~,yinitNoSOCE] =...
            computeSol(pPSO, yinits{1}, tSol, 2, freq(i), phosphateAccum);
    end
   
    % create plots
    colors = colororder;
    subplot(3,1,1)
    % final entry in color spec gives the opacity
    plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac,...
        'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac,...
        'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('Cytosolic calcium conc. (µM)')
    % title('Cytosolic calcium')
    legend('With SOCE', 'No SOCE')
    prettyGraph
    subplot(3,1,2)
    plot(Time_withSOCE, Y_withSOCE(:,2), 'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,2), 'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('SR calcium conc. (µM)')
    prettyGraph
    % title('SR calcium')
    subplot(3,1,3)
    plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge, 'Color', [colors(1,:), 0.8])
    hold on
    plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge, 'Color', [colors(2,:), 0.5])
    xlabel('Time (seconds)')
    ylabel('Relative force')
    prettyGraph
    % title('Force')
    drawnow
    forceRatio(i) = max(Y_withSOCE(:,crossbridgeIdx))/max(Y_noSOCE(:,crossbridgeIdx));
    peakForce(1,i) = max(Y_withSOCE(:,crossbridgeIdx));
    peakForce(2,i) = max(Y_noSOCE(:,crossbridgeIdx));
    endPeakRatio(1,i) = Y_withSOCE(end,crossbridgeIdx)/peakForce(1,i);
    endPeakRatio(2,i) = Y_noSOCE(end,crossbridgeIdx)/peakForce(2,i);
end
peakForce = peakForce / max(peakForce(:));
figure
semilogx(freq, peakForce)

function [t,y,fluxes,currents,maxCrossBridge,yinits] = computeSol(pPSO, yinit, tSol, expt, freq, phosphateAccum)
    juncLocLogic = true(1,31);
    juncLocLogic(17:21) = false; % cross bridges
    bulkLocLogic = true(1,31);
    bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
    if isempty(yinit)
        % load in initial condition starting estimate
        load yinit0.mat yinit0
        yinit0([24,26]) = pPSO(74:75);
        if yinit0(31) > pPSO(96)
            yinit0(31) = pPSO(96);
        end
        yinit = [yinit0(juncLocLogic); yinit0(bulkLocLogic)];
        
        % compute steady state solution without SOCE
        pPSO0 = pPSO;
        pPSO0(95) = 0; % no SOCE
        [~,~,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, phosphateAccum);
        pPSO(12) = pPSO(12)*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        if any(expt == [2,8]) || pPSO(95) == 0
            yinf = ySSFinal_noSOCE;
            yinits = {yinf};
        else
            [~,~,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, phosphateAccum);
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
    yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
    yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
    [~, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa, pPSO, tic, 10, phosphateAccum);
    crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
    maxCrossBridge = Y_maxCa(end,crossbridgeIdx);
    
    % test single frequency dynamics
    [t,y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, yinf, pPSO, tic, expt, phosphateAccum); % compute time-dependent solution
end