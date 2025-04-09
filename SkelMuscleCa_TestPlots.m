%% Test SOCE vs no SOCE steady states after adjusting parameters
load p0Struct.mat p0Struct
p0 = p0Struct.data;
load PSOpenaltyEst3_07-Apr-2025.mat pSol
VOnlyIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,...
            28,30,32,33,40,76,77,78,79,80,81,82,83,91];
if length(pVec) == 35 % then VOnly
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
ECVals = [1300;  % c_o (µM)
    147000.0;   % Na_o (µM)
    4000.0;      % K_o (µM)
    128000.0];  % Cl_o (µM)

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
pPSO(12) = SOCEFrac*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
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
ylabel('Applied current (nA)')
prettyGraph
subplot(4,1,2)
plot(Time_withSOCE, Y_withSOCE(:,5))
ylabel('SL voltage (mV)')
prettyGraph
subplot(4,1,3)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
ylabel('Myo calcium (uM)')
prettyGraph
subplot(4,1,4)
plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
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
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
xlim([0 .01])
ylabel('Myo calcium (uM)')
prettyGraph
subplot(3,1,2)
semilogy(Time_withSOCE, 1e-4*(fluxes(:,6)*JFrac+fluxes(:,8+6)*BFrac))
hold on
semilogy(Time_withSOCE, 1e-4*(fluxes(:,7)*JFrac+fluxes(:,8+7)*BFrac))
semilogy(Time_withSOCE, 1e-4*(fluxes(:,8)*JFrac+fluxes(:,8+8)*BFrac))
legend('RyR', 'SERCA', 'Leak')
ylabel('SR fluxes (uM/s) x1e4')
prettyGraph
xlim([0 .01])
subplot(3,1,3)
semilogy(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,8+1)*BFrac)
hold on
semilogy(Time_withSOCE, fluxes(:,2)*JFrac+fluxes(:,8+2)*BFrac)
semilogy(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,8+5)*BFrac)
legend('SOCE', 'SL leak', 'PMCA')
xlim([0 .01])
ylabel('PM fluxes (uM/s)')
xlabel('Time (s)')
prettyGraph

%% Dissect SOCE in baseline case
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
ylabel('Myo calcium (uM)')
prettyGraph
subplot(3,1,2)
plot(Time_withSOCE, Y_withSOCE(:,2)*JSRFrac + Y_withSOCE(:,27)*BSRFrac)
ylabel('SR calcium (uM)')
yyaxis right
plot(Time_withSOCE, Y_withSOCE(:,1))
ylabel('SOCE open probability')
prettyGraph
subplot(3,1,3)
plot(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,8+1)*BFrac)
ylabel('SOCE flux (uM/s)')
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
soceFactor = logspace(-1, 1, 6);
crefFactor = logspace(-1, 1, 6);
[soceFactor, crefFactor] = meshgrid(soceFactor, crefFactor);
ciSS = zeros(size(pmcaFactor));
cSRSS = zeros(size(pmcaFactor));
soceRef = pPSO(95);
crefRef = pVec(12)*p0(12);
tSol = [0, 30*60];

% pPSO(41) = 0.01*pVec(41)*p0(41); % tau SOCE
% pPSO(96) = 0.25*pVec(96)*p0(96); % total SR buffer

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
[TimeSS_noSOCE,ySS_noSOCE,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, phosphateAccum);
yinf_noSOCE = ySSFinal_noSOCE;
[Time_noSOCE,Y_noSOCE,~,fluxes2,currents2] = SkelMuscleCa_dydt(tSol, 0, yinf_noSOCE, pPSO, tic, 8, phosphateAccum); % compute time-dependent solution
figure
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
hold on
finalCa = zeros(size(soceFactor));
finalCaSR = zeros(size(soceFactor));

for i = 1:size(soceFactor,1)
    for j = 1:size(soceFactor,2)
        pPSO(12) = crefFactor(i,j)*crefRef*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        pPSO(95) = soceFactor(i,j)*soceRef; % gSOCE
        [TimeSS_withSOCE,ySS_withSOCE,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, phosphateAccum);
        yinf_withSOCE = ySSFinal_withSOCE;

        % test single frequency and plot calcium
        [Time_withSOCE,Y_withSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, 0, yinf_withSOCE, pPSO, tic, 7, phosphateAccum); % compute time-dependent solution
        plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
        drawnow
        finalCa(i,j) = Y_withSOCE(end,8)*JFrac + Y_withSOCE(end,32)*BFrac;
        finalCaSR(i,j) = Y_withSOCE(end,2)*JSRFrac + Y_withSOCE(end,27)*BSRFrac;
        fprintf('Done with soceFactor=%.2f, crefFactor=%.2f: final ca = %.2f and final SR Ca = %.2f\r\n',...
            soceFactor(i,j), crefFactor(i,j), finalCa(i,j), finalCaSR(i,j))
    end
end

%% plot SOCE vs no SOCE over time for baseline conditions (Fig 4)
pPSO(96) = 0.25*pVec(96)*p0(96); % total SR buffer
pPSO(12) = 0.5;%1.35*pVec(12)*p0(12); % set cref ratio
pPSO(95) = 1.35*pVec(95)*p0(95); % gSOCE
pPSO(41) = 1.0*pVec(41)*p0(41); % tau SOCE
tSol = [0,0.5];
[Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
    computeSol(pPSO, [], tSol, 1, 100, phosphateAccum);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[Time_noSOCE,Y_noSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100, phosphateAccum);

figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
hold on
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)
ylabel('Myo calcium (uM)')
prettyGraph
subplot(3,1,2)
plot(Time_withSOCE, Y_withSOCE(:,2))
hold on
plot(Time_noSOCE, Y_noSOCE(:,2))
ylabel('SR calcium (uM)')
prettyGraph
subplot(3,1,3)
plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
hold on
plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge)
ylabel('Rel force')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% test range of frequencies for SOCE vs no SOCE - Fig 5
tSol = [0, 0.5]; % 420 is max time for HIIT
freq = [1, 25, 50, 75, 100, 125, 150, 175, 200, 250];
forceRatio = zeros(size(freq));
peakForce = zeros(2, length(freq));
endPeakRatio = zeros(2, length(freq));
pPSO(96) = 0.25*pVec(96)*p0(96); % total SR buffer
pPSO(12) = 0.5; % set cref ratio
pPSO(95) = 5*pVec(95)*p0(95); % gSOCE
pPSO(41) = 1.0*pVec(41)*p0(41); % tau SOCE
phosphateAccum = true;
figure
for i = 1:length(freq)
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

function [t,y,fluxes,currents,maxCrossBridge,yinits] = computeSol(pPSO, yinit, tSol, expt, freq, phosphateAccum)
    juncLocLogic = true(1,31);
    juncLocLogic(17:21) = false; % cross bridges
    bulkLocLogic = true(1,31);
    bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
    if isempty(yinit)
        % load in initial condition starting estimate
        load yinit0.mat yinit0
        yinit0([24,26]) = pPSO(74:75);
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