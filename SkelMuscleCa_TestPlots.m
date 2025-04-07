%% Test SOCE vs no SOCE steady states after adjusting parameters
load p0Struct.mat p0Struct
p0 = p0Struct.data;
pSol = pVec;
if length(pVec) == 35 % then VOnly
    highSensIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    VOnly = true;
elseif length(pVec) == 95
    highSensIdx = [1:75,84:105];%1:105;%[12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
else % then fitting to both calcium and V using ca sens indices
    highSensIdx = 1:105;%[3,15,18,23,32,35,40,42,43,69,72,77,79,80,81,83,86,90,91];
    VOnlyIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    VOnlyStruct = load('pVec_VOnly.mat', 'pVec');
    pVecVOnly = VOnlyStruct.pVec;
    pVecVOnly = pVecVOnly(:); % be sure it is a column vector
    [VOnlyOnly,onlyIdx] = setdiff(VOnlyIdx, highSensIdx); % non overlapping indices
    p0(VOnlyOnly) = pVecVOnly(onlyIdx).*p0(VOnlyOnly); % set p0 according to previous estimation
end
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);

phosphateAccum = true;

pmcaFactor = 1.0;%logspace(-1, 1, 6);
sercaFactor = 1.0;%logspace(-1, 1, 6);
ncxFactor = 1.0;%logspace(-1, 1, 6);
[pmcaFactor, sercaFactor, ncxFactor] = meshgrid(pmcaFactor, sercaFactor, ncxFactor);
ciSS = zeros(size(pmcaFactor));
cSRSS = zeros(size(pmcaFactor));
sercaRef = pPSO(42);
pmcaRef = pPSO(43);
ncxRef = pPSO(20);
SOCEFrac = 0.5;%pPSO(12);

for i = 1:size(pmcaFactor,1)
    for j = 1:size(pmcaFactor,2)
        for k = 1:size(pmcaFactor,3)

            % pPSO(45) = 0.002*.5676;
            % load pBest.mat pVec
            % pVec(45) = 0;
            % pPSO = pVec;

            % fitIdx = [13,17,19,77,79,80,4,5,6,8,9,10];
            % pPSO(fitIdx) = pFit(:) .* pPSO(fitIdx);
            % pPSO(13) = 2*pPSO(13);

            % if desired, perturb parameters to test effects
            % pPSO(45) = 0*pPSO(45); % Na leak
            % pPSO(19) = 1*pPSO(19); % Na channel
            % number 69 is the rate of phosphate transport into SR
            % pPSO(69) = pPSO(69)*100
            % pPSO(21) = 0.1*pPSO(21); % NK pump

            % pPSO(13) = 10*pPSO(13); % stim current
            % pPSO(19) = pPSO(19); % sodium current
            % pPSO(77) = -50; % V_h (deactivation thresh)

            % try changing PMCA and SERCA systematically
            % pPSO(34) = 1*pPSO(34); % K_SERCA
            pPSO(42) = sercaFactor(i,j,k)*sercaRef; % nu_SERCA
            % pPSO(31) = 1*pPSO(31); % K_PMCA
            pPSO(43) = pmcaFactor(i,j,k)*pmcaRef; % g_PMCA
            pPSO(20) = ncxFactor(i,j,k)*ncxRef; % g_NCX
            % pPSO(44) = 2*pPSO(44); % SR calcium leak
            % pPSO(100) = 1.0*pPSO(100); % D_ion
            %
            pPSO(91) = 0.0;%*pPSO(91); % g0 DHPR
            % pPSO(89) = 0.1*pPSO(89); % SL leak calcium
            % pPSO(92) = 10*pPSO(92); % j0_RyR
            %
            % pPSO(52) = 0.5*pPSO(52); %Parv conc
            % pPSO(41) = 0.01*pPSO(41); % tau SOCE
            pPSO(96) = 0.25*pPSO(96); % total SR buffer
            % pPSO(97) = 2*pPSO(97); % K_CSQ
            % yinit(sum(juncLocLogic(1:31))) = pPSO(96)*0.5;
            % yinit(sum(juncLocLogic) + sum(bulkLocLogic(1:31))) = pPSO(96)*0.5;
            %
            % pPSO(65) = 1*pPSO(65); % h0

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

            % compute steady state solution without SOCE or phosphate accumulation
            pPSO0 = pPSO;
            pPSO0(95) = 0; % no SOCE
            [TimeSS_noSOCE,ySS_noSOCE,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, phosphateAccum);
            % pPSO(12) = 0.5;
            pPSO(12) = SOCEFrac*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
            pPSO(95) = pPSO(95)*5; % gSOCE
            [TimeSS_withSOCE,ySS_withSOCE,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, phosphateAccum);
            ciSS(i,j,k) = ySSFinal_withSOCE(8);
            cSRSS(i,j,k) = ySSFinal_withSOCE(2);
            yinf_noSOCE = ySSFinal_noSOCE;
            yinf_withSOCE = ySSFinal_withSOCE;
            fprintf('Done with pmca=%.2fx, serca=%.2fx, and ncx=%.2fx: ci = %.2f and cSR = %.2f\r\n', ...
                pmcaFactor(i,j,k), sercaFactor(i,j,k), ncxFactor(i,j,k), ciSS(i,j,k), cSRSS(i,j,k))
        end
    end
end

%% compute max cross bridge engagement - expt 10 is only crossbridge testing
% (all other dydt = 0)
yinit_maxCa = ySS_withSOCE(end,:);
yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
[Time_maxCa, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa, pPSO, tic, 10, phosphateAccum);
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
maxCrossBridge = Y_maxCa(end,crossbridgeIdx);

vol_Fiber = pi * (20 ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003 ;
vol_myo = 0.95 * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
diffusion_length = pPSO(99);%0.05;
vol_myoJ = SA_TT*diffusion_length;
JFrac = vol_myoJ / vol_myo;
BFrac = 1 - JFrac;

%% test single frequency and plot calcium with vs. without SOCE
tSol = [0, 10];%30*60]; % 420 is max time for HIIT
% yinf_withSOCE([2,27]) = 500;
% yinf_noSOCE([2,27]) = 500;
% pPSO(12) = pVec(12)*2000;
% yinf_withSOCE([5,29]) = -80;
% yinf_noSOCE([5,29]) = -80;
[Time_withSOCE,Y_withSOCE,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, 50, yinf_withSOCE, pPSO, tic, 1, phosphateAccum); % compute time-dependent solution
[Time_noSOCE,Y_noSOCE,~,fluxes2,currents2] = SkelMuscleCa_dydt(tSol, 50, yinf_noSOCE, pPSO, tic, 2, phosphateAccum); % compute time-dependent solution
figure
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
hold on
plot(Time_noSOCE, Y_noSOCE(:,8)*JFrac + Y_noSOCE(:,32)*BFrac)

%% Fig 3 plots
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
hold on
plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge)
ylabel('Rel force')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% Dissect currents for Fig 3
% current in pA/s = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                %    I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% plot sodium, KDR, and KIR currents
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,5))
xlim([0 .02])
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
xlim([0 .02])
% Potassium currents
subplot(3,1,3)
plot(Time_withSOCE, (currents(:,4)+currents(:,13+4))/1000)
hold on
plot(Time_withSOCE, (currents(:,5)+currents(:,13+5))/1000)
plot(Time_withSOCE, (currents(:,8)+currents(:,13+8))/1000)
legend('KDR', 'KIR', 'NaK')
xlim([0 .02])
ylabel('Potassium current (nA)')
xlabel('Time (s)')
prettyGraph

%% Dissect fluxes for Fig 3
% fluxes in µM/s = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
%                   LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];% 
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,8)*JFrac + Y_withSOCE(:,32)*BFrac)
% xlim([0 .02])
ylabel('Myo calcium (uM)')
prettyGraph
subplot(3,1,2)
plot(Time_withSOCE, 1e-4*(fluxes(:,6)*JFrac+fluxes(:,8+6)*BFrac))
hold on
plot(Time_withSOCE, 1e-4*(fluxes(:,7)*JFrac+fluxes(:,8+7)*BFrac))
plot(Time_withSOCE, 1e-4*(fluxes(:,8)*JFrac+fluxes(:,8+8)*BFrac))
legend('RyR', 'SERCA', 'Leak')
ylabel('SR fluxes (uM/s)')
prettyGraph
% xlim([0 .02])
subplot(3,1,3)
plot(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,8+1)*BFrac)
hold on
plot(Time_withSOCE, fluxes(:,2)*JFrac+fluxes(:,8+2)*BFrac)
plot(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,8+5)*BFrac)
legend('SOCE', 'SL leak', 'PMCA')
% xlim([0 .02])
ylabel('PM fluxes (uM/s)')
xlabel('Time (s)')
prettyGraph

%% SOCE vs no SOCE for figure 4
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