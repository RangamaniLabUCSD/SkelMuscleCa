%% SkelMuscleCa_MakeFigs script
% This file contains all the necessary code to generate figures for the paper 
% "Systems modeling reveals that store-operated calcium entry modulates
% force and fatigue during exercise"

%% Plot results of estimation
load pSol_fullWithMito/PSO_14-Jul-2026_withMito.mat pSol
% pSol([46,76]) = 1;
% pSol(95) = 2;
% load pSol_fullWithMito/pBest_withMito_Jul14.mat pVec
% load pSol_fullWithMito/pBest_withMito_VOnly.mat pVec
% load pSol_fullWithMito/PSO_05-Jul-2026_withMito_VOnly.mat pSol
% highSensIdx = [1,3,4,5,8,9,11,13,14,16,18,19,22,23,24,25,26,28,33,40,62,77,79,80,81,82,99];
% highSensIdxPrev = [1,3,5,8,9,11,13,14,16,18,22,23,24,25,26,28,33,40,46,76,77,79,80,81,82];
% pSolNew = ones(106,1);
% pSolNew(highSensIdxPrev) = pSol;
% pSolNew = pSolNew(highSensIdx);
% pSol = pVec;
% load PSO_24-Jun-2026_newRange.mat pSol
% load Data/pSol_allParam.mat pSol
% pSol(72) = 0; % no phosphate deg
% pSol(12) = 1;
pSol = [pSol, 1, 1, 1];
tic
SkelMuscleObj(pSol, true, '', true)
toc

%% Morris plots for Fig 2 and Fig S2
try uqlab
catch
    error('must first add UQLab to path before starting this section of code')
end
% load Data/MorrisResults.mat MorrisAnalysis
load MorrisResults_4e4_12-Jul-2026.mat MorrisAnalysis
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
Crossbridge_Cycle = [55:57,59:67,73,84,93,94,106,110:112];
Diffusion = 99:105; % parameters for diffusion between junctional and bulk compartments
Mito = 107:109;
% qois are [cSRSS, SLVoltSS, NaSS, ClSS, CaSS, KSS, FSS, MaxCa, MaxV, MaxPost, AvgCa, AvgPost, AvgVolt, VoltWidth] 
qoiList = [28+9,28+13, 8, 11]; % voltage qois from expt 8 (MaxV, AvgVolt), 
% calcium qois from expt 3 (MaxCa, AvgCa) 
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
    scatter(MorrisAnalysis.Results.MuStar(Mito,k),...
        MorrisAnalysis.Results.Std(Mito,k),'filled',...
        'MarkerFaceColor',[0,0,0])
    xlabel('mu*')
    ylabel('sigma')
    title(expt_names{k})
    hold on
    ten_percent = max(MorrisAnalysis.Results.MuStar(:,k))*0.2;
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
        'Diffusion','Mito','10% Threshold')
    prettyGraph
end

%% Test baseline model behavior - run baseline simulation for plots in Figs 3-4
load Data/p0Struct.mat p0Struct % load in default parameter values
% load PSO_24-Jun-2026_newRange.mat pSol
% load pBest_newWithMito.mat pVec
% pSol = pVec;
% load Data/pSol_allParam.mat pSol
% load pSol_fullWithMito/pBest_withMito_Jul14.mat pVec
% pSol = pVec;
load pSol_fullWithMito/PSO_14-Jul-2026_withMito.mat pSol
p0 = p0Struct.data;
% load Data/pSol_allParam.mat pSol
pSol(91) = 0; % no DHPR flux
% pSol = pVec(1:106);
pSol(12) = 1; % Keep cref at default value
% VOnlyIdx is the list of parameters to which SL voltage is sensitive
VOnlyIdx = [1,3,4,5,8,9,11,13,14,16,18,19,22,23,24,25,26,28,33,40,62,77,79,80,81,82,99];
% VOnlyIdx = [1,3,5,8,9,11,13,14,16,18,22,23,24,25,26,28,33,40,46,76,77,79,80,81,82];
% VOnlyIdx = [1,4,5,6,9,11,13,14,16,18,19,22,23,24,25,26,30,33,40,47,76,77,79,80,81,82,98];
% VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
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
pPSO(59) = 1*pSol(59)*p0(59); % kon1 (trop binding)
pPSO(75) = 1*pSol(75)*p0(75); % myo phosphate
pPSO(66) = 1*pSol(66)*p0(66); % phos-mediated reverse rate
pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion

% juncLocLogic and bulkLocLogic dictate whether a given state variable
% appears in junctional and/or bulk regions (see SkelMuscleCa_dydt)
juncLocLogic = true(1,40);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,40);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions

% define geometric quantities
diffusion_length = pPSO(99);
volFraction_mitoJ = 0.05; 
volFraction_mitoB = 0.05;
% exerciseParam = [6, 3];
exerciseParam = [0.65, 0.25*0.65];
withMito = true;
juncLocLogic_noMito = juncLocLogic(1:31);
bulkLocLogic_noMito = bulkLocLogic(1:31);
R_fiber = 20;
vol_Fiber = pi * (R_fiber ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003;
volFraction_SR = 0.05;
volFraction_myo = 1-volFraction_TT-volFraction_SR;
vol_myo = volFraction_myo * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
vol_myoJ = SA_TT*diffusion_length;
vol_myo_corrected = vol_myoJ*(1-volFraction_mitoJ) + (vol_myo-vol_myoJ)*(1-volFraction_mitoB);
vol_mito_tot = vol_myo - vol_myo_corrected;
JFrac = vol_myoJ*(1-volFraction_mitoJ) / vol_myo_corrected;
BFrac = 1 - JFrac;
vol_SA_SR = 1/10;
SRJ_occupancy = 0.5;
vol_SR = volFraction_SR*vol_Fiber;

geomParam = [vol_SA_ratio, volFraction_myo, volFraction_SR,...
             volFraction_TT, R_fiber, vol_SA_SR, SRJ_occupancy,...
             volFraction_mitoJ,volFraction_mitoB,3.002e-3];
geomParam(11:12) = [1,1];
% run simulation for 0.5 s at 100 Hz
tSol = [0, 1];%[0:.0001:0.2];%[0, 10.0];
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
SSPhosAccum = false; % whether phosphate accumulation is allowed in SS estimation (always false in tests for this paper)
% [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
%     computeSol(pPSO, [], tSol, [3,exerciseParam], 100, phosphateAccum, geomParam, SSPhosAccum, withMito);
mitoVolFractions = [0.05];%0.001, 0.003, 0.01, 0.03, 0.1, 0.3];
mitoVolSA = 3.002e-3*ones(size(mitoVolFractions));%,10.0];
CaOut = cell(length(mitoVolFractions),1);
CaMitoOut = cell(length(mitoVolFractions),1);
CaSROut = cell(length(mitoVolFractions),1);
phosOut = cell(length(mitoVolFractions),1);
ATPout = cell(length(mitoVolFractions),1);
forceOut = cell(length(mitoVolFractions),1);
for i = 1:length(mitoVolFractions)
    if mitoVolFractions(i) == 0
        [JFrac, BFrac, JSRFrac, BSRFrac, JFracMito, BFracMito] = getFracs(0,0,diffusion_length);
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
            computeSol(pPSO, [], tSol, 1, 100, phosphateAccum, geomParam, SSPhosAccum, false);
        [CaJuncIdx, CaBulkIdx] = getJuncBulkIdx(8, juncLocLogic_noMito, bulkLocLogic_noMito);
        CaOut{i} = JFrac*Y_withSOCE(:,CaJuncIdx) + BFrac*Y_withSOCE(:,CaBulkIdx);
        [CaMitoJuncIdx, CaMitoBulkIdx] = getJuncBulkIdx(33, juncLocLogic_noMito, bulkLocLogic_noMito);
        CaMitoOut{i} = JFracMito*Y_withSOCE(:,CaMitoJuncIdx) + BFracMito*Y_withSOCE(:,CaMitoBulkIdx);
        [CaSRJuncIdx, CaSRBulkIdx] = getJuncBulkIdx(2, juncLocLogic_noMito, bulkLocLogic_noMito);
        CaSROut{i} = JSRFrac*Y_withSOCE(:,CaSRJuncIdx) + BSRFrac*Y_withSOCE(:,CaSRBulkIdx);
        [PiJuncIdx, PiBulkIdx] = getJuncBulkIdx(26, juncLocLogic_noMito, bulkLocLogic_noMito);
        phosOut{i} = JFrac*Y_withSOCE(:,PiJuncIdx) + BFrac*Y_withSOCE(:,PiBulkIdx);
        [ATPJuncIdx, ATPBulkIdx] = getJuncBulkIdx(23, juncLocLogic_noMito, bulkLocLogic_noMito);
        ATPout{i} = JFrac*Y_withSOCE(:,ATPJuncIdx) + BFrac*Y_withSOCE(:,ATPBulkIdx);
        [~, crossbridgeIdx] = getJuncBulkIdx(21, juncLocLogic_noMito, bulkLocLogic_noMito);
        forceOut{i} = Y_withSOCE(:,crossbridgeIdx);
    else
        [JFrac, BFrac, JSRFrac, BSRFrac, JFracMito, BFracMito] =...
                    getFracs(mitoVolFractions(i), mitoVolFractions(i),diffusion_length);
        geomParamCur = geomParam;
        geomParamCur(8:9) = mitoVolFractions(i);
        geomParamCur(10) = mitoVolSA(i);
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
            computeSol(pPSO, [], tSol, 1, 100, phosphateAccum, geomParamCur, SSPhosAccum, true);
        [CaJuncIdx, CaBulkIdx] = getJuncBulkIdx(8, juncLocLogic, bulkLocLogic);
        CaOut{i} = JFrac*Y_withSOCE(:,CaJuncIdx) + BFrac*Y_withSOCE(:,CaBulkIdx);
        [CaMitoJuncIdx, CaMitoBulkIdx] = getJuncBulkIdx(33, juncLocLogic, bulkLocLogic);
        CaMitoOut{i} = JFracMito*Y_withSOCE(:,CaMitoJuncIdx) + BFracMito*Y_withSOCE(:,CaMitoBulkIdx);
        [CaSRJuncIdx, CaSRBulkIdx] = getJuncBulkIdx(2, juncLocLogic, bulkLocLogic);
        CaSROut{i} = JSRFrac*Y_withSOCE(:,CaSRJuncIdx) + BSRFrac*Y_withSOCE(:,CaSRBulkIdx);
        [PiJuncIdx, PiBulkIdx] = getJuncBulkIdx(26, juncLocLogic, bulkLocLogic);
        phosOut{i} = JFrac*Y_withSOCE(:,PiJuncIdx) + BFrac*Y_withSOCE(:,PiBulkIdx);
        [ATPJuncIdx, ATPBulkIdx] = getJuncBulkIdx(23, juncLocLogic, bulkLocLogic);
        ATPout{i} = JFrac*Y_withSOCE(:,ATPJuncIdx) + BFrac*Y_withSOCE(:,ATPBulkIdx);
        [~, crossbridgeIdx] = getJuncBulkIdx(21, juncLocLogic, bulkLocLogic);
        forceOut{i} = Y_withSOCE(:,crossbridgeIdx);
    end
end

%% "Effects of mito" plots
figure
legendEntries = cell(length(mitoVolFractions),1);
for i = 1:length(mitoVolSA)
    volCorrect = (1-mitoVolFractions(i));
    legendEntries{i} = sprintf('Mito vol frac = %.3f', mitoVolFractions(i));
    subplot(4,1,1)
    % yyaxis left
    % plot(tSol, CaOut{i}*volCorrect, 'LineWidth',1)%*(1-mitoVolFractions(i)))
    plot(tSol, CaSROut{i}, 'LineWidth',1)%*(1-mitoVolFractions(i)))
    hold on
    xlabel('Time (s)')
    ylabel('Myo calcium (µM)')
    % yyaxis right
    % plot(tSol, CaMitoOut{i}*mitoVolFractions(i), 'LineWidth', 1)
    % hold on
    prettyGraph
    subplot(4,1,2)
    plot(tSol, phosOut{i}*volCorrect, 'LineWidth',1)%*(1-mitoVolFractions(i)))
    xlabel('Time (s)')
    ylabel('Myo phosphate (µM)')
    prettyGraph
    hold on
    subplot(4,1,3)
    plot(tSol, ATPout{i}*volCorrect, 'LineWidth', 1)%*(1-mitoVolFractions(i)))
    xlabel('Time (s)')
    ylabel('Myo ATP (µM)')
    prettyGraph
    hold on
    subplot(4,1,4)
    plot(tSol, forceOut{i}*volCorrect/maxCrossBridge, 'LineWidth',1)%*(1-mitoVolFractions(i)))
    xlabel('Time (s)')
    ylabel('Rel force')
    prettyGraph
    hold on
end
legend(legendEntries)

%% Fig 3 plots
figure
subplot(4,1,1)
if withMito
    plot(Time_withSOCE, currents(:,14+13)/1000)
else
    plot(Time_withSOCE, currents(:,13+13)/1000)
end
[CaJuncIdx, CaBulkIdx] = getJuncBulkIdx(8, juncLocLogic, bulkLocLogic);
[CaSRJuncIdx, CaSRBulkIdx] = getJuncBulkIdx(2, juncLocLogic, bulkLocLogic);
[CaECJuncIdx, ~] = getJuncBulkIdx(27, juncLocLogic, bulkLocLogic);
[TTVoltIdx, SLVoltIdx] = getJuncBulkIdx(5, juncLocLogic, bulkLocLogic);
if withMito
    [CaMitoJuncIdx, CaMitoBulkIdx] = getJuncBulkIdx(33, juncLocLogic, bulkLocLogic);
end
hold on
xlim([0 0.2])
ylabel('Applied current (nA)')
ylim([0 110])
prettyGraph
subplot(4,1,2)
plot(Time_withSOCE, Y_withSOCE(:,TTVoltIdx))
hold on
ylabel('SL voltage (mV)')
xlim([0 0.2])
ylim([-85 25])
prettyGraph
subplot(4,1,3)
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
hold on
xlim([0 0.2])
ylim([0 30])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(4,1,4)
plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
hold on
xlim([0 0.2])
ylim([0 1])
ylabel('Rel force')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% Calcium and ATP plots
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac+Y_withSOCE(:,CaBulkIdx)*BFrac)
hold on
ylabel('Myo calcium (μM)')
prettyGraph
subplot(3,1,2)
ATPJuncIdx = sum(juncLocLogic(1:23));
ATPBulkIdx = sum(juncLocLogic)+sum(bulkLocLogic(1:23));
ATPvals = Y_withSOCE(:,ATPJuncIdx)*JFrac + Y_withSOCE(:,ATPBulkIdx)*BFrac; 
plot(Time_withSOCE, ATPvals)
hold on
ylabel('ATP')
prettyGraph
subplot(3,1,3)
phosJuncIdx = sum(juncLocLogic(1:26));
phosBulkIdx = sum(juncLocLogic)+sum(bulkLocLogic(1:26));
phosvals = Y_withSOCE(:,phosJuncIdx)*JFrac+Y_withSOCE(:,phosBulkIdx)*BFrac;
plot(Time_withSOCE, phosvals)
hold on
ylabel('Pi')
xlabel('Time (s)')
prettyGraph
set(gcf,'Renderer','painters')

%% ATP and phosphate mass conservation
tCur = Time_withSOCE;
yCur = Y_withSOCE;
[CaATPJuncIdx, CaATPBulkIdx] = getJuncBulkIdx(16, juncLocLogic, bulkLocLogic);
[MgATPJuncIdx, MgATPBulkIdx] = getJuncBulkIdx(22, juncLocLogic, bulkLocLogic);
[ATPMitoJuncIdx, ATPMitoBulkIdx] = getJuncBulkIdx(34, juncLocLogic, bulkLocLogic);
[~,DIdx] = getJuncBulkIdx(19, juncLocLogic, bulkLocLogic);
[~,PreIdx] = getJuncBulkIdx(20, juncLocLogic, bulkLocLogic);
[~,PostIdx] = getJuncBulkIdx(21, juncLocLogic, bulkLocLogic);
[ADPJuncIdx, ADPBulkIdx] = getJuncBulkIdx(32, juncLocLogic, bulkLocLogic);
[PiJuncIdx, PiBulkIdx] = getJuncBulkIdx(26, juncLocLogic, bulkLocLogic);
[PiMitoJuncIdx, PiMitoBulkIdx] = getJuncBulkIdx(35, juncLocLogic, bulkLocLogic);
[PiSRJuncIdx, PiSRBulkIdx] = getJuncBulkIdx(24, juncLocLogic, bulkLocLogic);
[PiCaJuncIdx, PiCaBulkIdx] = getJuncBulkIdx(25, juncLocLogic, bulkLocLogic);
[PiCaMitoJuncIdx, PiCaMitoBulkIdx] = getJuncBulkIdx(38, juncLocLogic, bulkLocLogic);
% Be careful to use correct volume factors here - JFrac, BFrac, JSRFrac,
% BSRFrac, JFracMito, BFracMito AND weight the relative volumes
JSRFracRel = JSRFrac * vol_SR / vol_myo_corrected;
BSRFracRel = BSRFrac * vol_SR / vol_myo_corrected;
JFracMitoRel = JFracMito * vol_mito_tot / vol_myo_corrected;
BFracMitoRel = BFracMito * vol_mito_tot / vol_myo_corrected;
ATPtot_myo = (yCur(:,CaATPJuncIdx) + yCur(:,MgATPJuncIdx) + yCur(:,ATPJuncIdx))*JFrac +...
             (yCur(:,CaATPBulkIdx) + yCur(:,MgATPBulkIdx) + yCur(:,ATPBulkIdx) + yCur(:,DIdx))*BFrac;
ATPtot_Mito = yCur(:,ATPMitoJuncIdx)*JFracMitoRel + yCur(:,ATPMitoBulkIdx)*BFracMitoRel;
ATPtot = ATPtot_myo + ATPtot_Mito;
ADPtot_myo = yCur(:,ADPJuncIdx)*JFrac + (yCur(:,PreIdx) + yCur(:,PostIdx) + yCur(:,ADPBulkIdx))*BFrac;
load  Data/pMitoNew.mat pMito
ANPtot_Mito = 15000*pMito(38);
ADPtot_Mito = (ANPtot_Mito - yCur(:,ATPMitoJuncIdx))*JFracMitoRel + (ANPtot_Mito - yCur(:,ATPMitoBulkIdx))*BFracMitoRel;
ADPtot = ADPtot_myo + ADPtot_Mito;
Pitot_myo = yCur(:,PiJuncIdx)*JFrac + (yCur(:,PreIdx) + yCur(:,PiBulkIdx))*BFrac;
Pitot_Mito = (yCur(:,PiMitoJuncIdx) + yCur(:,PiCaMitoJuncIdx))*JFracMitoRel +...
             (yCur(:,PiMitoBulkIdx) + yCur(:,PiCaMitoBulkIdx))*BFracMitoRel;
Pitot_SR = (yCur(:,PiSRJuncIdx) + yCur(:,PiCaJuncIdx))*JSRFracRel +...
           (yCur(:,PiSRBulkIdx) + yCur(:,PiCaBulkIdx))*BSRFracRel;
Pitot = Pitot_myo + Pitot_Mito + Pitot_SR;
figure
plot(tCur, ATPtot+ADPtot)
hold on
plot(tCur, ATPtot+Pitot)

%% Dissect currents for Fig 4A
% current in pA/s = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                %    I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% plot sodium, KDR, and KIR currents
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,SLVoltIdx))
xlim([0 .01])
ylabel('SL voltage (mV)')
prettyGraph
% sodium currents
subplot(3,1,2)
if withMito
    IbulkOffset = 14;
else
    IbulkOffset = 13;
end
plot(Time_withSOCE, (currents(:,10)+currents(:,IbulkOffset+10))/1000)
hold on
plot(Time_withSOCE, (currents(:,7)+currents(:,IbulkOffset+7))/1000)
plot(Time_withSOCE, (currents(:,9)+currents(:,IbulkOffset+9))/1000)
legend('Na', 'NCX', 'NaK')
ylabel('Sodium current (nA)')
prettyGraph
xlim([0 .01])
% Potassium currents
subplot(3,1,3)
plot(Time_withSOCE, (currents(:,4)+currents(:,IbulkOffset+4))/1000)
hold on
plot(Time_withSOCE, (currents(:,5)+currents(:,IbulkOffset+5))/1000)
plot(Time_withSOCE, (currents(:,8)+currents(:,IbulkOffset+8))/1000)
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
tMax = 0.01;
plotMaxes = false;
if withMito
    JbulkOffset = 13;
else
    JbulkOffset = 10;
end
figure
if plotMaxes
    subplot(3,1,1)
else
    subplot(5,1,1)
end
% if plotMaxes
%     [tMaxes,CaMaxes] = getMaxes(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac, 100);
%     plot(tMaxes, CaMaxes)
% else
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
% end
xlim([0 tMax])
ylabel('Myo calcium (μM)')
prettyGraph
if plotMaxes
    subplot(3,1,2)
    for fIdx = [6,8]%,12]
        [tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, fluxes(:,fIdx)*JFrac+fluxes(:,JbulkOffset+fIdx)*BFrac,100);
        semilogy(tMaxesCur, fluxMaxesCur)
        hold on
    end
else
    subplot(5,1,2)
    semilogy(Time_withSOCE, fluxes(:,6)*JFrac+fluxes(:,JbulkOffset+6)*BFrac)
    hold on
    semilogy(Time_withSOCE, fluxes(:,8)*JFrac+fluxes(:,JbulkOffset+8)*BFrac)
    semilogy(Time_withSOCE, fluxes(:,12)*JFrac+fluxes(:,JbulkOffset+12)*BFrac)
    legend('RyR', 'Leak','Mito NCX')
end
ylabel('SR fluxes (μM/s)')
prettyGraph
ylim([1e-5 5e5])
xlim([0 tMax])
mitoJuncFluxes = fluxes(:,11) - fluxes(:,12) + fluxes(:,13);
mitoBulkFluxes = fluxes(:,JbulkOffset+11) - fluxes(:,JbulkOffset+12) + fluxes(:,JbulkOffset+13);
if plotMaxes
    subplot(3,1,3)
    for fIdx = [7]%,11,13]
        [tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, fluxes(:,fIdx)*JFrac+fluxes(:,JbulkOffset+fIdx)*BFrac,100);
        semilogy(tMaxesCur, fluxMaxesCur)
        hold on
    end
    [tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, mitoJuncFluxes*JFrac+mitoBulkFluxes*BFrac,100);
    semilogy(tMaxesCur, fluxMaxesCur)
else
    subplot(5,1,3)
    semilogy(Time_withSOCE, fluxes(:,7)*JFrac+fluxes(:,JbulkOffset+7)*BFrac)
    hold on
    % semilogy(Time_withSOCE, mitoJuncFluxes*JFrac+mitoBulkFluxes*BFrac)
    semilogy(Time_withSOCE, fluxes(:,11)*JFrac+fluxes(:,JbulkOffset+11)*BFrac)
    semilogy(Time_withSOCE, fluxes(:,13)*JFrac+fluxes(:,JbulkOffset+13)*BFrac)
    legend('SERCA', 'MCU', 'mPTP')
end
xlim([0 tMax])
ylim([1e-5 5e5])
% ylim([-.5e4, 10e4])
set(gca, 'ydir', 'reverse')
% legend('RyR', 'Leak','SERCA')
% ylabel('SR fluxes (μM/s)')
prettyGraph
if plotMaxes
    subplot(3,1,2)
    for fIdx = [1,2]
        [tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, fluxes(:,fIdx)*JFrac+fluxes(:,JbulkOffset+fIdx)*BFrac,100);
        semilogy(tMaxesCur, fluxMaxesCur)
        hold on
    end
    legend('RyR', 'Leak', 'SOCE', 'SL leak')
else
    subplot(5,1,4)
    semilogy(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,JbulkOffset+1)*BFrac)
    hold on
    semilogy(Time_withSOCE, fluxes(:,2)*JFrac+fluxes(:,JbulkOffset+2)*BFrac)
    legend('SOCE', 'SL leak')
end
xlim([0 tMax])
% ylim([1e-0 1e5])
ylim([1e-3 1e2])
ylabel('PM fluxes (μM/s)')
prettyGraph
if plotMaxes
    subplot(3,1,3)
    for fIdx = [5,3]
        [tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, fluxes(:,fIdx)*JFrac+fluxes(:,JbulkOffset+fIdx)*BFrac,100);
        semilogy(tMaxesCur, fluxMaxesCur)
        hold on
    end
    legend('SERCA', 'Mito','PMCA','NCX')
else
    subplot(5,1,5)
    semilogy(Time_withSOCE, fluxes(:,5)*JFrac+fluxes(:,JbulkOffset+5)*BFrac)
    hold on
    semilogy(Time_withSOCE, fluxes(:,3)*JFrac+fluxes(:,JbulkOffset+3)*BFrac)
    legend('PMCA','NCX')
end
xlim([0 tMax])
% ylim([1e-0 1e5])
ylim([1e-3 1e2])
% ylabel('PM fluxes (μM/s)')
set(gca, 'ydir', 'reverse')
xlabel('Time (s)')
prettyGraph

%% Plot SOCE flux separately
figure
subplot(3,1,1)
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
prettyGraph
ylabel('Myo calcium')
hold on
subplot(3,1,2)
plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx)*JFrac + Y_withSOCE(:,CaSRBulkIdx)*BFrac)
prettyGraph
ylabel('SR calcium')
yyaxis right
[tMaxesCur,fluxMaxesCur] = getMaxes(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,JbulkOffset+1)*BFrac,100);
plot(tMaxesCur, fluxMaxesCur)
hold on
subplot(3,1,3)
plot(tMaxesCur, fluxMaxesCur)
hold on
plot(Time_withSOCE, fluxes(:,1)*JFrac+fluxes(:,JbulkOffset+1)*BFrac)
% xlim([0 1])
prettyGraph
ylabel('SOCE flux')
xlabel('Time (s)')

%% Plot calcium in different compartments (Fig S3)
tMax = 0.5;
figure
subplot(1,4,1)
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx))
hold on
plot(Time_withSOCE, Y_withSOCE(:,CaBulkIdx))
xlim([0 tMax])
legend('Junc','Bulk')
ylabel('Myo calcium (μM)')
prettyGraph
subplot(1,4,2)
plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx))
hold on
plot(Time_withSOCE, Y_withSOCE(:,CaSRBulkIdx))
ylabel('SR calcium (μM)')
xlabel('Time (s)')
xlim([0 tMax])
prettyGraph
subplot(1,4,3)
plot(Time_withSOCE, Y_withSOCE(:,CaECJuncIdx))
hold on
plot(Time_withSOCE, 1300*ones(size(Time_withSOCE)))
ylabel('EC calcium (μM)')
xlabel('Time (s)')
xlim([0 tMax])
prettyGraph
subplot(1,4,4)
plot(Time_withSOCE, Y_withSOCE(:,CaMitoJuncIdx))
hold on
plot(Time_withSOCE, Y_withSOCE(:,CaMitoBulkIdx))
ylabel('Mito calcium (μM)')
xlabel('Time (s)')
xlim([0 tMax])
prettyGraph

%% Plot mito potential changes
figure
MitoVJuncIdx = sum(juncLocLogic(1:37));
MitoVBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:37));
subplot(2,1,1)
plot(Time_withSOCE, Y_withSOCE(:,MitoVJuncIdx))
hold on
plot(Time_withSOCE, Y_withSOCE(:,MitoVBulkIdx))
xlim([0 1])
prettyGraph
subplot(2,1,2)
plot(Time_withSOCE, currents(:,14)+currents(:,IbulkOffset+14))
xlim([0 1])
% ylim([5e-2,5e6])
legend('Mito current')
prettyGraph

%% Plot mito calcium dynamics
figure
MitoCaJuncIdx = sum(juncLocLogic(1:33));
MitoCaBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:33));
subplot(2,1,1)
plot(Time_withSOCE, JFracMito*Y_withSOCE(:,MitoCaJuncIdx)+BFracMito*Y_withSOCE(:,MitoCaBulkIdx))
xlim([0 0.2])
prettyGraph
subplot(2,1,2)
semilogy(Time_withSOCE, fluxes(:,11)*JFracMito+fluxes(:,JbulkOffset+11)*BFracMito)
% plot(Time_withSOCE, fluxes(:,11)*JFracMito+fluxes(:,JbulkOffset+11)*BFracMito)
hold on
plot(Time_withSOCE, fluxes(:,12)*JFracMito+fluxes(:,JbulkOffset+12)*BFracMito)
plot(Time_withSOCE, fluxes(:,13)*JFracMito+fluxes(:,JbulkOffset+13)*BFracMito)
xlim([0 0.2])
% ylim([-200, 200])
ylim([1e-3 1e4])
ylabel('Mito fluxes (μM/s)')
legend('MCU', 'NCX','mPTP')
xlabel('Time (s)')
prettyGraph

%% Plot crossbridge cycling for Fig S1
tMax = 0.2;
figure
subplot(4,1,1)
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
hold on
ylabel('Myo calcium (μM)')
xlim([0 tMax])
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
xlim([0 tMax])
prettyGraph
subplot(4,1,3)
plot(Time_withSOCE, CaCaTrop)
hold on
ylabel('CaCaTrop (μM)')
xlim([0 tMax])
prettyGraph
subplot(4,1,4)
plot(Time_withSOCE, Apost)
hold on
ylabel('Apost (μM)')
prettyGraph
xlim([0 tMax])
xlabel('Time (s)')

%% SOCE vs no SOCE for Fig 5 TG test (Fig 5A)
% code allows for testing a range of SOCE conductance (soceFactor) and
% STIM1 calcium sensitivity (crefFactor), currently just keeps crefFactor
% at 1 for all tests
soceFactor = [0.3,1,3];
crefFactor = 1; % could be vector to test range
[soceFactor, crefFactor] = meshgrid(soceFactor, crefFactor);
ciSS = zeros(size(soceFactor));
cSRSS = zeros(size(soceFactor));
soceRef = p0(95);%pSol(95)*p0(95);
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
yinit0 = [yinit0; 100; 0.1; 5000; 5000; 50; 150; 0.0; 0.0; 0.0];
juncLocLogic = true(1,40);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,40);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
yinit = [yinit0(juncLocLogic); yinit0(bulkLocLogic)];

% define geometric quantities
vol_Fiber = pi * (20 ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003;
volFraction_SR = 0.05;
volFraction_myo = 1 - volFraction_SR - volFraction_TT; 
vol_myo = volFraction_myo * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
diffusion_length = pPSO(99);
vol_myoJ = SA_TT*diffusion_length;
vol_mitoJ = vol_myoJ * 0.05;
vol_mitoB = (vol_myo-vol_myoJ) * 0.05;
JFrac = (vol_myoJ-vol_mitoJ) / (vol_myo-vol_mitoJ-vol_mitoB);
BFrac = 1 - JFrac;
vol_SR = volFraction_SR*vol_Fiber;
SRJ_occupancy = 0.5;
SA_SRJ = SA_TT * SRJ_occupancy;
vol_SRJ = SA_SRJ*diffusion_length;
JSRFrac = vol_SRJ / vol_SR;
BSRFrac = 1 - JSRFrac;

% compute steady state solution without SOCE
phosphateAccum = false;
pPSO0 = pPSO;
pPSO0(95) = 0; % no SOCE
[TimeSS_noSOCE,ySS_noSOCE,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt_withMito(0:5:5000,0, yinit, pPSO0, tic, 2, false);%phosphateAccum);
yinf_noSOCE = ySSFinal_noSOCE;
[Time_noSOCE,Y_noSOCE,~,fluxes2,currents2] = SkelMuscleCa_dydt_withMito(tSol, 0, yinf_noSOCE, pPSO, tic, 8, phosphateAccum); % compute time-dependent solution
figure
subplot(2,1,1)
plot(Time_noSOCE, Y_noSOCE(:,CaJuncIdx)*JFrac + Y_noSOCE(:,CaBulkIdx)*BFrac)
hold on
prettyGraph
ylabel('Myo calcium (μM)')
subplot(2,1,2)
plot(Time_noSOCE, Y_noSOCE(:,CaSRJuncIdx)*JSRFrac + Y_noSOCE(:,CaSRBulkIdx)*BSRFrac)
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
        pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion
        [TimeSS_withSOCE,ySS_withSOCE,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt_withMito([0 1000], 0, yinit, pPSO, tic, 1, false);%phosphateAccum);
        yinf_withSOCE = ySSFinal_withSOCE;

        % test single frequency and plot calcium
        [Time_withSOCE,Y_withSOCE] = SkelMuscleCa_dydt_withMito(tSol, 0, yinf_withSOCE, pPSO, tic, 7, phosphateAccum); % compute time-dependent solution
        subplot(2,1,1)
        plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
        subplot(2,1,2)
        plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx)*JSRFrac + Y_withSOCE(:,CaSRBulkIdx)*BSRFrac)
        drawnow
        finalCa(i,j) = Y_withSOCE(end,8)*JFrac + Y_withSOCE(end,33+6)*BFrac;
        finalCaSR(i,j) = Y_withSOCE(end,2)*JSRFrac + Y_withSOCE(end,33+1)*BSRFrac;
        fprintf('Done with soceFactor=%.2f, crefFactor=%.2f: final ca = %.2f and final SR Ca = %.2f\r\n',...
            soceFactor(i,j), crefFactor(i,j), finalCa(i,j), finalCaSR(i,j))
    end
end

%% plot SOCE vs no SOCE over time for baseline conditions (Fig 5B-C)
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 0.002;%20*p0(95);%pSol(95)*p0(95); % SOCE conductance
pPSO(100) = 1.0*pSol(100)*p0(100); % ion diffusion
pPSO(59) = 1*pSol(59)*p0(59); % kon1 (trop binding)
geomParamCur = geomParam;
geomParamCur(8:9) = 0.05;
geomParamCur(10) = 3.002e-3;
geomParamCur(11:12) = [1,1];
phosphateAccum = true;

% test 1 s of stimulus here
tSol = [0, 1];
SSPhosAccum = false;
[Time_withSOCE,Y_withSOCE,~,~,maxCrossBridge, yinits] =...
    computeSol(pPSO, [], tSol, 1, 100, phosphateAccum, geomParamCur, SSPhosAccum, withMito);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[Time_noSOCE,Y_noSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100,...
    phosphateAccum, geomParamCur, SSPhosAccum, withMito);

[CaJuncIdx, CaBulkIdx] = getJuncBulkIdx(8, juncLocLogic, bulkLocLogic);
[CaSRJuncIdx, CaSRBulkIdx] = getJuncBulkIdx(2, juncLocLogic, bulkLocLogic);
figure
subplot(3,1,1)
plot(Time_noSOCE, Y_noSOCE(:,CaJuncIdx)*JFrac + Y_noSOCE(:,CaBulkIdx)*BFrac)
hold on
plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
% xlim([0 2])
ylabel('Myo calcium (μM)')
prettyGraph
subplot(3,1,2)
plot(Time_noSOCE, Y_noSOCE(:,CaSRJuncIdx)*JSRFrac + Y_noSOCE(:,CaSRBulkIdx)*BSRFrac)
hold on
plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx)*JSRFrac + Y_withSOCE(:,CaSRBulkIdx)*BSRFrac)
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

%% Try population behavior for soce vs no soce
Crossbridge_Cycle = [55:57,59:67,73,84,93,94,106];
pRef = pSol(:).*p0(:);
numTests = 100;
popParam = exp(log(2)*randn(numTests,length(Crossbridge_Cycle)));
tSol = 0:0.001:1;
forceDiff = zeros(numTests,length(tSol));
figure
hold on
for i = 1:numTests
    fprintf('Running case %d\r\n', i)
    pCur = pRef;
    pCur(Crossbridge_Cycle) = pCur(Crossbridge_Cycle) .* popParam(i,:)';
    pCur(12) = 1*p0(12); % set cref ratio
    pCur(95) = 20*pSol(95)*p0(95); % SOCE conductance
    pCur(100) = 1*pSol(100)*p0(100); % ion diffusion
    phosphateAccum = true; 
    % test 1 s of stimulus here
    SSPhosAccum = false;
    [Time_withSOCE,Y_withSOCE,~,~,maxCrossBridge, yinits] =...
        computeSol(pCur, [], tSol, 1, 100, phosphateAccum, geomParam, SSPhosAccum, withMito);
    pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
    [Time_noSOCE,Y_noSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100,...
        phosphateAccum, geomParam, SSPhosAccum, withMito);
    lastIdx = min([length(Y_withSOCE(:,1)),length(Y_noSOCE(:,1))]);
    forceDiff(i,1:lastIdx) = Y_withSOCE(1:lastIdx,crossbridgeIdx) - Y_noSOCE(1:lastIdx,crossbridgeIdx); 
    plot(tSol, forceDiff(i,:))
    drawnow
    if min(forceDiff(i,:)) < -.001
        fprintf('pause here\r\n')
    end
end

%% plot phosphate accumulation and force w/o phophate accum (Fig S4)
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 1*pSol(95)*p0(95); % gSOCE
pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion
tSol = [0, 2];
[tWithPhos,yWithPhos,fluxes,currents,maxCrossBridge,yinits] = computeSol(pPSO, [], tSol, 1, 100, true, geomParam, false, withMito);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[tNoPhos,yNoPhos] = computeSol(pPSO, yinits{2}, tSol, 1, 100, false, geomParam, false, withMito);
geomParamAlt = geomParam;
geomParamAlt(11) = 10;
% geomParamAlt(10) = geomParam(10);
% pPSO(68) = 1e6; % PP (threshold for phosphate precip formation)
[tMaxPiC,yMaxPiC] = computeSol(pPSO, yinits{2}, tSol, 1, 100, true, geomParamAlt, false, withMito);
% pPSO(68) = pSol(68)*p0(68); % reset PP (threshold for phosphate precip formation)
phosJuncIdx = sum(juncLocLogic(1:26));
phosBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:26));
phosSRJuncIdx = sum(juncLocLogic(1:24));
phosSRBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:24));
CaPhosJuncIdx = sum(juncLocLogic(1:25));
CaPhosBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:25));
figure
subplot(4,1,1)
plot(tWithPhos, yWithPhos(:,crossbridgeIdx)/maxCrossBridge)
hold on
plot(tNoPhos, yNoPhos(:,crossbridgeIdx)/maxCrossBridge)
plot(tMaxPiC, yMaxPiC(:,crossbridgeIdx)/maxCrossBridge)
prettyGraph
subplot(4,1,2)
plot(tWithPhos, yWithPhos(:,phosJuncIdx)*JFrac + yWithPhos(:,phosBulkIdx)*BFrac)
hold on
plot(tNoPhos, yNoPhos(:,phosJuncIdx)*JFrac + yNoPhos(:,phosBulkIdx)*BFrac)
plot(tMaxPiC, yMaxPiC(:,phosJuncIdx)*JFrac + yMaxPiC(:,phosBulkIdx)*BFrac)
ylabel("Myo phos")
prettyGraph
subplot(4,1,3)
plot(tWithPhos, yWithPhos(:,PiMitoJuncIdx)*JFracMito + yWithPhos(:,PiMitoBulkIdx)*BFracMito)
hold on
plot(tNoPhos, yNoPhos(:,PiMitoJuncIdx)*JFracMito + yNoPhos(:,PiMitoBulkIdx)*BFracMito)
plot(tMaxPiC, yMaxPiC(:,PiMitoJuncIdx)*JFracMito + yMaxPiC(:,PiMitoBulkIdx)*BFracMito)
ylabel("Mito phos")
prettyGraph
subplot(4,1,4)
plot(tWithPhos, yWithPhos(:,ATPJuncIdx)*JFrac + yWithPhos(:,ATPBulkIdx)*BFrac)
hold on
plot(tNoPhos, yNoPhos(:,ATPJuncIdx)*JFrac + yNoPhos(:,ATPBulkIdx)*BFrac)
plot(tMaxPiC, yMaxPiC(:,ATPJuncIdx)*JFrac + yMaxPiC(:,ATPBulkIdx)*BFrac)
ylabel("Myo ATP")
prettyGraph

%% plot hydrolysis rates (Fig S4)
if withMito
    JbulkOffset = 13;
    IbulkOffset = 14;
else
    JbulkOffset = 10;
    IbulkOffset = 13;
end
F = 96485.3321;
Jhyd = fluxes(:,9)*JFrac + fluxes(:,JbulkOffset+9)*BFrac;
J_NKX_ATP = 1e9*((currents(:,8) + currents(:,IbulkOffset+8)))/(2*F*vol_myo);  % one ATP per two potassium pumped, uM/s
J_SERCA_ATP = (fluxes(:,7)*JFrac+fluxes(:,JbulkOffset+7)*BFrac)/2; % uM/s, 1 ATP transports 2 calcium ions
J_PMCA_ATP = fluxes(:,5)*JFrac+fluxes(:,JbulkOffset+5)*BFrac; % uM/s
J_XB_ATP = pPSO(63)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:19)))*BFrac - ...
    pPSO(64)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:20)))*BFrac;
J_XB_ATP2 = pPSO(65)*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:20)))*BFrac -...
    (pPSO(66)/3000)*yWithPhos(:,PiBulkIdx).*yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:21)))*BFrac;
ATPjunc = yWithPhos(:,sum(juncLocLogic(1:23)));
ATPbulk = yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:23)));
J_baselinehyd = (pPSO(55)*(ATPjunc./(pPSO(56)+ATPjunc))*JFrac +...
    pPSO(55)*(ATPbulk./(pPSO(56)+ATPbulk))*BFrac);

figure
semilogy(tWithPhos, Jhyd)
hold on
plot(tWithPhos, J_NKX_ATP)
plot(tWithPhos, J_SERCA_ATP)
plot(tWithPhos, J_PMCA_ATP)
plot(tWithPhos, J_XB_ATP)
% plot(tWithPhos, J_XB_ATP2)
plot(tWithPhos, J_baselinehyd)
xlim([0 2])
ylim([1e0 5e3])
prettyGraph
xlabel('Time (s)')
ylabel('Hydrolysis rate (μM/s)')
legend('tot','NKX','SERCA','PMCA','XB','baseline')

%% plot ATP gen rates
ATPjunc = yWithPhos(:,sum(juncLocLogic(1:23)));
ATPbulk = yWithPhos(:,sum(juncLocLogic(:)) + sum(bulkLocLogic(1:23)));
[PiCaJuncIdx, PiCaBulkIdx] = getJuncBulkIdx(25, juncLocLogic, bulkLocLogic);
[PiCaMitoJuncIdx, PiCaMitoBulkIdx] = getJuncBulkIdx(38, juncLocLogic, bulkLocLogic);
[ADPMitoJuncIdx, ADPMitoBulkIdx] = getJuncBulkIdx(39, juncLocLogic, bulkLocLogic);
[VMitoJuncIdx, VMitoBulkIdx] = getJuncBulkIdx(37, juncLocLogic, bulkLocLogic);
% Be careful to use correct volume factors here - JFrac, BFrac, JSRFrac,
% BSRFrac, JFracMito, BFracMito AND weight the relative volumes
ATPtot_myo_junc = yWithPhos(:,CaATPJuncIdx) + yWithPhos(:,MgATPJuncIdx) + yWithPhos(:,ATPJuncIdx);
ATPtot_myo_bulk = yWithPhos(:,CaATPBulkIdx) + yWithPhos(:,MgATPBulkIdx) + yWithPhos(:,ATPBulkIdx);
ATP_Mito_junc = yWithPhos(:,ATPMitoJuncIdx);
ATP_Mito_bulk = yWithPhos(:,ATPMitoBulkIdx);
ADP_myo_junc = yWithPhos(:,ADPJuncIdx);
ADP_myo_bulk = yWithPhos(:,ADPBulkIdx);
ADP_Mito_junc = yWithPhos(:,ADPMitoJuncIdx);
ADP_Mito_bulk = yWithPhos(:,ADPMitoBulkIdx);
mitoStruct = load('PSOMito_11-Jul-2026.mat', 'pSol');
pMito = mitoStruct.pSol;
V_F1FO = 35000*pMito(22); % uM/s
p6 = 190*pMito(23); %mV
p7 = 8.5*pMito(24); %mV
K_iATP = 10000*pMito(25); %uM
V_ANT = 5000*pMito(26); %uM/s
alpha_c = 0.111*pMito(27);
alpha_m = 0.139*pMito(28);
f_ANT = 0.5*pMito(29);

mitoMyoRef = 0.05/(1-0.05);
mitoMyoCur = mitoMyoRef; % change if changing vol ratio
ATPset = pSol(106)*p0(106);
k_GLY_junc = 450*pMito(12)* ((ATPset-ATPjunc)/10) .* (yWithPhos(:,phosJuncIdx)/1000) * (mitoMyoRef/mitoMyoCur);
k_GLY_bulk = 450*pMito(12)* ((ATPset-ATPbulk)/10) .* (yWithPhos(:,phosBulkIdx)/1000) * (mitoMyoRef/mitoMyoCur);
NADH_ATP_conv = 1;
J_Gly = NADH_ATP_conv*(k_GLY_junc*JFracMito + k_GLY_bulk*BFracMito)*mitoMyoCur;
F = 96485.3321;
R = 8314.46261815;
T = 295.15;

V_Mito_junc = yWithPhos(:,VMitoJuncIdx);
V_Mito_bulk = yWithPhos(:,VMitoBulkIdx);
R_c_junc = alpha_c*ATPtot_myo_junc./ADP_myo_junc;
R_c_bulk = alpha_c*ATPtot_myo_bulk./ADP_myo_bulk;
R_m_junc = ADP_Mito_junc./(alpha_m*ATP_Mito_junc);
R_m_bulk = ADP_Mito_bulk./(alpha_m*ATP_Mito_bulk);
ANT_num_junc = (1-R_c_junc.*R_m_junc.*exp(-F*V_Mito_junc/(R*T)));
ANT_num_bulk = (1-R_c_bulk.*R_m_bulk.*exp(-F*V_Mito_bulk/(R*T)));
ANT_denom_junc = (1+R_c_junc.*exp(-f_ANT*F*V_Mito_junc/(R*T))).*(1+R_m_junc);
ANT_denom_bulk = (1+R_c_bulk.*exp(-f_ANT*F*V_Mito_bulk/(R*T))).*(1+R_m_bulk);
mitoSAFactor = 1.0;
J_ANT_junc = mitoSAFactor*V_ANT*ANT_num_junc./ANT_denom_junc;
J_ANT_bulk = mitoSAFactor*V_ANT*ANT_num_bulk./ANT_denom_bulk;
J_ANT = (J_ANT_junc*JFracMito + J_ANT_bulk*BFracMito)*mitoMyoCur;

J_F1FO_junc = mitoSAFactor*V_F1FO*(1./(1+exp((p6-V_Mito_junc)/p7))).*(K_iATP./(K_iATP+ATP_Mito_junc)).*(yWithPhos(:,PiMitoJuncIdx)/1000);
J_F1FO_bulk = mitoSAFactor*V_F1FO*(1./(1+exp((p6-V_Mito_bulk)/p7))).*(K_iATP./(K_iATP+ATP_Mito_bulk)).*(yWithPhos(:,PiMitoBulkIdx)/1000);
J_F1FO = (J_F1FO_junc*JFracMito + J_F1FO_bulk*BFracMito)*mitoMyoCur;

figure
plot(tWithPhos, J_Gly)
hold on
% plot(tWithPhos, J_ANT)
plot(tWithPhos, J_F1FO)
xlim([0 2])
prettyGraph
xlabel('Time (s)')
ylabel('ATP gen rate (μM/s)')
legend('Glycolysis','OxPhos')

%% plot calcium, SR calcium, phosphate, and force after init for soce vs no soce (Fig S5)
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 20*pSol(95)*p0(95); % gSOCE
pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion
tSol = [0, 10.5];
[tWithSOCE,yWithSOCE,fluxes,currents,maxCrossBridge,yinits] = computeSol(pPSO, [], tSol, 5, 100, true, geomParam, false, withMito);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[tNoSOCE,yNoSOCE] = computeSol(pPSO, yinits{1}, tSol, 6, 100, true, geomParam, false, withMito);
geomParamAlt = geomParam;
geomParamAlt(11) = 0.1;
% geomParamAlt(10) = geomParam(10);
[tLowPiC,yLowPiC] = computeSol(pPSO, yinits{1}, tSol, 6, 100, true, geomParamAlt, false, withMito);
phosJuncIdx = sum(juncLocLogic(1:26));
phosBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:26));
figure
subplot(4,1,1)
plot(tWithSOCE, yWithSOCE(:,CaJuncIdx)*JFrac + yWithSOCE(:,CaBulkIdx)*BFrac)
hold on
plot(tNoSOCE, yNoSOCE(:,CaJuncIdx)*JFrac + yNoSOCE(:,CaBulkIdx)*BFrac)
plot(tLowPiC, yLowPiC(:,CaJuncIdx)*JFrac + yLowPiC(:,CaBulkIdx)*BFrac)
ylabel("Myo calcium")
prettyGraph
subplot(4,1,2)
plot(tWithSOCE, yWithSOCE(:,CaSRJuncIdx)*JFrac + yWithSOCE(:,CaSRBulkIdx)*BFrac)
hold on
plot(tNoSOCE, yNoSOCE(:,CaSRJuncIdx)*JFrac + yNoSOCE(:,CaSRBulkIdx)*BFrac)
plot(tLowPiC, yLowPiC(:,CaSRJuncIdx)*JFrac + yLowPiC(:,CaSRBulkIdx)*BFrac)
ylabel("SR calcium")
prettyGraph
subplot(4,1,3)
plot(tWithSOCE, yWithSOCE(:,crossbridgeIdx)/maxCrossBridge)
hold on
plot(tNoSOCE, yNoSOCE(:,crossbridgeIdx)/maxCrossBridge)
plot(tLowPiC, yLowPiC(:,crossbridgeIdx)/maxCrossBridge)
prettyGraph
subplot(4,1,4)
plot(tWithSOCE, yWithSOCE(:,phosJuncIdx)*JFrac + yWithSOCE(:,phosBulkIdx)*BFrac)
hold on
plot(tNoSOCE, yNoSOCE(:,phosJuncIdx)*JFrac + yNoSOCE(:,phosBulkIdx)*BFrac)
plot(tLowPiC, yLowPiC(:,phosJuncIdx)*JFrac + yLowPiC(:,phosBulkIdx)*BFrac)
ylabel("Myo phos")
prettyGraph

%% plot calcium with vs without SOCE for different initial conditions
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 20*pSol(95)*p0(95); % gSOCE
pPSO(100) = 0.1*pSol(100)*p0(100); % ion diffusion
tSol = [0, 1];
[tWithSOCE,yWithSOCE,~,~,maxCrossBridge,yinits] = computeSol(pPSO, [], tSol, 1, 100, true, geomParam, false, withMito);
pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
[tNoSOCE,yNoSOCE] = computeSol(pPSO, yinits{1}, tSol, 2, 100, true, geomParam, false, withMito);
[tNoSOCE_altInit,yNoSOCE_altInit] = computeSol(pPSO, yinits{2}, tSol, 2, 100, true, geomParam, false, withMito);
figure
subplot(2,1,1)
[t1,y1] = getMaxes(tNoSOCE, yNoSOCE(:,CaJuncIdx)*JFrac + yNoSOCE(:,CaBulkIdx)*BFrac, 100);
% plot(t1,y1)
plot(tNoSOCE, yNoSOCE(:,CaJuncIdx)*JFrac + yNoSOCE(:,CaBulkIdx)*BFrac)
hold on
[t2,y2] = getMaxes(tWithSOCE, yWithSOCE(:,CaJuncIdx)*JFrac + yWithSOCE(:,CaBulkIdx)*BFrac, 100);
% plot(t2,y2)
plot(tWithSOCE, yWithSOCE(:,CaJuncIdx)*JFrac + yWithSOCE(:,CaBulkIdx)*BFrac)
[t3,y3] = getMaxes(tNoSOCE_altInit, yNoSOCE_altInit(:,CaJuncIdx)*JFrac + yNoSOCE_altInit(:,CaBulkIdx)*BFrac, 100);
% plot(t3, y3)
plot(tNoSOCE_altInit, yNoSOCE_altInit(:,CaJuncIdx)*JFrac + yNoSOCE_altInit(:,CaBulkIdx)*BFrac)
ylabel('Myo calcium (μM)')
prettyGraph
subplot(2,1,2)
[t4,y4] = getMaxes(tNoSOCE, yNoSOCE(:,CaSRJuncIdx)*JSRFrac + yNoSOCE(:,CaSRBulkIdx)*BSRFrac, 100);
% plot(t4, y4)
plot(tNoSOCE, yNoSOCE(:,CaSRJuncIdx)*JSRFrac + yNoSOCE(:,CaSRBulkIdx)*BSRFrac)
hold on
[t5,y5] = getMaxes(tWithSOCE, yWithSOCE(:,CaSRJuncIdx)*JSRFrac + yWithSOCE(:,CaSRBulkIdx)*BSRFrac, 100);
% plot(t5,y5)
plot(tWithSOCE, yWithSOCE(:,CaSRJuncIdx)*JSRFrac + yWithSOCE(:,CaSRBulkIdx)*BSRFrac)
[t6,y6] = getMaxes(tNoSOCE_altInit, yNoSOCE_altInit(:,CaSRJuncIdx)*JSRFrac + yNoSOCE_altInit(:,CaSRBulkIdx)*BSRFrac, 100);
% plot(t6,y6)
plot(tNoSOCE_altInit, yNoSOCE_altInit(:,CaSRJuncIdx)*JSRFrac + yNoSOCE_altInit(:,CaSRBulkIdx)*BSRFrac)
ylabel('SR calcium (μM)')
prettyGraph

%% test range of frequencies for SOCE vs no SOCE - Fig 6
tSol = [0, 10.5];
freq = 5:5:200;%5:5:200;
inclAltPhos = true; % whether to include additional test with 150% resting phosphate in SOCE KOs
peakForce = zeros(2+inclAltPhos, length(freq));
endPeakRatio = zeros(2+inclAltPhos, length(freq));
pPSO(12) = 1*p0(12); % set cref ratio
pPSO(95) = 0.002;%20*pSol(95)*p0(95); % gSOCE
pPSO(100) = 1.0*pSol(100)*p0(100); % ion diffusion
pPSO(75) = 1*pSol(75)*p0(75); % myo phosphate
geomParam(8:9) = 0.05;
geomParam(10) = 3.002e-3;
geomParam(11:12) = [1,1];
if inclAltPhos
    pPSOAlt = pPSO;
    % pPSOAlt(75) = 1.2*pSol(75)*p0(75); % myo phosphate
    geomParamAlt = geomParam;
    % geomParamAlt(8:9) = [0.02, 0.02];
    % geomParamAlt(10) = 10*geomParam(10);
    geomParamAlt(11) = 0.1;
end
phosphateAccum = true;
SSPhosAccum = false;
for i = 1:length(freq)
    % expt 1: HIIT stim, with SOCE, expt 2: HIIT stim, no SOCE
    if i == 1
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
            computeSol(pPSO, [], tSol, 5, freq(i), phosphateAccum, geomParam, SSPhosAccum, withMito);
        pPSO(12) = pPSO(12)*yinits{1}(2); % define cref in terms of resting SR calcium without SOCE
        [Time_noSOCE,Y_noSOCE,~,~,~,yinitNoSOCE] =...
            computeSol(pPSO, [], tSol, 6, freq(i), phosphateAccum, geomParam, SSPhosAccum, withMito);
        if inclAltPhos
            [Time_noSOCEAlt,Y_noSOCEAlt,~,~,~,yinitNoSOCEAlt] =...
                computeSol(pPSOAlt, [], tSol, 6, freq(i), phosphateAccum, geomParamAlt, SSPhosAccum, withMito);
            pPSOAlt(12) = pPSOAlt(12)*yinitNoSOCEAlt{1}(2); % define cref in terms of resting SR calcium without SOCE
        end
    else
        [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge] =...
            computeSol(pPSO, yinits{2}, tSol, 5, freq(i), phosphateAccum, geomParam, SSPhosAccum, withMito);
        [Time_noSOCE,Y_noSOCE] =...
            computeSol(pPSO, yinitNoSOCE{1}, tSol, 6, freq(i), phosphateAccum, geomParam, SSPhosAccum, withMito);
        if inclAltPhos
            [Time_noSOCEAlt,Y_noSOCEAlt] =...
                computeSol(pPSOAlt, yinitNoSOCEAlt{1}, tSol, 6, freq(i), phosphateAccum, geomParamAlt, SSPhosAccum, withMito);
        end
    end
   
    % create plots
    if any(freq(i)==[20, 60])
        figure
        colors = colororder;
        subplot(3,1,1)
        % final entry in color spec gives the opacity
        plot(Time_noSOCE, Y_noSOCE(:,CaJuncIdx)*JFrac + Y_noSOCE(:,CaBulkIdx)*BFrac)
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
        if inclAltPhos
            plot(Time_noSOCEAlt, Y_noSOCEAlt(:,CaJuncIdx)*JFrac + Y_noSOCEAlt(:,CaBulkIdx)*BFrac)
        end
        ylabel('Myo calcium (µM)')
        % title('Cytosolic calcium')
        legend('No SOCE', '20x SOCE', 'No SOCE Alt')
        prettyGraph
        subplot(3,1,2)
        plot(Time_noSOCE, Y_noSOCE(:,CaSRJuncIdx))
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx))
        if inclAltPhos
            plot(Time_noSOCEAlt, Y_noSOCEAlt(:,CaSRJuncIdx))
        end
        ylabel('SR calcium (µM)')
        prettyGraph
        % title('SR calcium')
        subplot(3,1,3)
        plot(Time_noSOCE, Y_noSOCE(:,crossbridgeIdx)/maxCrossBridge)
        hold on
        plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx)/maxCrossBridge)
        if inclAltPhos
            plot(Time_noSOCEAlt, Y_noSOCEAlt(:,crossbridgeIdx)/maxCrossBridge)
        end
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
hold on
subplot(1,2,2)
errorbar(freqExp, maxForceExpMean_withSOCE*169.4, maxForceExpSEM_withSOCE*169.4)
hold on
errorbar(freqExp, maxForceExpMean_noSOCE*169.4, maxForceExpSEM_noSOCE*169.4)
ylim([0 200])
xlim([0 200])
prettyGraph
ylabel('Specific force (mN/mm2)')
xlabel('Freq (Hz)')

%% plot SOCE vs no SOCE over time for fatiguing conditions (Fig 7)
exercise = 'HIIT'; % set this to either 'HIIT' or 'resistance'
pPSO = pSol(:).*p0(:);
% Crossbridge_Cycle = [55:57,59:67,73,84,93,94,106];
% pPSO(Crossbridge_Cycle) = pPSO(Crossbridge_Cycle) .* popParam(2,:)';
pPSO(12) = 0.5;%4*pSol(12)*p0(12); % set cref ratio
soceVals = 10.^[-.5,.5,1.5]*p0(95);%pSol(95)*p0(95); % gSOCE
pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion
phosphateAccum = true;
pPSO(59) = 3*pSol(59)*p0(59); % kon1 (trop binding)
pPSO(75) = 1*pSol(75)*p0(75); % myo phosphate
pPSO(66) = 1*pSol(66)*p0(66); % phos-mediated reverse rate
% phosVec = [0,0,0];

SSPhosAccum = false;
volFracMito = 0.05;
geomParam(8:9) = volFracMito;
geomParam(10) = 1*3.002e-3;
geomParam(11:12) = [1,1];
mitoFactor = [3,3,3];
tCell = cell(length(soceVals),1);
yCell = cell(length(soceVals),1);
switch exercise
    case 'HIIT'
        exerciseParam = [0.65, 0.25*0.65]; % sprinting
        freq = 120;
        tSol = [0, 20];
        % tSol = [0, 60]; % 20 on, 20 off, 20 on
    case 'resistance'
        exerciseParam = [6, 3]; % resistance
        freq = 40;
        tSol = [0, 60];
        % tSol = [0, 150]; % 60 on, 30 off, 60 on
    otherwise
        error('exercise "%s" not recognized, must be either "HIIT" or "resistance"',exercise)
end
withMito = true;
[juncLocLogic, bulkLocLogic] = getLocLogic(withMito);
[CaJuncIdx,CaBulkIdx] = getJuncBulkIdx(8, juncLocLogic, bulkLocLogic);
[CaSRJuncIdx,CaSRBulkIdx] = getJuncBulkIdx(2, juncLocLogic, bulkLocLogic);
[~,crossbridgeIdx] = getJuncBulkIdx(21, juncLocLogic, bulkLocLogic);
[phosJuncIdx,phosBulkIdx] = getJuncBulkIdx(26, juncLocLogic, bulkLocLogic);
[CaPhosJuncIdx,CaPhosBulkIdx] = getJuncBulkIdx(25, juncLocLogic, bulkLocLogic);
figure
for i = 1:length(tCell)
    % pPSO(75) = phosVec(i)*pSol(75)*p0(75); % myo phosphate
    geomParamCur = geomParam;
    geomParamCur(11) = mitoFactor(i);
    pPSO(95) = soceVals(i);
    [Time_withSOCE,Y_withSOCE,fluxes,currents,maxCrossBridge, yinits] =...
        computeSol(pPSO, [], tSol, [3, exerciseParam], freq, phosphateAccum, geomParamCur, SSPhosAccum, withMito);
    tCell{i} = Time_withSOCE;
    yCell{i} = Y_withSOCE;
    subplot(6,1,1)
    plot(Time_withSOCE, Y_withSOCE(:,CaJuncIdx)*JFrac + Y_withSOCE(:,CaBulkIdx)*BFrac)
    hold on
    ylabel('Myo calcium (μM)')
    prettyGraph
    subplot(6,1,2)
    plot(Time_withSOCE, Y_withSOCE(:,CaSRJuncIdx)*JSRFrac + Y_withSOCE(:,CaSRBulkIdx)*BSRFrac)
    hold on
    ylabel('SR calcium (μM)')
    prettyGraph
    subplot(6,1,3)
    maxCrossBridge
    plot(Time_withSOCE, Y_withSOCE(:,crossbridgeIdx))%/maxCrossBridge)
    hold on
    ylabel('Rel force')
    prettyGraph
    subplot(6,1,4)
    plot(Time_withSOCE, Y_withSOCE(:,phosJuncIdx)*JFrac + Y_withSOCE(:,phosBulkIdx)*BFrac)
    hold on
    ylabel('Myo phosphate (μM)')
    prettyGraph
    subplot(6,1,5)
    plot(Time_withSOCE, Y_withSOCE(:,CaPhosJuncIdx)*JSRFrac + Y_withSOCE(:,CaPhosBulkIdx)*BSRFrac)
    hold on
    ylabel('SR calcium phosphate (μM)')
    xlabel('Time (s)')
    prettyGraph
    subplot(6,1,6)
    plot(Time_withSOCE, Y_withSOCE(:,ATPJuncIdx)*JFrac + Y_withSOCE(:,ATPBulkIdx)*BFrac)
    hold on
    ylabel('Myoplasmic ATP (μM)')
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
    phosVec = yCell{i}(:,phosJuncIdx)*JFrac + yCell{i}(:,phosBulkIdx)*BFrac;
    subplot(3,1,2)
    plot(tCell{i},phosVec)
    hold on
    % xlim([0 20])
    prettyGraph
    ylabel('Myo phosphate (μM)')
    forceVec = yCell{i}(:,crossbridgeIdx);
    subplot(3,1,3)
    % plot(tCell{i},cumtrapz(tCell{i},forceVec))
    plot(tCell{i},forceVec/maxCrossBridge)
    hold on
    % [tMaxes,forceMaxes] = getMaxes(tCell{i},forceVec/maxCrossBridge,1/exerciseParam(1));
    % plot(tMaxes, forceMaxes)
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

%% Isolines for phase diagrams
soceFactor = logspace(-1,2,16);
crefFactor = logspace(-0.5,0.5,11);%linspace(0,10,21);
[soceFactor, crefFactor] = meshgrid(soceFactor, crefFactor);
Jrel = logspace(-6.01,-1.01,10);
cSR = 500;
figure
hold on
for i = 1:length(Jrel)
    soceFactorCur = logspace(-1,2,32);
    JrelCur = Jrel(i)./soceFactorCur;
    crefFactorCur = (JrelCur./(1-JrelCur)).^(1/4);
    plot(log10(soceFactorCur), log10(crefFactorCur))
end

%% Phase diagrams for Fig 7
exercise = 'RESISTANCE';%, 'RUNNING', 'RUNNING'};
% exerciseNames = {'RESISTANCE', 'RESISTANCE', 'RESISTANCE'};
switch exercise
    case 'RUNNING'
        freqVals = 120;%[70,120, 70, 120];
        konFactorVals = [3,3,1,1];
    case 'RESISTANCE'
        freqVals = 40;%[40,70, 40, 70];
        konFactorVals = [3,3,1,1];
end
mitoFactor = 4.0;
PiCFactor = 1.0;
% konFactor = 1.0;
crefFactor = 2.0;
figure
turboMap = colormap('turbo');
colormap(turboMap(10:240,:))
numRows = ceil(length(freqVals)/2);
for k = 1:length(freqVals)
    konFactor = konFactorVals(k);
    subplot(numRows,numRows,k)%1,length(exerciseNames),k)
    % exercise = exerciseNames{k};
    % nameStart = sprintf("sweepSOCEWithMito_%.1fmito_%.1fPiC_%.1fTrop_%s_freq%d",...
    %             mitoFactor, PiCFactor, konFactor, exerciseNames{k}, freqVals(k));
    % curFile = sprintf('sweepSOCE_0.1mito_0.1Trop_%s_24-Jun-2026.mat',exercise);
    nameStart = sprintf("sweepSOCEWithMitoFixed_%.1fTrop_%.1fcrefNotCal_%s_freq%d_",...
                konFactor, crefFactor, exercise, freqVals(k));
    try 
        % load(sprintf("%s%s.mat", nameStart, "08-Jul-2026"),...
        %     'results','freq','soceFactor','crefFactor','tLims')
         load(sprintf("%s%s.mat", nameStart, "14-Jul-2026"),...
            'results','freq','soceFactor','mitoFactor','tLims')
    catch
        % load(sprintf("%s%s.mat", nameStart, "03-Jul-2026"),...
        % 'results','freq','soceFactor','crefFactor','tLims')
        load(sprintf("%s%s.mat", nameStart, "011118-Jul-2026"),...
        'results','freq','soceFactor','mitoFactor','tLims')
    end
    idx1 = find(soceFactor==1);
    maxCrossBridge = results{freq==freqVals(k)}{1};
    dynCell = results{freq==freqVals(k)}{end-1};
    finalForceVals = zeros(size(dynCell));
    integratedForceVals = results{freq==freqVals(k)}{end};
    % if contains(exercise,'RUNNING')
    %     testIdx = find(tLims<20,1,'last');
    % else
    testIdx = length(dynCell{1});
    for idx = 1:length(finalForceVals(:))
        finalForceVals(idx) = mean(dynCell{idx}(testIdx))/maxCrossBridge(idx);
        % integratedForceVals(idx) = sum(dynCell{idx}(tLims<20));
    end
    % normIdx = find(soceFactor(:) == 1 & crefFactor(:) == 1);
    % contourf(log10(soceFactor), log10(mitoFactor), finalForceVals)
    contourf(log10(soceFactor), log10(mitoFactor), finalForceVals,20)
    % surf(log10(soceFactor), log10(mitoFactor), finalForceVals, 'FaceColor', 'interp')
    % surf(log10(soceFactor), log10(crefFactor), log10(finalForceVals./finalForceVals(normIdx)), 'FaceColor', 'interp')
    % surf(log10(soceFactor), log10(crefFactor), log10(integratedForceVals), 'FaceColor', 'interp')
    xlabel('Log-fold SOCE conductance')
    ylabel('Log-fold PiC activity')
    colorbar
    view([0 90])
    title(exercise)
    % switch exercise
    %     case 'RUNNING'
    %         clim([0.48,0.72]) 
    %     case 'RESISTANCE'
    %         clim([.25 .8])
    % end
    prettyGraph
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
        withMito = false;
    elseif length(varargin)==2
        geomParam = varargin{1};
        SSPhosAccum = varargin{2};
        withMito = false;
    elseif length(varargin)==3
        geomParam = varargin{1};
        SSPhosAccum = varargin{2};
        withMito = varargin{3};
    elseif isempty(varargin)
        geomParam = [0.01, 0.95, 0.05, 0.003, 20, 1/10, 0.5, 0.02, 0.02, 3.002e-3];
        SSPhosAccum = false;
        withMito = false;
    else
        error('Unknown argument(s)')
    end
    
    if withMito
        juncLocLogic = true(1,40);
        bulkLocLogic = true(1,40);
    else
        juncLocLogic = true(1,31);
        bulkLocLogic = true(1,31);
    end
    juncLocLogic(17:21) = false; % cross bridges
    bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
    if isempty(yinit)
        % load in initial condition starting estimate
        load Data/yinit0.mat yinit0
        yinit0([24,26]) = pPSO(74:75);
        if yinit0(31) > pPSO(96)
            yinit0(31) = pPSO(96);
        end
        if withMito
            yinit0(23) = pPSO(106);
            % yinit0(16) = 1000;
            % yinit0(22) = 1000;
            yinit0 = [yinit0; 100; 0.1; 5000; 5000; 50; 150; 0.0; 0.0; 0.0];
        end
        yinit = [yinit0(juncLocLogic); yinit0(bulkLocLogic)];
        
        % compute steady state solution without SOCE
        pPSO0 = pPSO;
        pPSO0(95) = 0; % no SOCE
        tSS = 0:5:5000;
        if withMito
            [t,y,~,fluxes,currents,ySSFinal_noSOCE] = SkelMuscleCa_dydt_withMito(tSS,0, yinit,...
                pPSO0, tic, 2, SSPhosAccum, geomParam);
        else
            [t,y,~,fluxes,currents,ySSFinal_noSOCE] = SkelMuscleCa_dydt(tSS,0, yinit,...
                pPSO0, tic, 2, SSPhosAccum, geomParam);
        end
        pPSO(12) = pPSO(12)*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        if any(expt(1) == [2,4,6]) || pPSO(95) == 0
            yinf = ySSFinal_noSOCE;
            yinits = {yinf};
        else
            if withMito
                [~,~,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt_withMito(tSS, 0, yinit,...
                    pPSO, tic, 1, SSPhosAccum, geomParam);
            else
                [~,~,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt(tSS, 0, yinit,...
                    pPSO, tic, 1, SSPhosAccum, geomParam);
            end
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
    % yinit_maxCa([21,50]) = pPSO(75); % assuming no phos accum
    yinit_maxCa(sum(juncLocLogic(1:26))) = pPSO(75); % assuming no phos accum
    yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:26))) = pPSO(75);
    yinit_maxCa(sum(juncLocLogic(1:8))) = 1e4;
    yinit_maxCa(sum(juncLocLogic)+sum(bulkLocLogic(1:8))) = 1e4;
    if withMito
        [~, Y_maxCa] = SkelMuscleCa_dydt_withMito([0 1], 0, yinit_maxCa,...
            pPSO, tic, 10, phosphateAccum, geomParam);
    else
        [~, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa,...
            pPSO, tic, 10, phosphateAccum, geomParam);
    end
    crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
    maxCrossBridge = Y_maxCa(end,crossbridgeIdx);
    
    % compute time-dependent solution
    % yinf(sum(juncLocLogic(1:2))) = 1000;
    % yinf(sum(juncLocLogic) + sum(bulkLocLogic(1:2))) = 1000;
    % yinf(sum(juncLocLogic(1:25))) = 0.0;
    % yinf(sum(juncLocLogic(:))+sum(bulkLocLogic(1:25))) = 0.0;
    if withMito
        [t,y,~,fluxes,currents] = SkelMuscleCa_dydt_withMito(tSol, freq, yinf, pPSO,...
            tic, expt, phosphateAccum, geomParam);
    else
        [t,y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, yinf, pPSO,...
            tic, expt, phosphateAccum, geomParam);
    end
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
    tMax = (tMax(1:end-1)+tMax(2:end))/2;
end

function [juncIdx,bulkIdx] = getJuncBulkIdx(idx, juncLocLogic, bulkLocLogic)
    juncIdx = sum(juncLocLogic(1:idx));
    bulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:idx));
end

function [juncLocLogic, bulkLocLogic] = getLocLogic(withMito)
    juncLocLogic = true(1,40);
    juncLocLogic(17:21) = false; % cross bridges
    bulkLocLogic = true(1,40);
    bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
    if ~withMito
        juncLocLogic = juncLocLogic(1:31);
        bulkLocLogic = bulkLocLogic(1:31);
    end
end

function [JFrac, BFrac, JSRFrac, BSRFrac, JFracMito, BFracMito] =...
                getFracs(volFraction_mitoJ, volFraction_mitoB, diffusion_length)
    R_fiber = 20;
    vol_Fiber = pi * (R_fiber ^ 2) * 100 ;
    vol_SA_ratio = 0.01;
    volFraction_TT = 0.003;
    volFraction_SR = 0.05;
    volFraction_myo = 1-volFraction_TT-volFraction_SR;
    vol_myo = volFraction_myo * vol_Fiber ;
    SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
    vol_myoJ = SA_TT*diffusion_length;
    vol_myo_corrected = vol_myoJ*(1-volFraction_mitoJ) + (vol_myo-vol_myoJ)*(1-volFraction_mitoB);
    JFrac = vol_myoJ*(1-volFraction_mitoJ) / vol_myo_corrected;
    BFrac = 1 - JFrac;
    if volFraction_mitoJ == 0
        JFracMito = 0.0;
    else
        JFracMito = vol_myoJ*volFraction_mitoJ/(vol_myo-vol_myo_corrected);
    end
    if volFraction_mitoB == 0
        BFracMito = 0.0;
    else
        BFracMito = 1 - JFracMito;
    end
    vol_SR = volFraction_SR*vol_Fiber;
    SRJ_occupancy = 0.5;
    SA_SRJ = SA_TT * SRJ_occupancy;
    vol_SRJ = SA_SRJ*diffusion_length;
    JSRFrac = vol_SRJ / vol_SR;
    BSRFrac = 1 - JSRFrac;
end