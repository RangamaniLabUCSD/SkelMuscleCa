load Data/p0Struct.mat p0Struct
p0 = p0Struct.data;
% p0(106) = 700;
% load Data/pSol_allParam.mat pSol
% load PSO_22-Jun-2026_newRange.mat pSol
load pSol_fullWithMito/PSO_14-Jul-2026_withMito.mat pSol
% pSol = pVec;
VOnlyIdx = [1,4,5,6,9,11,13,14,16,18,19,22,23,24,25,26,30,33,40,47,76,77,79,80,81,82,98];
% VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
highSensIdx = 1:106;
if length(highSensIdx) < length(p0) % only applicable is highSensIdx is restricted to not all of p0
    VOnlyStruct = load('Data/pVec_VOnlyNew.mat', 'pVec');
    pVecVOnly = VOnlyStruct.pVec;
    pVecVOnly = pVecVOnly(:); % be sure it is a column vector
    [VOnlyOnly,onlyIdx] = setdiff(VOnlyIdx, highSensIdx); % non overlapping indices
    p0(VOnlyOnly) = pVecVOnly(onlyIdx).*p0(VOnlyOnly); % set p0 according to previous estimation
end
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);
phosphateAccum = true;
withMito = true;
juncLocLogic = true(1,40);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,40);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
if ~withMito
    juncLocLogic = juncLocLogic(1:31);
    bulkLocLogic = bulkLocLogic(1:31);
end
crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
CaJuncIdx = sum(juncLocLogic(1:8));
CaBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:8));
CaSRJuncIdx = sum(juncLocLogic(1:2));
CaSRBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:2));
CaECJuncIdx = sum(juncLocLogic(1:27));
CaMitoJuncIdx = sum(juncLocLogic(1:33));

% define geometric quantities
R_fiber = 20;
vol_SA_SR = 1/10;
vol_Fiber = pi * (R_fiber ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003;
volFraction_SR = 0.05;
volFraction_myo = 1-volFraction_TT-volFraction_SR;
vol_myo = volFraction_myo * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
diffusion_length = pPSO(99);
vol_myoJ = SA_TT*diffusion_length;
if withMito
    volFraction_mitoJ = 0.05; 
    volFraction_mitoB = 0.05;
else
    volFraction_mitoJ = 0.0;
    volFraction_mitoB = 0.0;
end
vol_myo_corrected = vol_myoJ*(1-volFraction_mitoJ) + (vol_myo-vol_myoJ)*(1-volFraction_mitoB);
JFrac = vol_myoJ*(1-volFraction_mitoJ) / vol_myo_corrected;
BFrac = 1 - JFrac;
JFracMito = vol_myoJ*volFraction_mitoJ/(vol_myo-vol_myo_corrected);
BFracMito = 1 - JFracMito;
vol_SR = volFraction_SR*vol_Fiber;
SRJ_occupancy = 0.5;
SA_SRJ = SA_TT * SRJ_occupancy;
vol_SRJ = SA_SRJ*diffusion_length;
JSRFrac = vol_SRJ / vol_SR;
BSRFrac = 1 - JSRFrac;
geomParam = [vol_SA_ratio, volFraction_myo, volFraction_SR,...
             volFraction_TT, R_fiber, vol_SA_SR, SRJ_occupancy,...
             volFraction_mitoJ,volFraction_mitoB,3.002e-3];

username = getenv('USER');
foldername = sprintf('/tscc/lustre/ddn/scratch/%s/SkelMuscleCaData',username);
if isfolder(foldername) % this folder must exist to run on server
    onTSCC = true;
    parpool(60)
else % running locally, set smaller parpool
    foldername = "Data";
    onTSCC = false;
    % parpool(60)
end

soceFactor = logspace(-1,2,16);
% crefFactor = logspace(-0.5,0.5,11);%linspace(0,10,21);
mitoFactor = logspace(-1,1,11);
[soceFactor, mitoFactor] = meshgrid(soceFactor, mitoFactor);
zeroMat = zeros(size(soceFactor));
maxForce = zeroMat; peakForce = zeroMat; finalForce = zeroMat;
integratedForce = zeroMat;
cMax = zeroMat; cFinal = zeroMat;
cSRMax = zeroMat; cSRFinal = zeroMat;

pPSO(100) = 1*pSol(100)*p0(100); % ion diffusion 

exerciseParams = {[0.65, .25*.65], [6, 3]}; 
freqs = {[120,70], [40,70]};
exerciseNames = {'RUNNING', 'RESISTANCE'};

konFactor = [1,0.3];%[3, 3, 1, 1, 0.3, 0.3];
% mitoFactor = [0.01, 0.05, 0.01, 0.05];
% mitoFactor = [1, 4, 1, 4, 1, 4, 1, 4];
% PiCFactor = [1,1,1,1,5,5,5,5];
% MCUFactor = [1,1,1,1,1,1,1,1];
crefFactor = [2,2];%[2, 1, 2, 1, 2, 1];

for outeridx = 1:length(konFactor)
    pPSO(12) = crefFactor(outeridx)*0.25;%pSol(12)*p0(12); % set cref ratio
    pPSO(59) = konFactor(outeridx)*pSol(59)*p0(59); % kon1 (trop ca binding)
    % pPSO(72) = phosDegFactor(outeridx)*pSol(72)*p0(72); % phosphate degrad
    % geomParam(8:9) = mitoFactor(outeridx);
    % geomParam(11:12) = [PiCFactor(outeridx), MCUFactor(outeridx)];
    geomParam(11:12) = [1,1];

    % store reference for parfor loop
    pPSORef = pPSO;
    geomParamRef = geomParam;

    for exIdx = 1:length(exerciseParams)
    
        exerciseParam = exerciseParams{exIdx}; % period, activity width
        if strcmp(exerciseNames{exIdx}, 'RUNNING')
            tSol = [0,20];
        elseif strcmp(exerciseNames{exIdx}, 'RESISTANCE')
            tSol = [0,60];
        end
        tLims = 0:exerciseParam(1):tSol(end);
        if (tSol(end)-tLims(end)) > exerciseParam(2)
            tLims = [tLims, tSol(end)];
        end
        freq = freqs{exIdx};
        results = cell(size(freq));
        for i = 1:length(freq)
            nameStart = sprintf("sweepSOCEWithMitoFixed_%.1fTrop_%.1fcrefNotCal_%s_freq%d_",...
                konFactor(outeridx), crefFactor(outeridx), exerciseNames{exIdx}, freq(i));
            filename = nameStart + string(datetime("today")) +".mat";
            dynCell = cell(size(soceFactor));
            parfor idx = 1:length(soceFactor(:))
                pPSO = pPSORef;
                geomParamCur = geomParamRef;
                pPSO(95) = soceFactor(idx)*p0(95);%#ok<PFBNS> pSol(95)*p0(95); %#ok<PFBNS> % gSOCE
                geomParamCur(11) = mitoFactor(idx);
                % [tCur,yCur,fluxes,currents,maxCrossBridge, yinits] =...
                %     computeSol(pPSO, [], tSol, [3,exerciseParam], freq(i), phosphateAccum); %#ok<PFBNS>
                [tCur,yCur,fluxes,currents,maxCrossBridge, yinits] =...
                    computeSolNew(pPSO, [], tSol, [3,exerciseParam], freq(i), phosphateAccum, geomParamCur, false, withMito); %#ok<PFBNS>
                % collect stats
                maxForce(idx) = maxCrossBridge;
                peakForce(idx) = max(yCur(:,crossbridgeIdx));%/maxCrossBridge);
                finalForce(idx) = yCur(end,crossbridgeIdx);
                integratedForce(idx) = trapz(tCur,yCur(:,crossbridgeIdx));
                cCur = yCur(:,CaJuncIdx)*JFrac + yCur(:,CaBulkIdx)*BFrac;
                cSRCur = yCur(:,CaSRJuncIdx)*JSRFrac + yCur(:,CaSRBulkIdx)*BSRFrac;
                cMax(idx) = max(cCur);
                cFinal(idx) = cCur(end);
                cSRMax(idx) = max(cSRCur);
                cSRFinal(idx) = cSRCur(end);
                dynVals = zeros(length(tLims)-1,1);
                for k = 2:length(tLims)
                    curLogic = tCur > tLims(k-1) & tCur <= tLims(k);
                    dynVals(k-1) = max(yCur(curLogic,crossbridgeIdx));
                end
                dynCell{idx} = dynVals;
                fprintf('Done with case %d of %d for freq %.2f, case %s\r\n',...
                    idx, length(soceFactor(:)), freq(i), nameStart) %#ok<PFBNS>
            end
            results{i} = {maxForce, peakForce, finalForce, cMax, cFinal,...
                          cSRMax, cSRFinal, dynCell, integratedForce};
            save(fullfile(foldername,filename));
        end
    end
end


function [t,y,fluxes,currents,maxCrossBridge,yinits] =...
                computeSol(pPSO, yinit, tSol, expt, freq, phosphateAccum)
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
    % Outputs:
    %   - t is the vector of times
    %   - y is the vector of state variables
    %   - fluxes: calcium fluxes at each time point, each row consists:
    %        [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, 
    %         LumpedJ_SERCA, J_CaLeak_SR, Jhydtot, total calcium flux]
    %        in units of µM/s for junctional myoplasm, then the same
    %        entities for the bulk myoplasm
    %   - currents: Ionic and total current at each time point, each row consists of:
    %            [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N,
    %             I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL] in units of pA
    %             for the T-tubules and then the same entities for the SL
    %   - maxCrossBridge: post power stroke cross bridge density at
    %     saturating calcium concentration for these parameters
    %   - yinits: Cell vector with each entry a vector containing the initial values
    %     for all 51 state variables in the model.
    %     Cell length 1 if only SOCE-free condition is tested, or if
    %     initial condition is provided.
    %     Otherwise first cell entry is for case without SOCE and second
    %     entry is for case with SOCE.

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
        [~,~,~,~,~,ySSFinal_noSOCE] = SkelMuscleCa_dydt([0 1000],0, yinit, pPSO0, tic, 2, false);%phosphateAccum);
        pPSO(12) = pPSO(12)*ySSFinal_noSOCE(2); % set c_ref according to SS c_SR
        if any(expt(1) == [2,8]) || pPSO(95) == 0
            yinf = ySSFinal_noSOCE;
            yinits = {yinf};
        else
            [~,~,~,~,~,ySSFinal_withSOCE] = SkelMuscleCa_dydt([0 1000], 0, yinit, pPSO, tic, 1, false);%phosphateAccum);
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
    [~, Y_maxCa] = SkelMuscleCa_dydt([0 1], 0, yinit_maxCa, pPSO, tic, 10, phosphateAccum);
    crossbridgeIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
    maxCrossBridge = Y_maxCa(end,crossbridgeIdx);
    
    % test single frequency dynamics
    [t,y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, yinf, pPSO, tic, expt, phosphateAccum); % compute time-dependent solution
end

function [t,y,fluxes,currents,maxCrossBridge,yinits] = ...
    computeSolNew(pPSO, yinit, tSol, expt, freq, phosphateAccum, varargin)
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
            % yinit0 = [yinit0; 100; 0.1; 10000; 100; 50; 150; 0.0; 0; 0];
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
        if any(expt(1) == [2,4]) || pPSO(95) == 0
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
    if withMito
        [t,y,~,fluxes,currents] = SkelMuscleCa_dydt_withMito(tSol, freq, yinf, pPSO,...
            tic, expt, phosphateAccum, geomParam);
    else
        [t,y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, yinf, pPSO,...
            tic, expt, phosphateAccum, geomParam);
    end
end
