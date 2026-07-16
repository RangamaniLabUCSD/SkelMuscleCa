function [Time,Y,currtime,fluxes,currents,ySSFinal] =... 
    SkelMuscleCa_dydt_withMito(tSpan,freq, yinit, p,StartTimer,expt,phosphateAccum,varargin)
% inputs:
%     - tSpan is a vector of start and stop times
%       freq is a vector of test frequencies
%     - yinit is a vector of initial conditions for the state variables
%     - p is a vector of selected parameters to test
%     - StartTimer starts counting run time
%     - expt is the experimental value used for calculation. Options:
%           (-2): estimation, Rincon et al 2021, Fig 1B calcium data (5 peaks 100 Hz, IIb muscle)
%           (-3): estimation, Baylor and Hollingworth 2003, Fig 2A (fast twitch curve)
%           (-8): estimation, Miranda et al 2020, Fig 2B (WT) membrane voltage data
%           (1): Standard conditions with SOCE
%           (2): Standard conditions without SOCE
%           (3): Exercise test with SOCE
%           (4): Exercise test without SOCE
%           (10): Crossbridge-only test
%         In cases (3) and (4), expt may be a vector with three elements,
%         [expt_n, cycle_time, stim_time], where expt_n is 3 or 4,
%         cycle_time is the total time per repitition/stride and stim_time
%         is the time of activation for each repitition/stride. Default
%         cycle_time and stim_time if not specified are 6 s and 3 s.
%     - phosphateAccum is a logical variable determining if phosphate
%       accumulation is accounted for
% optional inputs (stored in varargin):
%     - geomParam (varargin{1}): vector with stored geometric parameters.
%     If not provided, geomParam is assigned default values (see code
%     below)
%
% outputs:
%     - T is the vector of times
%     - Y is the vector of state variables
%     - currtime is the total runtime
%     - fluxes: calcium fluxes at each time point, each row consists:
%            [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, 
%             LumpedJ_SERCA, J_CaLeak_SR, Jhydtot, total calcium flux]
%            in units of µM/s for junctional myoplasm, then the same 
%            entities for the bulk myoplasm
%     - currents: Ionic and total current at each time point, each row consists of:
%            [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N,
%             I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL] in units of pA
%             for the T-tubules and then the same entities for the SL
%     - ySSFinal: Only returns in case of steady-state (SS) estimation -
%       variable values associated with minimum dydt over tested simulation
%       (handles cases of oscillatory y)
% -------------------------------------------------------------------------

% first load geometric parameters if they are provided as an optional arg
if isscalar(varargin) % length 1, that is
    geomParam = varargin{1};
elseif isempty(varargin)
    geomParam = [0.01, 0.95, 0.05, 0.003, 20, 1/10, 0.5, 0.05, 0.05, 3.002e-3];
else
    error('Unknown argument(s)')
end

if length(expt) > 1
    if any(expt(1) == [3,4]) && length(expt) == 3
        exerciseParam = expt(2:3);
        expt = expt(1);
    else
        error('Unrecognized extra expt arguments')
    end
elseif any(expt == [3,4])
    exerciseParam = [6, 3];
end
        
% each compartment has junctional/terminal portion and bulk portion
% extracellular: TTvol and EC
% plasma membrane: TTM and SL
% myoplasm: JM and BM
% SRM: SRJM and SRBM
% SR: SRJ and SRB
% juncLocLogic stores whether each variable is present in the
% junctional portion of a compartment
% bulkLocLogic stores whether each variable is present in the bulk
% portion of a compartment

juncLocLogic = true(1,40);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,40);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions

% ode15s settings
% pull out non-negative variables (everything except voltage)
nonNegJ = [1:sum(juncLocLogic(1:4)),sum(juncLocLogic(1:6)):sum(juncLocLogic(1:36)),sum(juncLocLogic(1:38)):sum(juncLocLogic(:))];
nonNegB = [1:sum(bulkLocLogic(1:4)),sum(bulkLocLogic(1:6)):sum(bulkLocLogic(1:36)),sum(bulkLocLogic(1:38)):sum(bulkLocLogic(:))];
nonNegB = nonNegB + nonNegJ(end);
if freq == 0 && ~any(expt==[7,8]) % Steady state condition
    options = odeset('RelTol',1e-3,'MaxStep', 5.0, 'NonNegative',[nonNegJ, nonNegB]);
elseif freq == 0 && any(expt==[7,8]) % SOCE test without extra stimulus
    options = odeset('RelTol',1e-6,'MaxStep',0.1,'NonNegative', [nonNegJ, nonNegB]);
else
    options = odeset('RelTol',1e-6,'MaxStep',min([.001,0.1/freq]),'NonNegative', [nonNegJ, nonNegB]);
end
if expt < 0 % then corresponds to a case of estimation
    isEst = true;
    expt = abs(expt);
    %Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C]
            %[K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_t MJ_V],};
    Ca_o_exp = [1000, 1000, 2000, 2000, 2000,...
                2500, 1800, 5000, 2000, 2000, 1300];                                      %uM
    Na_o_exp = [138100, 138100, 150000, 150000, 150000,...
                 143800, 118240, 140000, 151000, 151000, 147000];                           %uM
    K_o_exp = [3900, 3900, 2000, 2000, 2000,...
               5000, 5330, 4000, 5000, 5000, 4000];                               %uM
    Cl_o_exp = [143700, 143700, 158000, 158000, 158000,...
                124000, 126170, 157000, 146000, 146000, 127400];                           %uM
    Temp = [(273.15+22),(273.15+22),(273.15+20),(273.15+22),(273.15+22),...
            (273.15+26),(273.15+22),(273.15+22),(273.15+35),(273.15+22),(273.15+30)]; %K
    T = Temp(expt);        %K
    if expt == 10
        expt_n = 5;
    else
        expt_n = expt;
    end
    Ca_EC = Ca_o_exp(expt_n);
    Na_EC = Na_o_exp(expt_n);
    K_EC = K_o_exp(expt_n);
    Cl_EC = Cl_o_exp(expt_n);
    yinit(sum(juncLocLogic(1:27))) = Ca_EC;    %µM
    yinit(sum(juncLocLogic(1:28))) = Na_EC;   %µM
    yinit(sum(juncLocLogic(1:29))) = K_EC;      %µM
    yinit(sum(juncLocLogic(1:30))) = Cl_EC;   %µM
else
    isEst = false;
    T = 273.15 + 22;%295.15;
    Ca_EC = 1300;
    Na_EC = 147000.0;
    K_EC = 4000.0;
    Cl_EC = 128000.0;
end
% load Data/pMito.mat pMito
% load Data/pMitoNew.mat pMito
% load pSol_fullWithMito/bestMitoParam.mat mitoParams
% pMito = mitoParams;
mitoStruct = load('PSOMito_11-Jul-2026.mat', 'pSol');
pMito = mitoStruct.pSol;
yinit(sum(juncLocLogic(1:38))) = yinit(sum(juncLocLogic(1:33)))*yinit(sum(juncLocLogic(1:35)))/(500*pMito(1));
yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:38))) = yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:33)))...
    *yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:35)))/(500*pMito(1));
yinit(sum(juncLocLogic(1:39))) = pMito(38)*15000 - yinit(sum(juncLocLogic(1:34)));
yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:39))) =...
    pMito(38)*15000 - yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:34)));
yinit(sum(juncLocLogic(1:40))) = pMito(39)*250 - yinit(sum(juncLocLogic(1:36)));
yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:40))) =...
    pMito(39)*250 - yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:36)));
if any(yinit(nonNegJ) < 0) || any(yinit(nonNegB) < 0)
    print('This will throw an error\r\n')
end
% pMito = ones(42,1);

% test integration
% tvals = 0:.00001:1;
% yvals = zeros(length(tvals),length(yinit));
% yvals(1,:) = yinit;
% for i = 1:length(tvals)-1
%     dydtTest = f(tvals(i),yvals(i,:),p,freq);
%     yvals(i+1,:) = yvals(i,:) + .00001*dydtTest;
% end

% call ode15s
[Time,Y] = ode15s(@f,tSpan,yinit,options,p,freq); %pass extra arguments at the end

% return fluxes and currents at select times
StartTimer = tic;
fluxes = zeros(length(Time), 2*13);
currents = zeros(length(Time), 2*14);
% Flux and ionic current through different pathways.
for i = 1:length(Time)
    [~, fluxesCur, currentsCur] = f(Time(i), Y(i,:), p, freq);
    fluxes(i,:) = [fluxesCur(1,:),fluxesCur(2,:)];
    currents(i,:) = [currentsCur(1,:), currentsCur(2,:)];
end

if freq == 0
    ssEst = true;
else
    ssEst = false;
end
if ssEst % find values with min dydt (handles cases of oscillatory y)
    StartTimer = tic;
    tSS = Time(Time > Time(end)/2);
    ySS = Y(Time > Time(end)/2, :);
    yMean = mean(ySS,1);
    yMean = yMean(:);
    dyVals = zeros(size(tSS));
    for i = 1:length(tSS)
        dyCur = f(tSS(i),ySS(i,:),p,freq);
        dyCur = dyCur(:);
        testIdx = find(yMean~=0);
        dyVals(i) = sqrt(sum((dyCur(testIdx)./yMean(testIdx)).^2));
    end
    [~,idx] = min(dyVals);
    ySSFinal = ySS(idx,:);
end


% -------------------------------------------------------
% ode rate
    function [dydt, fluxes, currents] = f(t,y,p,freq)
        currtime = toc(StartTimer);
        if isEst && currtime > 120
            error('too long to compute!') % catches cases of bad parameters that get stuck in long integration
        end

        %% State Variables - junctional and bulk
        Ca_EC_cur = Ca_EC;
        if ~isEst
            if any(expt == [7,8]) && (t > 6*60 && t < 17*60)
                Ca_EC_cur = Ca_EC*(0.001 + 0.999*exp(-(t-6*60)));
            elseif any(expt == [7,8]) && (t > 17*60)
                Ca_EC_cur = Ca_EC*(1-exp(-(t-17*60)));
            end
        end
        yAll = zeros(length(juncLocLogic), 2); % first column - junc variables, second column - bulk variables
        yAll(juncLocLogic,1) = y(1:sum(juncLocLogic));
        yAll(bulkLocLogic,2) = y(sum(juncLocLogic)+1:end);
        yAll(27:30,2) = [Ca_EC_cur; Na_EC; K_EC; Cl_EC];

        %% Global constants
        F = 96485.3321;
        R = 8314.46261815;
        KMOLE = 0.001660538783162726;

        pulsewidth = 0.001;%s

        %% Model Geometry    
        L_fiber = 100;  % arbitrary (cancels out) µm
        vol_SA_ratio = geomParam(1); %0.01 µm
        volFraction_myo = geomParam(2); %0.95;
        volFraction_SR = geomParam(3); %0.05 ;
        volFraction_TT = geomParam(4); %0.003 ;
        R_fiber = geomParam(5); %20 µm
        vol_SA_SR = geomParam(6); %1/10;
        SRJ_occupancy = geomParam(7); %0.5;
        volFraction_mitoJ = geomParam(8);
        volFraction_mitoB = geomParam(9);
        vol_SA_mito = geomParam(10);
        
        volFraction_myo = 1 - volFraction_SR - volFraction_TT; % correct volFraction_myo for consistency
        vol_Fiber = pi * (R_fiber ^ 2) * L_fiber;
        vol_myo = volFraction_myo * vol_Fiber;
        SA_SL = 2 * pi * R_fiber * L_fiber;
        vol_SR = volFraction_SR * vol_Fiber;
        SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
        L_TT = R_fiber;
        
        SA_SR = vol_SR / vol_SA_SR;
        
        SA_SRJ = SA_TT * SRJ_occupancy;
        SA_SRB = SA_SR - SA_SRJ;
        diffusion_length = p(99);%0.05;

        vol_myoJ = SA_TT*diffusion_length;
        vol_myoB = vol_myo - vol_myoJ;
        vol_SRJ = SA_SRJ*diffusion_length;
        vol_SRB = vol_SR - vol_SRJ;
        vol_TT = volFraction_TT * vol_Fiber;
        vol_mitoJ = vol_myoJ*volFraction_mitoJ;
        vol_mitoB = vol_myoB*volFraction_mitoB;
        SA_mitoJ = vol_mitoJ / vol_SA_mito;
        SA_mitoB = vol_mitoB / vol_SA_mito;
        % subtract mito from myo
        vol_myoJ = vol_myoJ - vol_mitoJ;
        vol_myoB = vol_myoB - vol_mitoB;

        volFactor = (vol_myo ./ (pi .* 0.26)); % size of our volume vs size of Senneff volume

        % define SA/vol ratios for junctional vs. bulk
        % [PM/EC, PM/myo, SRM/myo, SRM/lumen, MitoM/myo, MitoM/matrix]
        SAvols = [SA_TT/vol_TT, SA_TT/vol_myoJ, SA_SRJ/vol_myoJ, SA_SRJ/vol_SRJ, SA_mitoJ/vol_myoJ, SA_mitoJ/vol_mitoJ;
                  0.0         , SA_SL/vol_myoB, SA_SRB/vol_myoB, SA_SRB/vol_SRB, SA_mitoB/vol_myoB, SA_mitoB/vol_mitoB];
        vols = [vol_TT, vol_myoJ, vol_SRJ, vol_mitoJ;
                inf,    vol_myoB, vol_SRB, vol_mitoB];
        SAs = [SA_TT, SA_SRJ, SA_mitoJ;
               SA_SL, SA_SRB, SA_mitoB];
        
        %% define diffusive flux rates
        D_ionMyo = 1*p(100);
        D_ionSR = 1*p(100); 
        D_ionEC = 1000;%10*p(100);
        D_ATP = p(101);
        D_parv = p(102);
        D_CSQ = p(103);
        D_Pi = p(104);
        D_V = 1/p(105);
        SA_JBmyo = (1-SRJ_occupancy)*SA_TT;
        SA_TTtop = 2*vol_TT/L_TT;
        D_vec = [0, D_ionSR, 0, 0, D_V, D_ionMyo, D_ionMyo, D_ionMyo, 0, 0,...
                   0, 0, D_ionMyo, D_parv, D_parv, D_ATP, 0, 0, 0, 0,...
                   0, D_ATP, D_ATP, D_Pi, D_Pi, D_Pi, D_ionEC, D_ionEC, D_ionEC, D_ionEC, D_CSQ, D_ATP,...
                   0, 0, 0, 0, 0, 0, 0, 0];
        juncRatios = [0, SA_SRJ/vol_SRJ, 0, 0, SA_SL/SA_TT, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, 0, 0,...
                      0, 0, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, 0, 0, 0, 0,...
                      0, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_SRJ/vol_SRJ, SA_SRJ/vol_SRJ, SA_JBmyo/vol_myoJ,...
                      SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_SRJ/vol_SRJ, SA_JBmyo/vol_myoJ...
                      0, 0, 0, 0, 0, 0, 0, 0];
        bulkRatios = [0, SA_SRJ/vol_SRB, 0, 0, 1.0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, 0, 0,...
                      0, 0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, 0, 0, 0, 0,...
                      0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_SRJ/vol_SRB, SA_SRJ/vol_SRB, SA_JBmyo/vol_myoB,...
                      0, 0, 0, 0, SA_SRJ/vol_SRB, SA_JBmyo/vol_myoB...
                      0, 0, 0, 0, 0, 0, 0, 0];

        %% define different fluxes/currents based on junctional vs. bulk
        % fluxes = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, 
        %           LumpedJ_SERCA, J_CaLeak_SR];
        % currents = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, 
        %             I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
        juncFluxFlags = [1, 1, 1, 1, 1, 1, 0.1, 1];
        bulkFluxFlags = [0, 1, 1, 0, 1, 0, 1, 1];
        testBulkSOCE = false;
        if testBulkSOCE
            bulkFluxFlags(1) = 3.0;
        end
        juncCurrentFlags = [1, 0.1, 1, 0.45, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 0];
        bulkCurrentFlags = [1, 1. , 0, 1.  , 1, 1, 1, 1.0, 1.0 , 1. , 1, 0, 1];
        fluxFlags = {juncFluxFlags, bulkFluxFlags};
        currentFlags = {juncCurrentFlags, bulkCurrentFlags};
        varLogic = {juncLocLogic, bulkLocLogic};
        fluxes = zeros(2,13);
        currents = zeros(2,14);
        dydt = zeros(size(y));
        Rmag = zeros(size(y));
        curStartIdx = 1;
        Nf = [1;1000;1;1;100;1000;1000;0.1;1;1;1;1;100000;500;1000;...
            1;1;1;1;1;1;1;1;1;1;1;1300; 147000.0; 4000.0; 128000.0; 1000; 100;...
            1;1000;1000;100;100;100;1000;100]; %Normalization factor

        % Calculate ADP in the myoplasm and mito
        % ANPtot_myo = ANPtot;
        % ADP = 500;%ANPtot_myo - ATP - MgATP - CATP;

        for k = 1:2 % k = 1: junctional, k = 2: bulk
            fluxFlagsCur = fluxFlags{k};
            currentFlagsCur = currentFlags{k};
            varLogicCur = varLogic{k};
            SAvolsCur = SAvols(k,:);
            volsCur = vols(k,:);
            vol_myo = volsCur(2);
            vol_SR = volsCur(3);
            SAsCur = SAs(k,:);
            SA_SL = SAsCur(1);
            SRMFrac = SAsCur(2) / SA_SR;
            vol_mito = volsCur(4);
            SA_mito = SAsCur(3);
            SA_mito_ref = vol_mito / 3.002e-3;

            %% load variables for current compartment
            SOCEProb = yAll(1,k);
            c_SR = yAll(2,k);
            h_K = yAll(3,k);
            w_RyR = yAll(4,k);
            Voltage_SL = yAll(5,k);
            Na_i = yAll(6,k);
            Cl_i = yAll(7,k);
            c_i = yAll(8,k);
            n = yAll(9,k);
            m = yAll(10,k);
            h = yAll(11,k);
            S = yAll(12,k);
            K_i = yAll(13,k);
            CaParv = yAll(14,k);
            MgParv = yAll(15,k);
            CATP = yAll(16,k);
            CaTrop = yAll(17,k);
            CaCaTrop = yAll(18,k);
            D_2 = yAll(19,k);
            Pre_Pow = yAll(20,k);
            Post_Pow = yAll(21,k);
            MgATP = yAll(22,k);
            ATP = yAll(23,k);
            p_i_SR = yAll(24,k);
            PiCa_SR = yAll(25,k);
            p_i_Myo = yAll(26,k);
            c_o = yAll(27,k);
            Na_o = yAll(28,k);
            K_o = yAll(29,k);
            Cl_o = yAll(30,k);
            CSQ = yAll(31,k);
            ADP = yAll(32,k);
            % if ADP < 0
            %     error("ADP cannot be negative")
            % end
            % mito variables
            c_M = yAll(33,k); % mito calcium
            ATP_M = yAll(34,k);
            Pi_M = yAll(35,k);
            NADH = yAll(36,k);
            V_M = yAll(37,k); % IM potential mito
            PiCa_M = yAll(38,k); % PiCa precipitate in mito
            ADP_M = yAll(39,k);
            NAD = yAll(40,k);

            %% Variable Parameters
            A_a = p(1);
            A_hkinf = p(2);
            A_Sinf = p(3);
            alpha_h0 = p(4);
            alpha_m0 = p(5);
            alpha_n0 = p(6);
            alpha_w_r0 = p(7);
            beta_h0 = p(8);
            beta_m0 = p(9);
            beta_n0 = p(10);
            C_SL = p(11);
            c_ref = p(12);
            ClampCurrent = p(13)*currentFlagsCur(13);
            delta = p(14);
            f_RyR = p(15);
            g_Cl = p(16)*currentFlagsCur(2);
            g_K = p(17)*currentFlagsCur(4);
            G_K = p(18)*currentFlagsCur(5);
            g_Na = p(19)*currentFlagsCur(10);
            g_NCX = p(20)*currentFlagsCur(6);
            J_NaK_NKX = p(21)*currentFlagsCur(8);
            K_alphah = p(22);
            K_alpham = p(23);
            K_alphan = p(24);
            K_betah = p(25);
            K_betam = p(26);
            K_betan = p(27);
            K_K = p(28);
            K_mK_NKX = p(29);
            K_mNa_NKX = p(30);
            K_PMCA = p(31);
            K_RyR = p(32);
            K_S = p(33);
            K_SERCA = p(34);
            K_w_r0 = p(35);
            Kdact_NCX = p(36);
            Kmc_i_NCX = p(37);
            ksat_NCX = p(38);
            nu_NCX = p(39);
            S_i = p(40);
            tau_SOCEProb = p(41);
            nu_SERCA = p(42)*fluxFlagsCur(7);
            g_PMCA = p(43)*fluxFlagsCur(5);
            nu_leakSR = p(44)*fluxFlagsCur(8);
            g_leakNa = p(45)*currentFlagsCur(10);
            k_onATP = p(46);
            k_offATP = p(47);
            k_onParvCa = p(48);
            k_offParvCa = p(49);
            k_onParvMg = p(50);
            k_offParvMg = p(51);
            Parv_itot = p(52);
            k_onMA = p(53);               %(µM/s) Mg2+ binding to ATP
            k_offMA =p(54);
            kHYD =p(55);
            kH = p(56);
            kPROD = p(57);
            Mg = p(58);
            k_onTrop1 = p(59);
            k_offTrop1 = p(60);
            k_onCa = p(61);
            k_offCa = p(62);
            f0 = p(63);
            fP = p(64);
            h0 = p(65);
            hP = p(66);
            g0 = p(67);
            PP = p(68);
            kP = p(69);
            Ap = p(70);
            Bp = p(71);
            bP = 0*p(72);
            Trop_tot = p(73);   
            V_a = p(76);
            V_h = p(77);
            V_hkinf = p(78);
            V_m = p(79);
            V_n = p(80);
            V_Sinf = p(81);
            V_tau = p(82);
            VBar_RyR = p(83);
            kATP = p(84);
            KmNa_i_NCX = p(85);
            KmNa_EC_NCX = p(86);
            Kmc_EC_NCX = p(87);
            % Kmc_EC_NCX_C = p(88);
            g_SLLeak = p(89)*fluxFlagsCur(2);
            L_RyR = p(90);
            g0_DHPR = 0*p(91)*fluxFlagsCur(4);
            j0_RyR = p(92)*fluxFlagsCur(6);
            k_onTrop2 = p(93);
            k_offTrop2 = p(94);   
            B_CSQ = p(96);
            K_CSQ = p(97);
            k_offCSQ = p(98);
            ATPset = p(106);

    
            %% Carrier Valence
            z_Ca = 2;
            z_Na = 1;
            z_K = 1;
            z_Cl = -1;
    
            %% Temperature Coeff
            QCorr = (T-293)/10;
            Q10NCX = 1.57;
            Q10g_K = 1.5;
            Q10g_Na = 1.5;
            Q10g_KIR= 1.55;
            Q10g_NaK = 1 ;
            Q10g_Cl = 1.5;
            Q10C_SL = 1.02;
            Q10alpha_n = 2.5;
            Q10alpha_m = 2.3;
            Q10alpha_h = 2.5;
            Q10beta_n = 2.5;
            Q10beta_m = 2.3;
            Q10beta_h = 2.3;
            Q10SERCA = 2.6;
            Q10PMCA = 2.35;
    
            %% Cl channel
            A = (1.0 / (1.0 + exp(((Voltage_SL - V_a) / A_a))));
            E_Cl =  - (log((Cl_o ./ Cl_i)) .* R .* T ./ F);
            I_Cl = ((g_Cl .* (Q10g_Cl.^ QCorr) .* (A ^ 4.0) .* (E_Cl - Voltage_SL )));
            J_Cl = ((I_Cl ./ (z_Cl .* F)) .* 1E09);
    
            %% KIR Channel
            E_K = (log((K_o ./ K_i)) .* R .* T ./ F);
            K_R = (K_o .* exp(( - delta .* E_K .* F ./ (R .* T))));
            I_K_IR = ((G_K .* (Q10g_KIR .^ QCorr) .* ((K_R ^ 2.0) ./ ...
                (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* ...
                exp(((2.0 .* (1.0 - delta) .* Voltage_SL .* F) ./ (R .* T))))) .* ...
                (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K - Voltage_SL )));
            J_K_IR = ((I_K_IR ./ (z_K .* F)) .* 1.0E9);
    
            %% KDR Channel
            I_K_DR = ((g_K .* (Q10g_K^ QCorr) .* (n ^ 4.0) .* ...
                h_K .* (E_K - Voltage_SL )));
            J_K_DR =  - ((I_K_DR ./ (z_K .* F)) .* 1E09);
    
            tau_hK = (exp(( -(Voltage_SL + 40.0) ./ 25.75)));
            h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_SL - V_hkinf) ./ A_hkinf))));
            J_r5 = ((h_Kinf - h_K) ./ tau_hK);
    
            alpha_n = (alpha_n0 .* (Q10alpha_n .^ QCorr) .* (Voltage_SL - V_n) ...
                ./ (1.0 - exp(( - (Voltage_SL - V_n) ./ K_alphan))));
            beta_n = (beta_n0 .* (Q10beta_n .^ QCorr).* exp(( - (Voltage_SL - V_n) ./ K_betan)));
            J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
    
            %% Na-K Pump
            fPump_NKX = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_SL .* F ./ (R .* T)))) + ...
                ((0.04 ./ 7.0) .* (exp((Na_o ./ 67300.0)) - 1.0) .* ...
                exp(( - Voltage_SL .* F ./ (R .* T))))) ^  - 1.0);
            C_NKX = (fPump_NKX .* F .* Q10g_NaK .* J_NaK_NKX ./ ...
                (((1.0 + (K_mK_NKX ./ K_o)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0))) *...
                ATP / (kATP + ATP);
            I_NKX_K = 2 * C_NKX;
            I_NKX_N = -3 * C_NKX;
    
            J_NKX_N =  - ((I_NKX_N ./ (z_Na .* F)) .* 1E09);
            J_NKX_K = ((I_NKX_K ./ (z_K .* F)) .* 1E09);
    
            J_NKX_tot = ( J_NKX_N * (SA_SL ./ vol_myo) / 3 );
    
            %% Na channel
            E_Na = (log((Na_o ./ Na_i)) .* R .* T ./ F);
            I_Na = (((g_Na * Q10g_Na .* (m ^ 3.0) * h * S) * ...
                (1.0)) + g_leakNa ) * (E_Na - Voltage_SL); %
            J_Na = ((I_Na ./ (z_Na .* F)) .* 1E09);
    
    
            S_inf = (1.0 ./ (1.0 + exp(((Voltage_SL - V_Sinf) ./ A_Sinf))));
            tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_SL + V_tau) ./ 100.0) ^ 2.0))));
            J_r3 = ((S_inf - S) ./ tau_S);
    
            alpha_m = (alpha_m0 .* (Q10alpha_m .^ QCorr) .* (Voltage_SL - V_m) ./ ...
                (1.0 - exp(( - (Voltage_SL - V_m) ./ K_alpham))));
            beta_m = (beta_m0 .* (Q10beta_m.^ QCorr) .* exp(( - (Voltage_SL - V_m) ./ K_betam)));
            J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
    
            alpha_h = (alpha_h0 .* (Q10alpha_h.^ QCorr) .*exp(( - (Voltage_SL - V_h) ./ K_alphah)));
            beta_h = (Q10beta_h .^ QCorr) .*(beta_h0 ./ (1.0 + exp(( - (Voltage_SL - V_h) ./ K_betah))));
            J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
    
            %% Na-Calcium Exchanger
            s1_NCX = (exp((nu_NCX .* Voltage_SL .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_o);
            Ka_NCX = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
            Qcorr_NCX = ((T - 310.0) ./ 10.0);
            s2_NCX = (exp(((nu_NCX - 1.0) .* Voltage_SL .* F ./ (R .* T))) .* (Na_o ^ 3.0) .* c_i);
    
            s3_NCX = ((Kmc_i_NCX .* (Na_o ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ...
                ((KmNa_EC_NCX ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + ...
                (Kmc_EC_NCX .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_o) + ((Na_o ^ 3.0) .* c_i));
    
            if c_o == 0 % set to zero in absence of EC calcium
                g_NCX = 0;
            end
            I_NCX_N = ( - (3.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX .* ...
                (s1_NCX - s2_NCX) ./ s3_NCX ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* ...
                Voltage_SL ./ (R .* T ./ F))))))));
            J_NCX_N = ((I_NCX_N ./ (z_Na .* F)) .* 1E09);
    
            Qcorr_NCX = ((T - 310.0) ./ 10.0);
            
            I_NCX_C = ((2.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX .* ...
                (s1_NCX - s2_NCX) ./ s3_NCX ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* ...
                Voltage_SL ./ (R .* T ./ F))))))));
            J_NCX_C =  - ((I_NCX_C ./ (z_Ca .* F)) .* 1E09);
    
            %% SOCE
            if ~isEst
                if any(expt == [1,3,5,7])
                    g_SOCE = p(95)*fluxFlagsCur(1);
                    SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
                    J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
                else
                    g_SOCE = 0; % Expt 2,4,6,8 have no SOCE.
                    SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
                    J_r6 = 0.0;
                end
                %Expt 1,2 are provided cont stimulus.
                if any(expt == [1,2])
                    continuousStim = true;
                else
                    continuousStim = false;
                end
            else
                g_SOCE = p(95)*fluxFlagsCur(1);
                continuousStim = false;
                SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
                J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
            end
            
            

            E_Ca = (log((c_o ./ c_i)) .* (R .* T) ./ (2.0 .* F));
            if k == 2 && g_SOCE > 0
                SOCEProb = SOCEProb_inf;
            end
            if c_o == 0 % set to zero in absence of EC calcium
                I_SOCE = 0;
            else
                I_SOCE = ((g_SOCE .* (E_Ca - Voltage_SL ) .* SOCEProb));
            end
            J_SOCE = ((I_SOCE ./ (z_Ca .* F)) .* 1E09);
    
            
    
            %% SERCA
            if ~isEst
                if any(expt == [7,8]) && t > 11*60
                    nu_SERCA = (nu_SERCA .* volFactor)*(0.1 + 0.9*exp(-(t-11*60)/10));
                else
                    nu_SERCA = (nu_SERCA .* volFactor);
                end
            else
                nu_SERCA = (nu_SERCA .* volFactor);
            end
            SERCA_ATPdep = (ATP / (kATP + ATP) );
            LumpedJ_SERCA = (Q10SERCA .^ QCorr) * (602.214179 * nu_SERCA * ...
                c_i ./ (K_SERCA + c_i)) * SERCA_ATPdep * SRMFrac;
    
            %% PMCA
            scalePMCA = 1/(1+exp(-(c_i-0.05)/.005));
            g_PMCA =scalePMCA*( (g_PMCA * vol_myo) / (700 * SA_SL ) );
            I_PMCA = (Q10PMCA .^ QCorr) * - (g_PMCA .* c_i ./ (K_PMCA + c_i)) * (ATP / (kATP + ATP));
            J_PMCA =  - ((I_PMCA ./ (z_Ca .* F)) .* 1E09) ;
    
            %% SL Calcium leak
            if c_o == 0
                I_CaLeak_SL = 0;
            else
                I_CaLeak_SL = ((g_SLLeak .* (E_Ca - Voltage_SL )));
            end
            J_CaLeak_SL = ((I_CaLeak_SL ./ (z_Ca .* F)) .* 1.0E09);
    
            %% SR Calcium Leak
            nu_leakSR = nu_leakSR * volFactor ;
            J_CaLeak_SR = 602.214179 * nu_leakSR * (c_SR - c_i) * SRMFrac;
    
            %% DHPR - RyR
            w_DHPR = w_RyR;
            f_DHPR = f_RyR;
            K_DHPR = K_RyR;
            L_DHPR = L_RyR;
            VBar_DHPR = VBar_RyR;
            openDHPR = ((((1.0 + (exp(((Voltage_SL - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* ...
                (f_DHPR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_SL - VBar_DHPR) ./ ...
                (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ...
                ((1.0 + exp(((Voltage_SL - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0)))) .* w_DHPR);
            
            if c_o == 0
                I_DHPR = 0;
            else
                I_DHPR = ((g0_DHPR .* openDHPR .* (E_Ca - Voltage_SL )));
            end
            J_DHPR = ((I_DHPR ./ (z_Ca .* F)) .* 1E09);

            voltProb = (((1.0 + (exp(((Voltage_SL - VBar_RyR) ./ (4.0 .* K_RyR))) .* ...
                (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_SL - VBar_RyR) ./ ...
                (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ...
                ((1.0 + exp(((Voltage_SL - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
            openProb = voltProb * w_RyR;
            LumpedJ_RyR = 602.214179 * j0_RyR * openProb * (c_SR - c_i);
    
            tau_w_r0 = (1.0 / (alpha_w_r0 * ( 1 + (c_i / K_w_r0)))) ; %(100.0 .* (1.0 + c_i)));
            wInf_r0 = (1.0 ./ (1.0 + (c_i / K_w_r0)));
            J_r0 = ((wInf_r0 - w_RyR) ./ tau_w_r0);
            
            KFlux_SL_myo = SAvolsCur(2);%(SA_SL ./ vol_myo);
            Ctot = (C_SL .* SA_SL) .* (Q10C_SL .^ QCorr);
    
            %% Calcium buffering in the myoplasm and SR
            %ATP Addition
            Parv = Parv_itot - CaParv - MgParv;
            g0_prime = g0 / 700;
    
            %Crossbridge Cycling (Calcium and Troponin Binding Process)
            %Calcium system
            if k == 1
                Trop = 0;
                D_1 = 0;
                D_0 = 0;
                k_onTrop2_D = 0;
                k_offTrop2_D = 0;
                k_onTrop1_D = 0;
                k_offTrop1_D = 0;
            else
                k_onTrop2_D = 0;
                k_offTrop2_D = 0;
                k_onTrop1_D = 0;
                k_offTrop1_D = 0;
                D_1 = 0;
                D_0 = 0;
                Trop = Trop_tot - CaTrop - CaCaTrop - D_2 - Pre_Pow - Post_Pow - D_1 - D_0;
            end
            dCT = k_onTrop1*c_i*Trop - k_offTrop1*CaTrop - k_onTrop2*c_i*CaTrop + k_offTrop2*CaCaTrop + 0.15*1000*D_1; % Calcium buffering with Troponin
            dCA = k_onATP*c_i*ATP - k_offATP*CATP; % Calcium buffering with ATP
            dCP = k_onParvCa*c_i*Parv - k_offParvCa*CaParv; % Calcium buffering with Parvalbumin
            dMP = k_onParvMg*Mg * Parv - k_offParvMg*MgParv; % Mg buffering with Parvalbumin
            dCCT = k_onTrop2*c_i*CaTrop - k_offTrop2*CaCaTrop - (k_onCa*CaCaTrop*ATP/700) + k_offCa*D_2;
            %Crossbridge attach/detachment
            dD2 = (k_onCa*CaCaTrop*ATP/700) - k_offCa*D_2 - f0*D_2 + fP*Pre_Pow + (g0_prime*Post_Pow*ATP) + k_onTrop2_D*D_1;

            %Concentration of Pre/Post Power Stroke Filaments
            hP_prime = hP*p_i_Myo / 3000;
            % hP_prime = hP;

            dPre = f0*D_2 - fP*Pre_Pow + hP_prime*Post_Pow - h0*Pre_Pow;
            dPost = -hP_prime*Post_Pow + h0*Pre_Pow - (g0_prime*Post_Pow*ATP);
            %ATP
            J_SERCA_tot = LumpedJ_SERCA * ( KMOLE / vol_myo );
            J_PMCA_tot = J_PMCA * KFlux_SL_myo;
            J_NKX_tot2 = J_NKX_tot;
    
            Jhyd = J_NKX_tot2 + (J_SERCA_tot/2) + J_PMCA_tot +  kHYD*(ATP / (kH + ATP) ); %(D_2*f0)/kHYD  (D_2*f0*p_i_SR)/1000
            dMA = k_onMA*Mg*ATP - k_offMA*MgATP; 
            dATP = -Jhyd - (k_onATP*c_i*ATP - k_offATP*CATP) - (k_onMA*Mg*ATP - k_offMA*MgATP) -(g0_prime*Post_Pow*ATP)- (k_onCa*CaCaTrop*ATP/700) + k_offCa*D_2;
            dADP = Jhyd + g0_prime*Post_Pow*ATP;
            Jhydtot = Jhyd + f0*D_2 - fP*Pre_Pow;
            
            %% Input stimulus for different conditions at frequency freq - square pulses of width 1 ms
            I_SL = 0;
            tMax = inf;
            if isEst
                if expt == 2
                    tMax = 0.05;
                elseif expt == 5
                    tMax = 0.07;
                elseif expt == 6
                    tMax = 0.07;
                elseif expt == 11
                    tMax = 0.33;
                else
                    tMax = 0.001;
                end
            end

            if (freq > 0 && t > 0 && t < tMax)
                if continuousStim || isEst
                    if (mod(t,1/freq) < pulsewidth)
                        I_SL = - ClampCurrent;
                    end
                elseif any(expt == [3,4]) % Prescribed stimulus - cycle set by exerciseParam
                    % if freq >= 100 && t > 20 && t < 40
                    %     I_SL = 0; %rest
                    % elseif freq < 100 && t > 60 && t < 90
                    %     I_SL = 0; %rest
                    % else
                    period = exerciseParam(1);
                    numPeriods = floor(t / period);
                    timeInCurrentPeriod = t - numPeriods * period;
                    if timeInCurrentPeriod <= exerciseParam(2)
                        if (mod(timeInCurrentPeriod,1/freq) < pulsewidth)
                            I_SL = - ClampCurrent;
                        end
                    end
                    % end
                elseif any(expt == [5,6]) % soce tests with eq period first 10 s
                    if (mod(t,1/freq) < pulsewidth) && t >= 10
                        I_SL = - ClampCurrent;
                    end
                end
            end

            %SR Phosphate
            PC_tot = p_i_SR * c_SR;
            kP = kP*volFactor;
            
            % PP = 1000*100;
            if freq == 0
                PiCaTerm = Bp*PiCa_SR*(PP-PC_tot);
            else
                if PC_tot > PP
                    Aterm = - Ap* (PC_tot)*(PC_tot -PP);
                    Bterm = 0;
                else
                    Aterm = 0;
                    Bterm = Bp*PiCa_SR*(PP-PC_tot);
                end
                % Hfunc = 0.5*(1+tanh((PC_tot-PP)/(0.01*PP)));
                PiCaTerm = Aterm + Bterm;
            end
            dPi_SR = kP*(p_i_Myo - p_i_SR)*SRMFrac/vol_SR + PiCaTerm;
            dPiCa = -PiCaTerm;
    
            %Myoplasmic Phosphate
            dPi_Myo = Jhyd + (h0*Pre_Pow - hP*Post_Pow*(p_i_Myo / 3000)) - bP*p_i_Myo - kP* (p_i_Myo - p_i_SR)*SRMFrac/vol_myo;
    
            % Buffering with CSQ
            k_onCSQ = k_offCSQ/K_CSQ;
            dCSQ = k_offCSQ*(B_CSQ - CSQ) - k_onCSQ * CSQ * c_SR;

            I_ionic = I_CaLeak_SL + I_Cl + I_DHPR + I_K_DR + I_K_IR + I_NCX_C...
                + I_NCX_N + I_NKX_K + I_NKX_N + I_Na + I_PMCA + I_SOCE;
            J_SL_Ca = J_SOCE + J_CaLeak_SL - J_NCX_C + J_DHPR - J_PMCA;
            J_SR_Ca = LumpedJ_RyR - LumpedJ_SERCA + J_CaLeak_SR;
            CaBuffer = - dCP - dCA - (k_onTrop1*c_i*Trop - k_offTrop1*CaTrop + k_onTrop2*c_i*CaTrop - k_offTrop2*CaCaTrop) -...
                (k_onTrop1_D*c_i*D_0 - k_offTrop1_D*D_1 + k_onTrop2_D*c_i*D_1 - k_offTrop2_D*D_2);
            
            % dPiCa = dPiCa*phosphateAccum;
            NaStim = false; % account for sodium fluxes due to applied current
            if NaStim
                J_NaStim = 1e9*I_SL/F; %#ok<UNRCH>
                if k == 2 % bulk
                    NaBulkIdx = sum(juncLocLogic(:)) + sum(bulkLocLogic(1:6));
                    J_NaStim = J_NaStim + (yinit(NaBulkIdx) - Na_i) * D_ionMyo*(SA_JBmyo/vol_myo) / diffusion_length;
                end
            else
                J_NaStim = 0;
            end
            % if t > 2
            %     I_SL = 0;
            % end
            
            %% Now define mito portion
            if length(geomParam) == 12
                mitoPhosScale = geomParam(11);
                mitoCaScale = geomParam(12);
            else
                mitoPhosScale = 1;
                mitoCaScale = 1;
            end
            % define mito params
            % pMito = ones(42,1);
            % PP_M = 1e5*pMito(1);
            p0Mito = [500, 0.01, 0.0006, 6, 0.38, 50, 0.1, 0.35,...
                      0.016, 0.008, 0.05, 450, 0.1, 1, 600, 100,...
                      177, 5, 25, 0.14, 0.1, 35000, 190, 8.5,...
                      10000, 5000, 0.111, 0.139, 0.5, 2, -30, 5000,...
                      1.8, 20, 3.43, 1000, 1.0, 15000, 250, 10000];
            Kdb = 500*pMito(1);%0.1*pMito(1);
            b_M = 0.01*pMito(2);%(1 + Kdb*Pi_M/(Kdb+c_M)^2 + 50)^(-1);
            N_A = 6.022e11; % molecules per pmol
            V_MCU = 0.0006*pMito(3); % uM/s 
            K_trans = 6*pMito(4); %uM CORRECTED
            K_act = 0.38*pMito(5); %uM CORRECTED
            n_au = 2.3; %uniporter cooperativity factor
            L = 50*pMito(6);
            p1 = 0.1*pMito(7); %1/mV
            V_NCX_M = 0.35*pMito(8); % uM/s ** DOUBLE CHECK THIS
            p2 = .016*pMito(9); %1/mV ** DOUBLE CHECK THIS
            k_mPTP = 0.008*pMito(10); % 1/s
            p3 = 0.05*pMito(11); % 1/mV
            % scaleGLY = 450*pMito(12);
            % scaleGLY = (1 + 20/(ATPset + ATP)^2)*450*pMito(12);
            scaleGLY = 450*pMito(12)* (ATPset-ATP)/10;
            if scaleGLY < 0
                scaleGLY = 0;
            end
            % scale glycolysis according to mito-myo volume ratio
            mitoMyoRef = 0.05/(1-0.05);
            mitoMyoCur = vol_mito/vol_myo;
            k_GLY = scaleGLY*(mitoMyoRef/mitoMyoCur)*(p_i_Myo/1000); % uM/s
            K_MCa = 0.1*pMito(13); % uM
            K_mNAD = 1*pMito(14);
            V_O = 600*pMito(15); % uM/s
            K_O = 100*pMito(16); % uM IS THIS 100 OR 1000??
            p4 = 177*pMito(17); % mV
            p5 = 5*pMito(18); % mV
            V_AGC = 25*pMito(19); % uM/s 
            K_AGC = 0.14*pMito(20); % uM
            K_iCa = 0.1*pMito(21); % uM OR MAYBE 0.3
            V_F1FO = 35000*pMito(22); % uM/s
            p6 = 190*pMito(23); %mV
            p7 = 8.5*pMito(24); %mV
            K_iATP = 10000*pMito(25); %uM
            V_ANT = 5000*pMito(26); %uM/s
            alpha_c = 0.111*pMito(27);
            alpha_m = 0.139*pMito(28);
            f_ANT = 0.5*pMito(29);
            p8 = 2*pMito(30); %uM/s per mV
            p9 = -30*pMito(31); %uM/s
            V_Pi = mitoPhosScale*5000*pMito(32);%500.0*pMito(32); % UNKNOWN!!
            % C_p = 1e-5; % pC/mV per sq micron
            % SA_mito = vol_mito/V_SA_mito;
            C_pEff = 1.8*pMito(33); % uM/mV effective capacitance CHECK!!!
            a1 = 20*pMito(34);
            a2 = 3.43*pMito(35);
            % g0_anion = 0;%1000*pMito(36);
            % V_anion = 230*pMito(37);
            KPi_i = 1000*pMito(36);%10*pMito(36);
            % KPi_M = 1000*p(39);
            % g_anion = 0;%g0_anion*Pi_M/1000;%y(7);
            % ginf_anion = Pi_M/1000;%(Pi_M + KPi_M);
            % tau_anion = 100.0;%*p(42);
            koffB = 1.0*pMito(37);
            konB = koffB / Kdb;
            % Ap_M = 0.1*1e-6 * pMito(39);
            % Bp_M = 0.1*1e-7 * pMito(40);
            ANPtot = 15000*pMito(38); %uM
            % ADP_M = ANPtot - ATP_M;
            % if ADP_M < 0
            %     error('ADP cannot be negative')
            %     % ADP_M = 0;
            % end
            NADtot = 250*pMito(39); %uM
            if length(pMito) == 39
                pMito(40) = 1.0;
            end
            KPi_M = 10000*pMito(40);
            % NAD = NADtot - NADH; % NAD+ in mito is implicit due to mass cons
            % if NAD < 0
            %     error('NAD cannot be negative')
            %     % NAD = 0;
            % end
            if abs(NAD + NADH - NADtot) > 0.01*NADtot
                % warning("NAD and NADH are not staying consistent")
                if NADH > NADtot
                    NADH = NADtot;
                end
                NAD = NADtot - NADH;
            elseif abs(ATP_M + ADP_M - ANPtot) > 0.01*ANPtot
                % warning("ATP and ADP are not staying consistent in mito")
                if ATP_M > ANPtot
                    ATP_M = ANPtot;
                end
                ADP_M = ANPtot - ATP_M;
            end
            if ADP_M < -1e-3
                error("ADP cannot be negative!!!!")
                % ADP_M = 0;
                % ATP_M = ANPtot;
            elseif NAD < -1e-3
                error("NAD cannot be negative!!!!")
                % NAD = 0;
                % NADH = NADtot;
            end
            % now include mitochondria calcium and ATP and phosphate fluxes
            mitoSAFactor = SA_mito / SA_mito_ref;
            MCUnum = (c_i/K_trans)*(1+c_i/K_trans)^3;
            MCUdenom = (1+c_i/K_trans)^4 + L/(1+c_i/K_act)^n_au;
            J_MCU = mitoSAFactor*V_MCU*(MCUnum/MCUdenom)*exp(p1*V_M);
            J_NCX_M = mitoSAFactor*V_NCX_M*(c_M/c_i)*exp(p2*V_M);
            J_mPTP = mitoCaScale*mitoSAFactor*k_mPTP*(c_i-c_M)*exp(p3*V_M);
            J_PDH = mitoSAFactor*k_GLY*(c_M/(K_MCa + c_M))*(1/(K_mNAD + NADH/NAD));
            J_O = mitoSAFactor*V_O*(NADH/(K_O+NADH))/(1+exp((V_M-p4)/p5));
            J_AGC = mitoSAFactor*V_AGC*(c_i/(K_AGC+c_i))*(K_iCa/(K_iCa+c_M))*exp(V_M*0.01);
            J_F1FO = mitoSAFactor*V_F1FO*(1/(1+exp((p6-V_M)/p7)))*(K_iATP/(K_iATP+ATP_M))*(Pi_M/1000);%*(ADP_M/100);
            ATP_tot = ATP + CATP + MgATP;
            if ADP > 0
                R_c = alpha_c*ATP_tot/ADP;
            else
                R_c = alpha_c*1e6;
            end
            if ATP_M > 0
                R_m = ADP_M/(alpha_m*ATP_M);
            else
                R_m = 1e6/alpha_m;
            end
            ANT_num = (1-R_c*R_m*exp(-F*V_M/(R*T)));
            % ANT_denom = (1+(alpha_c*ATP_tot/ADP)*exp(-f_ANT*F*V_M/(R*T)))*(1+ADP_M/(ATP_M*alpha_m));
            ANT_denom = (1+R_c*exp(-f_ANT*F*V_M/(R*T)))*(1+R_m);
            J_ANT = mitoSAFactor*V_ANT*ANT_num/ANT_denom;
            % if J_ANT < 0
            %     fprintf('catch it')
            % end
            J_Hleak = mitoSAFactor*(p8*V_M + p9);
            a3 = 0;%.5; % negative charge transport due to PiP
            J_Pi_M = mitoSAFactor*V_Pi*(p_i_Myo/KPi_i - Pi_M/KPi_M)/(1+ p_i_Myo/KPi_i + Pi_M/KPi_M);%*(Pi_i - Pi_M);%*exp(V_M*0.01);%(V_M-190)/100;
            % define other charge transport variables
            J_anion = 0;%-0*mitoSAFactor*g_anion*(V_M - V_anion);%(Pi_M/(Pi_M+KPi_M));%H_Pi; % effective variable for anion transport
            J_MCU_charge = J_MCU;
            J_F1FO_charge = J_F1FO;
            % calcium phosphate buffering in mito
            % PC_tot_M = Pi_M * c_M;
            % if PC_tot_M >= PP_M
            %     dPrecip = Ap_M* (PC_tot_M)*(PC_tot_M - PP_M);
            % else
            %     dPrecip = -Bp_M*PiCa_M*(PP_M - PC_tot_M) ;
            % end
            dPrecip = Pi_M * c_M * konB - koffB * PiCa_M;
            C_pEff = mitoSAFactor*C_pEff; % C_pEff correction
            mitoCouple = true;
            J_mito_Ca = (J_MCU - J_NCX_M + J_mPTP)*vol_mito; % in uM*um^3/s
            I_mito_Ca = (-J_NCX_M - 2*J_MCU_charge - 2*J_mPTP);
            if mitoCouple
                NADH_ATP_conv = 1;%2/10;
                Jprod = J_ANT*vol_mito/vol_myo + NADH_ATP_conv*k_GLY*vol_mito/vol_myo;
                if ATP_tot < 10
                    fprintf('troubleshoot')
                end
                dADP = dADP - Jprod;
                dPi_Myo = dPi_Myo - NADH_ATP_conv*k_GLY*vol_mito/vol_myo - J_Pi_M*vol_mito/vol_myo;
                J_to_I_mito = 602.214179*vol_mito*F/N_A; % convert flux in uM/s to total current 
                I_mito = (a1*J_O - a2*J_F1FO_charge + a3*J_Pi_M - J_Hleak + I_mito_Ca - J_ANT - J_AGC + J_anion)*J_to_I_mito;
                % J_NCX_N = J_NCX_N - 3*J_NCX_M*vol_mito/vol_myo;
            else
                I_mito = 0;
                Jprod =  100 * (ATPset - ATP) ;  %ATP production rate  (Post_Pow*g0*ATP)/1000 +
            end
            dATP = dATP + Jprod;
            dydtMito = zeros(6,1);
            dydtMito(1) = b_M*(J_mito_Ca/vol_mito - dPrecip);%konPrecip*Pi_M*c_M + konPrecip*Kdb*PiCa_M); % mito ca
            dydtMito(2) = J_F1FO - J_ANT; % mito ATP
            dydtMito(3) = -J_F1FO + J_Pi_M - dPrecip;%konPrecip*Pi_M*c_M + konPrecip*Kdb*PiCa_M; % mito Pi
            dydtMito(4) = J_PDH - J_O + J_AGC; % NADH - AGC is indirect?
            % dydt(5) = J_to_I*(a1*J_O - a2*J_F1FO - J_Hleak - J_NCX_M - 2*J_MCU - 2*J_mPTP - J_ANT - J_AGC)/(C_p*SA_mito);
            dydtMito(5) = (a1*J_O - a2*J_F1FO_charge + a3*J_Pi_M - J_Hleak + I_mito_Ca - J_ANT - J_AGC + J_anion)/(C_pEff);
            dydtMito(6) = dPrecip;%konPrecip*Pi_M*c_M - konPrecip*Kdb*PiCa_M;
            dydtMito(7) = -dydtMito(2); % for ADP_mito
            dydtMito(8) = -dydtMito(4); % for NAD_mito
            % dydt(7) = (ginf_anion - y(7))/tau_anion;
            % fluxesMito = [J_MCU,J_NCX_M,J_mPTP,J_PDH,J_O,J_AGC,J_F1FO,J_ANT,J_Hleak,J_Pi];

            dydtCur = [
                J_r6;    % rate for SOCEProb(1)
                (KMOLE * (LumpedJ_SERCA - LumpedJ_RyR - J_CaLeak_SR))/vol_SR - dPiCa + dCSQ; %c_SR (2)
                J_r5;    % rate for h_K (3)
                J_r0;    % rate for w_RyR (4)
                (1000/Ctot) * (SA_SL * I_ionic + I_SL);    % rate for Voltage_SL (5)
                KFlux_SL_myo * (J_Na - J_NKX_N + J_NCX_N + J_NaStim);    % rate for Na_i (6)
                (J_Cl .* KFlux_SL_myo);    % rate for Cl_i (7)
                (KFlux_SL_myo * J_SL_Ca) + (J_SR_Ca * KMOLE / vol_myo) - (mitoCouple*J_mito_Ca / vol_myo) + CaBuffer; % rate for c_i (8)
                J_r4;    % rate for n (9)
                J_r2;    % rate for m (10)
                J_r1;    % rate for h (11)
                J_r3;    % rate for S (12)
                KFlux_SL_myo * ( J_K_IR - J_K_DR + J_NKX_K);    % rate for K_i  (13)
                dCP;      % Rate for Parvalbumin bound Ca2+ (14)
                dMP;      % Rate for Parvalbumin bound Mg2+ (15)
                dCA;      % Rate for ATP bound Ca2+ (16)
                dCT;      % Rate for Trop bound Ca2+ (17)
                dCCT;     % Rate for Ca2+ bound TropCa2+ (18)
                dD2;      % Rate for Tropomyo opening from CaCaT bound (19)
                dPre;     % Rate for Pre-Power Stroke from D_2 bound (20)
                dPost;    % Rate for Post-Power Stroke from A_1 bound (21)
                dMA;      % Rate for ATP bound Mg2+ (22)
                dATP;     % Rate for free ATP (23)
                dPi_SR;   % Rate for SR Phosphate (24)
                dPiCa;    % Rate for Cal-Phos Precipitate (25)
                dPi_Myo;  % Rate for Myoplasmic Phosphate (26)
                -SAvols(1) * J_SL_Ca; % c_o (27)
                -SAvols(1) * (J_Na - J_NKX_N + J_NCX_N + J_NaStim); % Na_o (28)
                -SAvols(1) * (J_K_IR - J_K_DR + J_NKX_K); % K_o (29)
                -SAvols(1) *  J_Cl; % Cl_o (30)
                dCSQ; % CSQ (31)
                dADP; % ADP (32)
                % Then mito vars: Ca_Mito (33), ATP_Mito (34), Pi_Mito (35), NADH_Mito (36), Voltage_Mito (37), CaPi_Mito (38)
                ];
            dydtCur = [dydtCur; dydtMito];
            % quick mass conservation check before diffusion - would need
            % volume weighting
            % dATPtot = sum(dydtCur([16,22,23,34,19]));
            % dADPtot = sum(dydtCur([20,21,32])) - dydtCur(34);
            % dPitot = sum(dydtCur([26,35,20,24,25,38]));
            

            % now add diffusive fluxes
            otherLocs = find(1:size(yAll,2)~=k);
            for j = otherLocs
                if k == 1
                    tau_vec = D_vec(:) .* juncRatios(:) / diffusion_length;
                elseif k == 2
                    tau_vec = D_vec(:) .* bulkRatios(:) / diffusion_length;
                end
                dydtCur = dydtCur + (yAll(:,j)-yAll(:,k)).*tau_vec(:);
            end
            if ~phosphateAccum && freq==0
                % dydtCur(24:26) = 0;
                dydtCur([23:26,32]) = 0;
                % dydtCur([35,38]) = 0;
                dydtCur([34:35,39]) = 0;
            elseif ~phosphateAccum
                dydtCur(26) = 0;
            end
            if expt == 10
                crossbridgeIdx = 17:21;
                crossbridgeLogic = false(size(dydtCur));
                crossbridgeLogic(crossbridgeIdx) = true;
                dydtCur(~crossbridgeLogic) = 0;
            end
            RmagCur = abs(dydtCur ./ Nf);
            % if c_i < 0
            %     error('oops')
            % end
            % if any(isnan(dydtCur)) || any(imag(dydtCur) ~= 0)
            %     error('what???')
            % end

            dydt(curStartIdx:curStartIdx+sum(varLogicCur)-1) = dydtCur(varLogicCur);
            Rmag(curStartIdx:curStartIdx+sum(varLogicCur)-1) = RmagCur(varLogicCur);
            curStartIdx = curStartIdx + sum(varLogicCur);

            fluxes(k,:) = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
                LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR, Jhydtot, dydtCur(8), J_MCU, J_NCX_M, J_mPTP];
            % convert all fluxes to µM/s
            fluxes(k,1:5) = fluxes(k,1:5) * KFlux_SL_myo;
            fluxes(k,6:8) = fluxes(k,6:8) * KMOLE / vol_myo;
            fluxes(k,11:13) = fluxes(k,11:13) * vol_mito / vol_myo;
            currents(k,:) = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL, I_mito];
            currents(k,1:end-2) = SA_SL * currents(k,1:end-2); % convert to pA

        end
        % if any(~isreal(y)) || any(~isreal(dydt))
        %     fprintf("But why?\r\n")
        % end`
        if dydt(sum(juncLocLogic(1:40))) ~= -dydt(sum(juncLocLogic(1:36)))
            fprintf("Stop right there")
        end
        if all(Rmag < 0.00001) && freq==0
            dydt = zeros(length(y),1);
            fluxes = zeros(2,13);
            currents = zeros(2,14);
            return
        end
    end
end