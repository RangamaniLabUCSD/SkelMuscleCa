function [Time,Y,currtime,fluxes,currents,ySSFinal] = SkelMuscleCa_dydt(tSpan,freq, yinit, p,StartTimer,expt,phosphateAccum)
% input:
%     - tSpan is a vector of start and stop times
%       freq is a vector of test frequencies
%     - yinit is a vector of initial conditions for the state variables
%     - p is a vector of selected parameters to test
%     - StartTimer starts counting run time
%     - expt is the experimental value used for calculation
%     - phosphateAccum is a logical variable determining if phosphate
%       accumulation is accounted for
%
% output:
%     - T is the vector of times
%     - Y is the vector of state variables
%     - currtime is the total runtime
%     - fluxes: calcium fluxes at each time point, each row consists:
%            [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR]
%            in units of µM/s (in terms of myoplasmic conc)
%     - currents: Ionic and total current at each time point, each row consists of:
%            [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N,
%             I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL] in units of pA
% -------------------------------------------------------------------------

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

juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
% yinit(sum(juncLocLogic(1:24))) = p(74);
% yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:24))) = p(74);
% yinit(sum(juncLocLogic(1:26))) = p(75);
% yinit(sum(juncLocLogic)+sum(bulkLocLogic(1:26))) = p(75);

% ode15s settings
% pull out non-negative variables (everything except voltage)
nonNegJ = [1:sum(juncLocLogic(1:4)),sum(juncLocLogic(1:6)):sum(juncLocLogic(:))];
nonNegB = [1:sum(bulkLocLogic(1:4)),sum(bulkLocLogic(1:6)):sum(bulkLocLogic(:))];
nonNegB = nonNegB + nonNegJ(end);
if freq == 0 && ~any(expt==[7,8]) % Steady state condition
    options = odeset('RelTol',1e-3,'NonNegative',[nonNegJ, nonNegB]);
elseif freq == 0 && any(expt==[7,8]) % SOCE test without extra stimulus
    options = odeset('RelTol',1e-3,'MaxStep',0.1,'NonNegative', [nonNegJ, nonNegB]);
else
    options = odeset('RelTol',1e-3,'MaxStep',.001,'NonNegative', [nonNegJ, nonNegB]);
end
if expt < 0 % then corresponds to a case of estimation
    isEst = true;
    expt = abs(expt);
    %Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C]
            %[K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_t MJ_V],};
    Ca_o_exp = [1000, 1000, 2000, 2000, 2000,...
                2500, 1800, 5000, 2000, 2000, 1300];                                      %mM
    Na_o_exp = [138100, 138100, 150000, 150000, 150000,...
                 143800, 118240, 140000, 151000, 151000, 147000];                           %mM
    K_o_exp = [3900, 3900, 2000, 2000, 2000,...
               5000, 5330, 4000, 5000, 5000, 4000];                               %mM
    Cl_o_exp = [143700, 143700, 158000, 158000, 158000,...
                124000, 126170, 157000, 146000, 146000, 127400];                           %mM
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
    T = 295.15;
    Ca_EC = 1300;
    Na_EC = 147000.0;
    K_EC = 4000.0;
    Cl_EC = 128000.0;
end

% call ode15s
[Time,Y] = ode15s(@f,tSpan,yinit,options,p,freq); %pass extra arguments at the end

% return fluxes and currents at select times
StartTimer = tic;
fluxes = zeros(length(Time), 2*8);
currents = zeros(length(Time), 2*13);
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
if ssEst
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
        if currtime > 480
            error('too long to compute!')
        end

        %% State Variables - junctional and bulk
        Ca_EC_cur = Ca_EC;
        if ~isEst
            if any(expt == [7,8]) && (t > 6*60 && t < 17*60)
                Ca_EC_cur = 0.0;
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

        %% Model Geometry
        vol_SA_ratio = 0.01 ;       %µm
        volFraction_myo = 0.95 ;
        volFraction_SR = 0.05 ;
        volFraction_TT = 0.003 ;
        pulsewidth = 0.001 ;        %s
        R_fiber = 20;               %µm
        L_fiber = 100;              %µm

        vol_Fiber = pi * (R_fiber ^ 2) * L_fiber;
        vol_myo = volFraction_myo * vol_Fiber;
        SA_SL = 2 * pi * R_fiber * L_fiber;
        vol_SR = volFraction_SR * vol_Fiber;
        SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
        L_TT = R_fiber;
        vol_SA_SR = 1/10;
        SA_SR = vol_SR / vol_SA_SR;
        SRJ_occupancy = 0.5;
        SA_SRJ = SA_TT * SRJ_occupancy;
        SA_SRB = SA_SR - SA_SRJ;
        diffusion_length = p(99);%0.05;

        vol_myoJ = SA_TT*diffusion_length;
        vol_myoB = vol_myo - vol_myoJ;
        vol_SRJ = SA_SRJ*diffusion_length;
        vol_SRB = vol_SR - vol_SRJ;
        vol_TT = volFraction_TT * vol_Fiber;

        % define SA/vol ratios for junctional vs. bulk
        % [PM/EC, PM/myo, SRM/myo, SRM/lumen]
        SAvols = [SA_TT/vol_TT, SA_TT/vol_myoJ, SA_SRJ/vol_myoJ, SA_SRJ/vol_SRJ;
                  0.0         , SA_SL/vol_myoB, SA_SRB/vol_myoB, SA_SRB/vol_SRB];
        vols = [vol_TT, vol_myoJ, vol_SRJ;
                inf,    vol_myoB, vol_SRB];
        SAs = [SA_TT, SA_SRJ;
               SA_SL, SA_SRB];
        
        %% define diffusive flux rates
        D_ion = p(100);
        D_ionSR = 1*p(100); % slower maybe?
        D_ATP = p(101);
        D_parv = p(102);
        D_CSQ = p(103);
        D_Pi = p(104);
        D_V = 1/p(105);
        SA_JBmyo = (1-SRJ_occupancy)*SA_TT;
        SA_TTtop = vol_TT/L_TT;
        D_vec = [0, D_ionSR, 0, 0, D_V, D_ion, D_ion, D_ion, 0, 0,...
                   0, 0, D_ion, D_parv, D_parv, D_ATP, 0, 0, 0, 0,...
                   0, D_ATP, D_ATP, D_Pi, D_Pi, D_Pi, D_ion, D_ion, D_ion, D_ion, D_CSQ];
        juncRatios = [0, SA_SRJ/vol_SRJ, 0, 0, SA_SL/SA_TT, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, 0, 0,...
                      0, 0, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, 0, 0, 0, 0,...
                      0, SA_JBmyo/vol_myoJ, SA_JBmyo/vol_myoJ, SA_SRJ/vol_SRJ, SA_SRJ/vol_SRJ, SA_JBmyo/vol_myoJ,...
                      SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_TTtop/vol_TT, SA_SRJ/vol_SRJ];
        bulkRatios = [0, SA_SRJ/vol_SRB, 0, 0, 1.0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, 0, 0,...
                      0, 0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, 0, 0, 0, 0,...
                      0, SA_JBmyo/vol_myoB, SA_JBmyo/vol_myoB, SA_SRJ/vol_SRB, SA_SRJ/vol_SRB, SA_JBmyo/vol_myoB,...
                      0, 0, 0, 0, SA_SRJ/vol_SRB];

        %% define different fluxes/currents based on junctional vs. bulk
        % fluxes = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, 
        %           LumpedJ_SERCA, J_CaLeak_SR];
        % currents = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, 
        %             I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
        juncFluxFlags = [1, 1, 1, 1, 1, 1, 0.1, 1];
        bulkFluxFlags = [0, 1, 1, 0, 1, 0, 1, 1];
        juncCurrentFlags = [1, 0.1, 1, 0.45, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 0];
        bulkCurrentFlags = [1, 1. , 0, 1.  , 1, 1, 1, 1.0, 1.0 , 1. , 1, 0, 1];
        fluxFlags = {juncFluxFlags, bulkFluxFlags};
        currentFlags = {juncCurrentFlags, bulkCurrentFlags};
        varLogic = {juncLocLogic, bulkLocLogic};
        fluxes = zeros(2,8);
        currents = zeros(2,13);
        dydt = zeros(size(y));
        Rmag = zeros(size(y));
        curStartIdx = 1;
        Nf = [1;1000;1;1;100;1000;1000;0.1;1;1;1;1;100000;500;1000;...
            1;1;1;1;1;1;1;1;1;1;1;1300; 147000.0; 4000.0; 128000.0; 1000]; %Normalization factor
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
            bP = p(72);
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
            Kmc_EC_NCX_N = p(87);
            Kmc_EC_NCX_C = p(88);
            g_SLLeak = p(89)*fluxFlagsCur(2);
            L_RyR = p(90);
            g0_DHPR = p(91)*fluxFlagsCur(4);
            j0_RyR = p(92)*fluxFlagsCur(6);
            k_onTrop2 = p(93);
            k_offTrop2 = p(94);   
            B_CSQ = p(96);
            K_CSQ = p(97);
            k_offCSQ = p(98);
            ATPset = 700;%p(106);

    
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
    
            s3_NCX_N = ((Kmc_i_NCX .* (Na_o ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ...
                ((KmNa_EC_NCX ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + ...
                (Kmc_EC_NCX_N .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_o) + ((Na_o ^ 3.0) .* c_i));
    
            if c_o == 0 % set to zero in absence of EC calcium
                g_NCX = 0;
            end
            I_NCX_N = ( - (3.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX .* ...
                (s1_NCX - s2_NCX) ./ s3_NCX_N ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* ...
                Voltage_SL ./ (R .* T ./ F))))))));
            J_NCX_N = ((I_NCX_N ./ (z_Na .* F)) .* 1E09);
    
            s3_NCX_C = ((Kmc_i_NCX .* (Na_o ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ...
                ((KmNa_EC_NCX ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + ...
                (Kmc_EC_NCX_C .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_o) + ((Na_o ^ 3.0) .* c_i));
            Ka_NCX_C = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
            Qcorr_NCX = ((T - 310.0) ./ 10.0);
            
            I_NCX_C = ((2.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_C .* ...
                (s1_NCX - s2_NCX) ./ s3_NCX_C ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* ...
                Voltage_SL ./ (R .* T ./ F))))))));
            J_NCX_C =  - ((I_NCX_C ./ (z_Ca .* F)) .* 1E09);
    
            %% SOCE
            if ~isEst
                if any(expt == [1,3,5,7])
                    g_SOCE = p(95)*fluxFlagsCur(1);
                else
                    g_SOCE = 0; % Expt 2,4,6,8 have no SOCE.
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
            end
    
            E_Ca = (log((c_o ./ c_i)) .* (R .* T) ./ (2.0 .* F));
            if c_o == 0 % set to zero in absence of EC calcium
                I_SOCE = 0;
            else
                I_SOCE = ((g_SOCE .* (E_Ca - Voltage_SL ) .* SOCEProb));
            end
            J_SOCE = ((I_SOCE ./ (z_Ca .* F)) .* 1E09);
    
            SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
            J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
    
            %% SERCA
            volFactor = (vol_myo ./ (pi .* 0.26));
            if ~isEst
                if any(expt == [7,8]) && t > 11*60
                    nu_SERCA = (nu_SERCA .* volFactor)*(0.1 + 0.9*exp(-(t-11*60)/10));
                else
                    nu_SERCA = (nu_SERCA .* volFactor);
                end
            else
                nu_SERCA = (nu_SERCA .* volFactor);
            end
            LumpedJ_SERCA = (Q10SERCA .^ QCorr) * (602.214179 * nu_SERCA * ...
                c_i ./ (K_SERCA + c_i)) * (ATP / (kATP + ATP) );% * SRMFrac;
    
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
            else
                Trop = Trop_tot - CaTrop - CaCaTrop - D_2 - Pre_Pow - Post_Pow;
            end
            dCT = k_onTrop1*c_i*Trop - k_offTrop1*CaTrop - k_onTrop2*c_i*CaTrop + k_offTrop2*CaCaTrop; % Calcium buffering with Troponin
            dCA = k_onATP*c_i*ATP - k_offATP*CATP; % Calcium buffering with ATP
            dCP = k_onParvCa*c_i*Parv - k_offParvCa*CaParv; % Calcium buffering with Parvalbumin
            dMP = k_onParvMg*Mg * Parv - k_offParvMg*MgParv; % Mg buffering with Parvalbumin
            dCCT = k_onTrop2*c_i*CaTrop - k_offTrop2*CaCaTrop - (k_onCa*CaCaTrop*ATP/700) + k_offCa*D_2;
            %Crossbridge attach/detachment
            dD2 = (k_onCa*CaCaTrop*ATP/700) - k_offCa*D_2 - f0*D_2 + fP*Pre_Pow + (g0_prime*Post_Pow*ATP);

            %Concentration of Pre/Post Power Stroke Filaments
            dPre = f0*D_2 - fP*Pre_Pow + hP*Post_Pow*(p_i_Myo /3000) - h0*Pre_Pow;
            dPost = -hP*Post_Pow*(p_i_Myo / 3000) + h0*Pre_Pow - (g0_prime*Post_Pow*ATP);
            %ATP
            J_SERCA_tot = LumpedJ_SERCA * ( KMOLE / vol_myo );
            J_PMCA_tot = J_PMCA * KFlux_SL_myo;
            J_NKX_tot2 = J_NKX_tot;
    
            Jhyd = J_NKX_tot2 + (J_SERCA_tot/2) + J_PMCA_tot +  kHYD*(ATP / (kH + ATP) ); %(D_2*f0)/kHYD  (D_2*f0*p_i_SR)/1000
            Jprod =  kPROD * (ATPset - ATP) ;  %ATP production rate  (Post_Pow*g0*ATP)/1000 +
            dMA = k_onMA*Mg*ATP - k_offMA*MgATP; %did not include diffusion terms from Supplemental
            dATP = Jprod -Jhyd - (k_onATP*c_i*ATP - k_offATP*CATP) - (k_onMA*Mg*ATP - k_offMA*MgATP) -(g0_prime*Post_Pow*ATP)- (k_onCa*CaCaTrop*ATP/700) + k_offCa*D_2;
    
            %SR Phosphate
            PC_tot = p_i_SR * c_SR;
    
            if PC_tot >= PP && freq~=0
                dPi_SR = kP*(p_i_Myo - p_i_SR)  - Ap* (PC_tot)*(PC_tot -PP);
            else
                dPi_SR = kP*(p_i_Myo - p_i_SR)  + Bp* PiCa_SR* (PP - PC_tot) ;
            end
    
            %Calcium-Phophate Precipitate (SR)
            if PC_tot >= PP && freq~=0
                dPiCa = Ap * (PC_tot)*(PC_tot - PP);
            else
                dPiCa = -Bp * PiCa_SR* (PP - PC_tot) ;
            end
    
            %Myoplasmic Phosphate
            dPi_Myo = Jhyd + (h0*Pre_Pow - hP*Post_Pow*(p_i_Myo / 3000)) - bP*p_i_Myo - kP* (p_i_Myo - p_i_SR);
    
            % Buffering with CSQ
            k_onCSQ = k_offCSQ/K_CSQ;
            dCSQ = k_offCSQ*(B_CSQ - CSQ) - k_onCSQ * CSQ * c_SR; 
            % f_SR = 1/(1 + B_SRtot*K_SRBuffer./((K_SRBuffer+c_SR).^2));
    
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
                    ClampCurrent = ClampCurrent - 5000; % higher current for this expt
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
                elseif any(expt == [3,4]) % Prescribed stimulus - 500ms at 50Hz every 2.5s per Wei La-Pierre et al.
                    period = 2.5;
                    numPeriods = floor(t / period);
                    timeInCurrentPeriod = t - numPeriods * period;
                    if timeInCurrentPeriod <= 0.5
                        if (mod(timeInCurrentPeriod,1/freq) < pulsewidth)
                            I_SL = - ClampCurrent;
                        end
                    end
                elseif any(expt == [5,6]) % Prescribed stimulus - 3 x 60s with 120s rest in between
                    if (t > 0 && t <= 60) || (t > 180 && t <= 240) || (t > 360  && t <= 420)
                        if (mod(t,1/freq) < pulsewidth)
                            I_SL = - ClampCurrent;
                        end
                    end
                end
            end

            I_ionic = I_CaLeak_SL + I_Cl + I_DHPR + I_K_DR + I_K_IR + I_NCX_C...
                + I_NCX_N + I_NKX_K + I_NKX_N + I_Na + I_PMCA + I_SOCE;
            J_SL_Ca = J_SOCE + J_CaLeak_SL - J_NCX_C + J_DHPR - J_PMCA;
            J_SR_Ca = LumpedJ_RyR - LumpedJ_SERCA + J_CaLeak_SR;
            CaBuffer = - dCP - dCA - (k_onTrop1*c_i*Trop - k_offTrop1*CaTrop + k_onTrop2*c_i*CaTrop - k_offTrop2*CaCaTrop);
            
            dPiCa = dPiCa*phosphateAccum;

            J_NaStim = 0;%1e9*I_SL/F;
            dydtCur = [
                J_r6;    % rate for SOCEProb(1)
                (KMOLE * (LumpedJ_SERCA - LumpedJ_RyR - J_CaLeak_SR))/vol_SR - dPiCa + dCSQ; %c_SR (2)
                J_r5;    % rate for h_K (3)
                J_r0;    % rate for w_RyR (4)
                (1000/Ctot) * (SA_SL * I_ionic + I_SL);    % rate for Voltage_SL (5)
                KFlux_SL_myo * (J_Na - J_NKX_N + J_NCX_N + J_NaStim);    % rate for Na_i (6)
                (J_Cl .* KFlux_SL_myo);    % rate for Cl_i (7)
                (KFlux_SL_myo * J_SL_Ca) + (J_SR_Ca * KMOLE / vol_myo) + CaBuffer % rate for c_i (8)
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
                -SAvols(1) * J_SL_Ca; % c_o
                -SAvols(1) * (J_Na - J_NKX_N + J_NCX_N + J_NaStim); % Na_o
                -SAvols(1) * (J_K_IR - J_K_DR + J_NKX_K); % K_o
                -SAvols(1) *  J_Cl; % Cl_o
                dCSQ; % CSQ
                ];

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
            if ~phosphateAccum
                dydtCur(24:26) = 0;
            end
            if expt == 10
                crossbridgeIdx = 17:21;
                crossbridgeLogic = false(size(dydtCur));
                crossbridgeLogic(crossbridgeIdx) = true;
                dydtCur(~crossbridgeLogic) = 0;
            end
            RmagCur = abs(dydtCur ./ Nf);

            dydt(curStartIdx:curStartIdx+sum(varLogicCur)-1) = dydtCur(varLogicCur);
            Rmag(curStartIdx:curStartIdx+sum(varLogicCur)-1) = RmagCur(varLogicCur);
            curStartIdx = curStartIdx + sum(varLogicCur);

            fluxes(k,:) = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA,...
                LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];
            % convert all fluxes to µM/s
            fluxes(k,1:5) = fluxes(k,1:5) * KFlux_SL_myo;
            fluxes(k,6:8) = fluxes(k,6:8) * KMOLE / vol_myo;
            currents(k,:) = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C,... 
                I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
            currents(k,1:end-1) = SA_SL * currents(k,1:end-1);

        end
        if all(Rmag < 0.00001) && freq==0
            dydt = zeros(length(y),1);
            fluxes = zeros(2,8);
            currents = zeros(2,13);
            return
        end
    end
end