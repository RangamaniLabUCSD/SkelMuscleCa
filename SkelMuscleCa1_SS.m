function [Time,Y,currtime] = SkelMuscleCa1_SS(tSpan,freq, lowATP, yinit, p,StartTimer,expt)
% [T,Y,yinit,param] = SkelMuscleCa_AK(argTimeSpan,argYinit,argParam)
%
% input:
%     tSpan is a vector of start and stop times
%     freq is a vector of test frequencies
%     lowATP is true for low ATP conditions
%     p is a vector of selected parameters to test
%     yinit is a vector of initial conditions for the state variables
%     expt is the exapiremental value used for calculation
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yInf is the predicted SS values
% -------------------------------------------------------------------------

timeSpan = [0.0 1.0];

if nargin >= 1
    if ~isempty(tSpan)
        %
        % TimeSpan overridden by function arguments
        %
        timeSpan = tSpan;
    end
end
%
% invoke the integrator
%

if freq == 0
    options = odeset('RelTol',1e-8,'NonNegative',[1:4,6:16]);
else
    options = odeset('RelTol',1e-8,'MaxStep',.001,'NonNegative',[1:4,6:16]);%,'OutputFcn',@odeplot);
end
[Time,Y] = ode15s(@f,tSpan,yinit,options,p,freq,lowATP); %pass extra arguments at the end



% -------------------------------------------------------
% ode rate
    function dydt = f(t,y,p,freq,lowATP)
        % State Variables
        SOCEProb = y(1);
        c_SR = y(2);
        h_K = y(3);
        w_RyR = y(4);
        %w_DHPR = y(5);
        Voltage_PM = y(5);
        Na_i = y(6);
        Cl_i = y(7);
        c_i = y(8);
        n = y(9);
        m = y(10);
        h = y(11);
        S = y(12);
        K_i = y(13);
        CaParv = y(14);
        MgParv = y(15);
        CATP = y(16);

        %% Parameters
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
        C_PM = p(11);
        c_ref = p(12);
        ClampCurrent = p(13);
        delta = p(14);
        f_RyR = p(15);
        g_Cl = p(16);
        g_K = p(17);
        G_K = p(18);
        g_Na = p(19);
        g_NCX = p(20);
        J_NaK_NKX = p(21);
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

        %% Optimization Parameters

        % ClampCurrent = p(1) ;
        % K_S = p(2) ;
        % delta = p(3) ;
        % beta_m0 = p(4) ;
        % K_betam = p(5) ;
        % alpha_m0 = p(6) ;
        % K_alpham = p(7) ;
        % K_RyR = p(8) ;
        % f_RyR = p(9) ;
        % 
        % % Non-optimized parameters
        % A_a = 150;
        % A_hkinf = 7.5;
        % A_Sinf = 5.8;
        % alpha_h0 = 8.1;
        % alpha_n0 = 13.1;
        % alpha_w_r0 = 100;
        % beta_h0 = 4380;
        % beta_n0 = 67;
        % C_PM = 0.01;
        % c_ref = 500;
        % g_Cl = 0.1;
        % g_K = 0.648;
        % G_K = 0.111;
        % g_Na = 8.04;
        % g_NCX = 0.129;
        % J_NaK_NKX = 0.00003542;
        % K_alphah = 14.7;
        % K_alphan = 7;
        % K_betah = 9;
        % K_betan = 40;
        % K_K = 950000000;
        % K_mK_NKX = 1000;
        % K_mNa_NKX = 13000;
        % K_PMCA = 0.17;
        % K_SERCA = 1;
        % K_w_r0 = 1;
        % Kdact_NCX = 0.05;
        % Kmc_i_NCX = 6.59;
        % ksat_NCX = 0.32;
        % nu_NCX = 0.27;
        % S_i = 10000;
        % tau_SOCEProb = 0.01;

        %% Constants
        vol_SA_ratio = 0.01 ;
        volFraction_cyto = 0.95 ;
        volFraction_SR = 0.05 ;
        volFraction_TT = 0.003 ;
        pulsewidth = 0.001 ;
        R_fiber = 20;
        
        %% Experimental Inputs
        %Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C]
        %[K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_t MJ_V],};
        Ca_o_exp = [1000, 1000, 2000, 2000,2000, 2500, 1800, 2000,  5000];% 2000]; 1300,             %mM
        Na_o_exp = [138100, 138100, 150000, 150000, 150000, 143800, 118240, 140000, 151000.0, 151000];%133000] ;147000.0 %mM
        K_o_exp = [3900, 3900, 2000, 2000, 2000, 5000, 5330, 4000, 5900 , 5000]; %5900]; %4000               %mM

        %K_o_exp = [3900, 4000, 4000, 4000, 5000, 5330];               %mM
        Cl_o_exp = [143700, 143700, 158000, 158000,158000, 124000, 126170, 157000, 128000, 146000];  %mM  128000
        Temp = [(273+22),(273+22),(273+20),(273+22), (273+22),(273+ 26),(273+22),(273+22),(273+35),(273+22)]; %K
        
        if expt == 2 
            expt_n = 1; 
        elseif expt < 4 || expt > 5
            expt_n = 3;
        else
            expt_n = expt;
         
        end
        % Constants
        c_EC_init_uM = Ca_o_exp(expt_n); %1300.0;
        Cl_EC_init_uM = Cl_o_exp(expt_n); %128000.0;
        K_EC_init_uM = K_o_exp(expt_n); %4000.0;
        Na_EC_init_uM = Na_o_exp(expt_n); %147000.0;

        carrierValence_CaLeak_PM = 2.0;
        carrierValence_Cl = -1.0;
        carrierValence_DHPR = 2.0;
        carrierValence_K_DR = 1.0;
        carrierValence_K_IR = 1.0;
        carrierValence_Na = 1.0;
        carrierValence_NCX_C = 2.0;
        carrierValence_NCX_N = 1.0;
        carrierValence_NKX = 1.0;
        carrierValence_PMCA = 2.0;
        carrierValence_RyR = 1.0;
        carrierValence_SOCE = 2.0;

        h_KUnit = 1.0;
        hUnit = 1.0;
        L_fiber = 100;
        mUnit = 1.0;
        nUnit = 1.0;
        SOCEProbUnit = 1.0;
        SUnit = 1.0;
        tauUnit = 1.0;
        wUnit_DHPR = 1.0;
        wUnit_RyR = 1.0;

        F = 96485.3321;
        PI = 3.141592653589793;
        R = 8314.46261815;
        T = Temp(expt); %293;
        KMOLE = 0.001660538783162726;

        V_a = 70.0;
        V_h = -45.0;
        V_hkinf = -40.0;
        V_m = -46.0;
        V_n = -40.0;
        V_Sinf = -78.0;
        V_tau = 90.0;
        %VBar_DHPR = -20;
        VBar_RyR = -20.0;

        QCorr = (T-293)/10; %(273+20))/10;
        Q10NCX = 1.57;
        Q10g_K = 1.5;
        Q10g_Na = 1.5;
        Q10g_KIR= 1.55;
        Q10g_NaK = 1 ;
        Q10g_Cl = 1.5;
        Q10C_PM = 1.02;
        Q10alpha_n = 2.5;
        Q10alpha_m = 2.3;
        Q10alpha_h = 2.5;
        Q10beta_n = 2.5;
        Q10beta_m = 2.3;
        Q10beta_h = 2.3;
        Q10SERCA = 2.6;
        Q10PMCA = 2.35;


        % Geometry
        vol_Fiber = PI * (R_fiber^2) * L_fiber;
        vol_cyto = volFraction_cyto* vol_Fiber;
        SA_PM = 2* PI * R_fiber * L_fiber ;
        vol_SR = volFraction_SR * vol_Fiber ;
        SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio;
        TTFrac = SA_TT / SA_PM;

        % Functions

        K_EC = K_EC_init_uM;
        E_K_K_IR = (log((K_EC ./ K_i)) .* R .* T ./ F);
        K_R = (K_EC .* exp(( - delta .* E_K_K_IR .* F ./ (R .* T))));
        I_K_IR = ((G_K .* (Q10g_KIR .^ QCorr) .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_PM .* F) ./ (R .* T))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_PM)) .* (1.0 + TTFrac));
        J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* F)) .* 1.0E9);
        Na_EC = Na_EC_init_uM;
        E_Na = (log((Na_EC ./ Na_i)) .* R .* T ./ F);
        KmNa_i_NCX = (12.29 .* 1000.0);
        volFactor = (vol_cyto ./ (PI .* 0.26));
        nu_SERCA = (4875.0 .* volFactor);
        LumpedJ_SERCA = (Q10SERCA .^ QCorr) * (602.214179 * nu_SERCA * c_i ./ (K_SERCA + c_i));
        I_Na = ((g_Na * Q10g_Na .* ((m ./ mUnit) ^ 3.0) * (h ./ hUnit) * (S ./ SUnit) * (E_Na - Voltage_PM)) * (1.0 + (0.1 .* TTFrac)));
        J_Na = ((I_Na ./ (carrierValence_Na .* F)) .* 1E09);
        fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
        I_NKX_K = ((2.0 .* (fPump_NKX_K .* F .* Q10g_NaK .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
        g_PMCA = (5.37 ./ SA_PM);
        I_PMCA = (Q10PMCA .^ QCorr) * ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
        fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
        I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* F .* Q10g_NaK .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
        g_SOCE = (0.01 ./ 210.44);
        c_EC = c_EC_init_uM;
        E_Ca = (log((c_EC ./ c_i)) .* (R .* T) ./ (2.0 .* F));
        I_SOCE = ((g_SOCE .* (E_Ca - Voltage_PM) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
        E_K_K_DR = (log((K_EC ./ K_i)) .* R .* T ./ F);
        I_K_DR = ((g_K .* (Q10g_K^ QCorr) .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_PM)) .* (1.0 + (0.45 .* TTFrac)));
        g_PMLeak = 5.0 .* 2.0E-6;
        I_CaLeak_PM = ((g_PMLeak .* (E_Ca - Voltage_PM)) .* (1.0 + TTFrac));
        s1_NCX_N = (exp((nu_NCX .* Voltage_PM .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_EC);
        Ka_NCX_N = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
        Qcorr_NCX = ((T - 310.0) ./ 10.0);
        s2_NCX_N = (exp(((nu_NCX - 1.0) .* Voltage_PM .* F ./ (R .* T))) .* (Na_EC ^ 3.0) .* c_i);
        KmNa_EC_NCX_N = (87.5 .* 1000.0);
        Kmc_EC_NCX_N = (1.3 .* 1000.0);
        s3_NCX_N = ((Kmc_i_NCX .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ((KmNa_EC_NCX_N ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + (Kmc_EC_NCX_N .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
        I_NCX_N = ( - (3.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_N .* (s1_NCX_N - s2_NCX_N) ./ s3_NCX_N ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* Voltage_PM ./ (R .* T ./ F))))))) .* (1.0 + TTFrac));
        A = (1.0 / (1.0 + exp(((Voltage_PM - V_a) / A_a))));
        Cl_EC = Cl_EC_init_uM;
        E_Cl =  - (log((Cl_EC ./ Cl_i)) .* R .* T ./ F);
        I_Cl = ((g_Cl .* (Q10g_Cl.^ QCorr) .* (A ^ 4.0) .* (E_Cl - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
        s2_NCX_C = (exp(((nu_NCX - 1.0) .* Voltage_PM .* F ./ (R .* T))) .* (Na_EC ^ 3.0) .* c_i);
        Kmc_EC_NCX_C = (1.6 .* 1000.0);
        KmNa_EC_NCX_C = (87.5 .* 1000.0);
        s3_NCX_C = ((Kmc_i_NCX .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ((KmNa_EC_NCX_C ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + (Kmc_EC_NCX_C .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
        Ka_NCX_C = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
        Qcorr_NCX = ((T - 310.0) ./ 10.0);
        s1_NCX_C = (exp((nu_NCX .* Voltage_PM .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_EC);
        I_NCX_C = ((2.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_C .* (s1_NCX_C - s2_NCX_C) ./ s3_NCX_C ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* Voltage_PM ./ (R .* T ./ F))))))) .* (1.0 + TTFrac));
        %L_RyR = (1000.0 ./ 0.002);
        %tau_w_r7 = (1.0 ./ (100.0 .* (1.0 + c_i)));
        % alpha_w_r0 = 100;
        % K_w_r0 = 1;

        tau_w_r0 = (1.0 / (alpha_w_r0 * ( 1 + (c_i / K_w_r0)))) ;%(100.0 .* (1.0 + c_i)));

        alpha_n = (alpha_n0 .* (Q10alpha_n .^ QCorr) .* (Voltage_PM - V_n) ./ (1.0 - exp(( - (Voltage_PM - V_n) ./ K_alphan))));
        alpha_m = (alpha_m0 .* (Q10alpha_m .^ QCorr) .* (Voltage_PM - V_m) ./ (1.0 - exp(( - (Voltage_PM - V_m) ./ K_alpham))));
        
        L_RyR = (1000.0 ./ 0.002);
        w_DHPR = w_RyR;
        f_DHPR = f_RyR;
        K_DHPR = K_RyR;
        L_DHPR = L_RyR;
        VBar_DHPR = VBar_RyR;
        openDHPR = ((((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ((1.0 + exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0)))) .* (w_DHPR ./ wUnit_DHPR));
        g0_DHPR = (9.39 ./ 100.0);
        I_DHPR = ((g0_DHPR .* openDHPR .* (E_Ca - Voltage_PM)) .* TTFrac);
        voltProb = (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ((1.0 + exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
                
        j0_RyR = 300.0 * volFactor;
        openProb = voltProb * w_RyR;
        LumpedJ_RyR = 602.214179 * j0_RyR * openProb * (c_SR - c_i);

        nu_leakSR = 0.2*volFactor;  %1.1338; %
        J_CaLeak_SR = 602.214179 * nu_leakSR * (c_SR - c_i);

        alpha_h = (alpha_h0 .* (Q10alpha_h.^ QCorr) .*exp(( - (Voltage_PM - V_h) ./ K_alphah)));
        J_CaLeak_PM = ((I_CaLeak_PM ./ (carrierValence_CaLeak_PM .* F)) .* 1.0E09);
        %wInf_r7 = (1.0 ./ (1.0 + c_i));
        %J_r7 = ((wInf_r7 - w_DHPR) ./ tau_w_r7);
        SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
        J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
        tau_hK = (exp(( - (Voltage_PM + 40.0) ./ 25.75)) .* tauUnit);
        h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_hkinf) ./ A_hkinf))));
        J_r5 = ((h_Kinf - h_K) ./ tau_hK);
        beta_n = (beta_n0 .* (Q10beta_n .^ QCorr).* exp(( - (Voltage_PM - V_n) ./ K_betan)));
        J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
        S_inf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_Sinf) ./ A_Sinf))));
        tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_PM + V_tau) ./ 100.0) ^ 2.0))));
        J_r3 = ((S_inf - S) ./ tau_S);
        beta_m = (beta_m0 .* (Q10beta_m.^ QCorr) .* exp(( - (Voltage_PM - V_m) ./ K_betam)));
        J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
        beta_h = (Q10beta_h .^ QCorr) .*(beta_h0 ./ (1.0 + exp(( - (Voltage_PM - V_h) ./ K_betah))));
        J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
        wInf_r0 = (1.0 ./ (1.0 + (c_i / K_w_r0)));
        J_r0 = ((wInf_r0 - w_RyR) ./ tau_w_r0);
        J_NKX_N =  - ((I_NKX_N ./ (carrierValence_NKX .* F)) .* 1E09);
        J_SOCE = ((I_SOCE ./ (carrierValence_SOCE .* F)) .* 1E09);
        J_NKX_K = ((I_NKX_K ./ (carrierValence_NKX .* F)) .* 1E09);
        J_K_DR =  - ((I_K_DR ./ (carrierValence_K_DR .* F)) .* 1E09);

        %KFlux_SRM_SR = (SA_SRM ./ vol_SR);
        %KFlux_SRM_SR = 1/vol_SR;
        KFlux_PM_cyto = (SA_PM ./ vol_cyto);
        J_DHPR = ((I_DHPR ./ (carrierValence_DHPR .* F)) .* 1E09);
        J_NCX_N = ((I_NCX_N ./ (carrierValence_NCX_N .* F)) .* 1E09);
        J_NCX_C =  - ((I_NCX_C ./ (carrierValence_NCX_C .* F)) .* 1E09);
        %KFlux_SRM_cyto = (SA_SRM ./ vol_cyto);
        %KFlux_SRM_cyto = (1/ vol_cyto);
        %I_PM = - ClampCurrent;
        J_PMCA =  - ((I_PMCA ./ (carrierValence_PMCA .* F)) .* 1E09);
        device_PM.Capacitance = (C_PM .* SA_PM) .* (Q10C_PM .^ QCorr);
        J_Cl = ((I_Cl ./ (carrierValence_Cl .* F)) .* 1E09);


        %% Emmet - added code to VCell function
        % calcium buffering in the cytosol and SR

        k_onATP = 0.01364*1000;
        k_offATP = 30*1000;
        k_onParvCa = 41.7; %(uM s)^-1
        k_offParvCa = 0.5; %s^-1
        k_onParvMg = 0.033;
        k_offParvMg = 3;
        Parv_itot = 1500;
        ATP_itot = 8000;

        ATP = ATP_itot - CATP;
        Parv = Parv_itot - CaParv - MgParv;
        Mg = 1000; %1mM constant concentration


        % if freq > 0
        %      %LumpedJ_SERCA = 0.8 * LumpedJ_SERCA;
        %      %LumpedJ_RyR = 2 * LumpedJ_RyR;
        %     % J_DHPR = 2* J_DHPR;
        %     % Parv = 0.5 * Parv;
        %     %J_PMCA = 2 * J_PMCA;
        % 
        % end


        dCP = k_onParvCa*c_i*Parv - k_offParvCa*CaParv;
        dMP = k_onParvMg*Mg * Parv - k_offParvMg*MgParv;
        dCA = k_onATP*c_i*ATP - k_offATP*CATP;

        %f_i = 1/(1 + (K_iATP*ATP_itot/((K_iATP+c_i).^2)) + (K_iParv*Parv_itot/((K_iParv+c_i).^2))); % Assuming rapid buffering

        B_SRtot = 31000;
        K_SRBuffer = 500;
        f_SR = 1/(1 + B_SRtot*K_SRBuffer./((K_SRBuffer+c_SR).^2));

        currtime = toc(StartTimer);     

        % prescribed stimulus (applied current I_PM at frequency freq - square
        % pulses of width 1 ms)
        I_PM = 0;

        % n_pulses = [1, 5, 1, 1, 5, 1, 1, 1, 1, 1];
        % 
        % if expt == 8
        %     pulsewidth = pulsewidth/2; %0.0005;
        % end
        % for p = 1 : n_pulses(expt)
        %     if t >= 2*(p-1)*pulsewidth && t <= (2*p - 1) * pulsewidth
        %         if expt == 6
        %             I_PM =  - ClampCurrent + 5000;
        %         else
        %             I_PM = - ClampCurrent;
        %         end
        %     end
        % end
        
        if expt == 2 
            if t > 0 && t < 0.05
                if (mod(t,1/freq) < pulsewidth)% || mod(t,1/freq) > ((1/freq)-pulsewidth/2)
                    I_PM = - ClampCurrent; 
                end
            end
        elseif expt == 5
            if t > 0 && t < 0.07
                if (mod(t,1/freq) < pulsewidth)% || mod(t,1/freq) > ((1/freq)-pulsewidth/2)
                    I_PM = - ClampCurrent; 
                end
            end
        elseif expt == 6 
            if t > 0 && t < 0.07
                if (mod(t,1/freq) < pulsewidth)% || mod(t,1/freq) > ((1/freq)-pulsewidth/2)
                    I_PM =  - ClampCurrent + 5000; %  Increased I_PM from -20k to -25kPA
                end
            end
        elseif expt == 8
            if t > 0 && t < 0.001
                pulsewidth = 0.05;
                I_PM = - ClampCurrent ;
            end
        elseif expt < 2 || expt > 2 || expt < 5 || expt > 6 
            if t > 0 && t < 0.001
                I_PM = - ClampCurrent;
            end
        end


        % switch lowATP condition on to test predictions in energy deficit (all
        % pumps go to half activity)
        % if lowATP
        %     LumpedJ_SERCA = LumpedJ_SERCA*0.5;
        %     J_PMCA = J_PMCA*0.5;
        %     I_PMCA = I_PMCA*0.5;
        %     J_NKX_K = J_NKX_K*0.5;
        %     J_NKX_N = J_NKX_N*0.5;
        %     I_NKX_K = I_NKX_K*0.5;
        %     I_NKX_N = I_NKX_N*0.5;
        % end

        %% Rates

        Nf = [1;1000;1;1;100;1000;1000;0.1;1;1;1;1;100000;500;1000;1]; %Normalization factor
        if freq == 0 && currtime > 600
            dydt = zeros(16,1);
            return;
        end
        dydt = [
            J_r6;    % rate for SOCEProb
            (f_SR * KMOLE * (LumpedJ_SERCA - LumpedJ_RyR - J_CaLeak_SR))/vol_SR; %c_SR
            J_r5;    % rate for h_K
            J_r0;    % rate for w_RyR
            %J_r7;    % rate for w_DHPR
            (1000 /device_PM.Capacitance) * (SA_PM * (I_CaLeak_PM + I_Cl +...
            I_DHPR + I_K_DR + I_K_IR + I_NCX_C + I_NCX_N + I_NKX_K + I_NKX_N ...
            + I_Na + I_PMCA + I_SOCE) + I_PM);    % rate for Voltage_PM
            KFlux_PM_cyto * (J_Na - J_NKX_N + J_NCX_N);    % rate for Na_i
            (J_Cl .* KFlux_PM_cyto);    % rate for Cl_i
            %f_i*((KFlux_PM_cyto * (J_SOCE + J_CaLeak_PM - J_NCX_C + J_DHPR - J_PMCA))...
            %+ ((LumpedJ_RyR - LumpedJ_SERCA + J_CaLeak_SR) * KMOLE /vol_cyto)); % rate for c_i assuming rapid buffering
            (KFlux_PM_cyto * (J_SOCE + J_CaLeak_PM - J_NCX_C + J_DHPR - J_PMCA))...
            + ((LumpedJ_RyR - LumpedJ_SERCA + J_CaLeak_SR) * KMOLE / vol_cyto) - dCP - dCA; % rate for c_i
            J_r4;    % rate for n
            J_r2;    % rate for m
            J_r1;    % rate for h
            J_r3;    % rate for S
            KFlux_PM_cyto * ( J_K_IR - J_K_DR + J_NKX_K);    % rate for K_i
            dCP;     % Rate for Parvalbumin bound Ca2+
            dMP;     % Rate for Parvalbumin bound Mg2+
            dCA;     % Rate for ATP bound Ca2+
            ];

        R = abs(dydt ./ Nf); %Normalize
        if all(R < 0.00001)
           dydt = zeros(16,1);
           %fprintf('Threshold met \n')
           return
        end
      
    end
end

