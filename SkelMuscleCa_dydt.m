function [Time,Y,currtime,fluxes,currents] = SkelMuscleCa_dydt(tSpan,freq, lowATP, yinit, p,StartTimer,expt)
% [T,Y,yinit,param] = SkelMuscleCa_AK(argTimeSpan,argYinit,argParam)
%
% input:
%     tSpan is a vector of start and stop times
%     freq is a vector of test frequencies
%     lowATP is true for low ATP conditions
%     p is a vector of selected parameters to test
%     yinit is a vector of initial conditions for the state variables
%     expt is the experimental value used for calculation
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     currtime:
%     fluxes: calcium fluxes at each time point, each row has: 
%            [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR]
%            in units of uM/s (in terms of myoplasmic conc)
% -------------------------------------------------------------------------

if freq == 0
    options = odeset('RelTol',1e-3,'NonNegative',[1:4,6:17]);
else
    options = odeset('RelTol',1e-3,'MaxStep',.001,'NonNegative',[1:4,6:17]);%,'OutputFcn',@odeplot);
end
[Time,Y] = ode15s(@f,tSpan,yinit,options,p,freq,lowATP); %pass extra arguments at the end

fluxes = zeros(length(Time), 8);
currents = zeros(length(Time), 13);
for i = 1:length(Time)
    [~, fluxes(i,:), currents(i,:)] = f(Time(i), Y(i,:), p, freq, lowATP);
end

% -------------------------------------------------------
% ode rate
    function [dydt, fluxes, currents] = f(t,y,p,freq,lowATP)
        % State Variables
        SOCEProb = y(1);
        c_SR = y(2);
        h_K = y(3);
        w_RyR = y(4);
        %w_DHPR = y(5);
        Voltage_SL = y(5);
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
        CaTrop = y(17);

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
        C_SL = p(11);
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
        nu_SERCA = p(42);
        g_PMCA = p(43);
        nu_leakSR = p(44);
        g_leakNa = p(45);

        %% Constants
        vol_SA_ratio = 0.01 ;       
        volFraction_myo = 0.95 ;
        volFraction_SR = 0.05 ;
        volFraction_TT = 0.003 ;
        pulsewidth = 0.001 ;        %s
        R_fiber = 20;              %um
       
        % Constants
        c_EC_init_uM = 1300; %2000;        %uM
        Cl_EC_init_uM = 128000.0;
        K_EC_init_uM = 4000.0;
        Na_EC_init_uM = 147000.0; 
        % Cl_EC_init_uM = 158000;     %uM
        % K_EC_init_uM = 2000;        %uM
        % Na_EC_init_uM = 150000;     %uM

        carrierValence_CaLeak_SL = 2.0;
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
        T = 293;        %K
        KMOLE = 0.001660538783162726;

        V_a = 70.0;
        V_h = -45.0;
        V_hkinf = -40.0;
        V_m = -46.0;
        V_n = -40.0;
        V_Sinf = -78.0;
        V_tau = 90.0;
        VBar_RyR = -20.0;

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


        % Geometry
        vol_Fiber = PI * (R_fiber ^ 2) * L_fiber ;
        vol_myo = volFraction_myo * vol_Fiber ;
        SA_SL = 2 * PI * R_fiber * L_fiber ;
        vol_SR = volFraction_SR * vol_Fiber ;
        SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio ;
        TTFrac = SA_TT / SA_SL ;

        % Functions

        K_EC = K_EC_init_uM;
        E_K_K_IR = (log((K_EC ./ K_i)) .* R .* T ./ F);
        K_R = (K_EC .* exp(( - delta .* E_K_K_IR .* F ./ (R .* T))));
        I_K_IR = ((G_K .* (Q10g_KIR .^ QCorr) .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_SL .* F) ./ (R .* T))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_SL )) .* (1.0 + TTFrac));
        J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* F)) .* 1.0E9);
        Na_EC = Na_EC_init_uM;
        E_Na = (log((Na_EC ./ Na_i)) .* R .* T ./ F);
        KmNa_i_NCX = (12.29 .* 1000.0);
        volFactor = (vol_myo ./ (PI .* 0.26));
        nu_SERCA = nu_SERCA .* volFactor;  
        LumpedJ_SERCA = (Q10SERCA .^ QCorr) * (602.214179 * nu_SERCA * c_i ./ (K_SERCA + c_i));
        I_Na = (((g_Na * Q10g_Na .* ((m ./ mUnit) ^ 3.0) * (h ./ hUnit) * (S ./ SUnit)) * (1.0 + (0.1 .* TTFrac))) + g_leakNa) * (E_Na - Voltage_SL);
        J_Na = ((I_Na ./ (carrierValence_Na .* F)) .* 1E09);
        fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_SL .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_SL .* F ./ (R .* T))))) ^  - 1.0);
        I_NKX_K = ((2.0 .* (fPump_NKX_K .* F .* Q10g_NaK .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
        g_PMCA = (g_PMCA * vol_myo) / (700 * SA_SL ); 
        I_PMCA = (Q10PMCA .^ QCorr) * ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
        fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_SL .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_SL .* F ./ (R .* T))))) ^  - 1.0);
        I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* F .* Q10g_NaK .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));

        if expt == 0
            g_SOCE = 0;
        elseif expt >= 1
            g_SOCE = (0.01 ./ 210.44);
        end

        if expt == 2
            continuousStim = true;
        else
            continuousStim = false;
        end

        c_EC = c_EC_init_uM;
        E_Ca = (log((c_EC ./ c_i)) .* (R .* T) ./ (2.0 .* F));
        I_SOCE = ((g_SOCE .* (E_Ca - Voltage_SL ) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
        E_K_K_DR = (log((K_EC ./ K_i)) .* R .* T ./ F);
        I_K_DR = ((g_K .* (Q10g_K^ QCorr) .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_SL )) .* (1.0 + (0.45 .* TTFrac)));
        g_SLLeak = 5.0 .* 2.0E-6;
        I_CaLeak_SL = ((g_SLLeak .* (E_Ca - Voltage_SL )) .* (1.0 + TTFrac));
        s1_NCX_N = (exp((nu_NCX .* Voltage_SL .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_EC);
        Ka_NCX_N = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
        Qcorr_NCX = ((T - 310.0) ./ 10.0);
        s2_NCX_N = (exp(((nu_NCX - 1.0) .* Voltage_SL .* F ./ (R .* T))) .* (Na_EC ^ 3.0) .* c_i);
        KmNa_EC_NCX_N = (87.5 .* 1000.0);
        Kmc_EC_NCX_N = (1.3 .* 1000.0);
        s3_NCX_N = ((Kmc_i_NCX .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ((KmNa_EC_NCX_N ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + (Kmc_EC_NCX_N .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
        I_NCX_N = ( - (3.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_N .* (s1_NCX_N - s2_NCX_N) ./ s3_NCX_N ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* Voltage_SL ./ (R .* T ./ F))))))) .* (1.0 + TTFrac));
        A = (1.0 / (1.0 + exp(((Voltage_SL - V_a) / A_a))));
        Cl_EC = Cl_EC_init_uM;
        E_Cl =  - (log((Cl_EC ./ Cl_i)) .* R .* T ./ F);
        I_Cl = ((g_Cl .* (Q10g_Cl.^ QCorr) .* (A ^ 4.0) .* (E_Cl - Voltage_SL )) .* (1.0 + (0.1 .* TTFrac)));
        s2_NCX_C = (exp(((nu_NCX - 1.0) .* Voltage_SL .* F ./ (R .* T))) .* (Na_EC ^ 3.0) .* c_i);
        Kmc_EC_NCX_C = (1.6 .* 1000.0);
        KmNa_EC_NCX_C = (87.5 .* 1000.0);
        s3_NCX_C = ((Kmc_i_NCX .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ((KmNa_EC_NCX_C ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + (Kmc_EC_NCX_C .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
        Ka_NCX_C = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
        Qcorr_NCX = ((T - 310.0) ./ 10.0);
        s1_NCX_C = (exp((nu_NCX .* Voltage_SL .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_EC);
        I_NCX_C = ((2.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_C .* (s1_NCX_C - s2_NCX_C) ./ s3_NCX_C ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* Voltage_SL ./ (R .* T ./ F))))))) .* (1.0 + TTFrac));
        
        %tau_w_r7 = (1.0 ./ (100.0 .* (1.0 + c_i)));
        % alpha_w_r0 = 100;
        % K_w_r0 = 1;

        tau_w_r0 = (1.0 / (alpha_w_r0 * ( 1 + (c_i / K_w_r0)))) ;%(100.0 .* (1.0 + c_i)));

        alpha_n = (alpha_n0 .* (Q10alpha_n .^ QCorr) .* (Voltage_SL - V_n) ./ (1.0 - exp(( - (Voltage_SL - V_n) ./ K_alphan))));
        alpha_m = (alpha_m0 .* (Q10alpha_m .^ QCorr) .* (Voltage_SL - V_m) ./ (1.0 - exp(( - (Voltage_SL - V_m) ./ K_alpham))));
        
        L_RyR = (1000.0 ./ 0.002);
        w_DHPR = w_RyR;
        f_DHPR = f_RyR;
        K_DHPR = K_RyR;
        L_DHPR = L_RyR;
        VBar_DHPR = VBar_RyR;
        openDHPR = ((((1.0 + (exp(((Voltage_SL - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_SL - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ((1.0 + exp(((Voltage_SL - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0)))) .* (w_DHPR ./ wUnit_DHPR));
        g0_DHPR = (9.39 ./ 100.0);
        I_DHPR = ((g0_DHPR .* openDHPR .* (E_Ca - Voltage_SL )) .* TTFrac);
        voltProb = (((1.0 + (exp(((Voltage_SL - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_SL - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ((1.0 + exp(((Voltage_SL - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
                
        j0_RyR = 300.0 * volFactor;
        openProb = voltProb * (w_RyR / wUnit_RyR);
        LumpedJ_RyR = 602.214179 * j0_RyR * openProb * (c_SR - c_i);

        nu_leakSR = nu_leakSR * volFactor;  %1.1338; % 0.2
        J_CaLeak_SR = 602.214179 * nu_leakSR * (c_SR - c_i);

        alpha_h = (alpha_h0 .* (Q10alpha_h.^ QCorr) .*exp(( - (Voltage_SL - V_h) ./ K_alphah)));
        J_CaLeak_SL = ((I_CaLeak_SL ./ (carrierValence_CaLeak_SL .* F)) .* 1.0E09);
        SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
        J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
        tau_hK = (exp(( - (Voltage_SL + 40.0) ./ 25.75)) .* tauUnit);
        h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_SL - V_hkinf) ./ A_hkinf))));
        J_r5 = ((h_Kinf - h_K) ./ tau_hK);
        beta_n = (beta_n0 .* (Q10beta_n .^ QCorr).* exp(( - (Voltage_SL - V_n) ./ K_betan)));
        J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
        S_inf = (1.0 ./ (1.0 + exp(((Voltage_SL - V_Sinf) ./ A_Sinf))));
        tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_SL + V_tau) ./ 100.0) ^ 2.0))));
        J_r3 = ((S_inf - S) ./ tau_S);
        beta_m = (beta_m0 .* (Q10beta_m.^ QCorr) .* exp(( - (Voltage_SL - V_m) ./ K_betam)));
        J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
        beta_h = (Q10beta_h .^ QCorr) .*(beta_h0 ./ (1.0 + exp(( - (Voltage_SL - V_h) ./ K_betah))));
        J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
        wInf_r0 = (1.0 ./ (1.0 + (c_i / K_w_r0)));
        J_r0 = ((wInf_r0 - w_RyR) ./ tau_w_r0);
        J_NKX_N =  - ((I_NKX_N ./ (carrierValence_NKX .* F)) .* 1E09);
        J_SOCE = ((I_SOCE ./ (carrierValence_SOCE .* F)) .* 1E09);
        J_NKX_K = ((I_NKX_K ./ (carrierValence_NKX .* F)) .* 1E09);
        J_K_DR =  - ((I_K_DR ./ (carrierValence_K_DR .* F)) .* 1E09);
       
        KFlux_SL_myo = (SA_SL ./ vol_myo);
        J_DHPR = ((I_DHPR ./ (carrierValence_DHPR .* F)) .* 1E09);
        J_NCX_N = ((I_NCX_N ./ (carrierValence_NCX_N .* F)) .* 1E09);
        J_NCX_C =  - ((I_NCX_C ./ (carrierValence_NCX_C .* F)) .* 1E09);
        J_PMCA =  - ((I_PMCA ./ (carrierValence_PMCA .* F)) .* 1E09);
        device_SL.Capacitance = (C_SL .* SA_SL) .* (Q10C_SL .^ QCorr);
        J_Cl = ((I_Cl ./ (carrierValence_Cl .* F)) .* 1E09);



        % Calcium buffering in the myoplasm and SR added code to VCell function

        k_onATP = 0.01364*1000;     %(uM s)^-1
        k_offATP = 30*1000;         %s^-1
        k_onParvCa = 41.7;          %(uM s)^-1
        k_offParvCa = 0.5;          %s^-1
        k_onParvMg = 0.033;         %(uM s)^-1
        k_offParvMg = 3;            %s^-1
        Parv_itot = 1500;           %(uM)
        ATP_itot = 8000;            %(uM)

        ATP = ATP_itot - CATP;
        Parv = Parv_itot - CaParv - MgParv;
        Mg = 1000;                  %(uM) constant concentration

        Trop_tot = 140;
        Trop = Trop_tot - CaTrop;
        k_offTrop = 0.115;
        k_onTrop = 0.0885;
        dCT = k_onTrop*c_i*Trop - k_offTrop*CaTrop;

        dCP = k_onParvCa*c_i*Parv - k_offParvCa*CaParv;
        dMP = k_onParvMg*Mg * Parv - k_offParvMg*MgParv;
        dCA = k_onATP*c_i*ATP - k_offATP*CATP;

        B_SRtot = 31000;
        K_SRBuffer = 800;           % k_off/k_on
        f_SR = 1/(1 + B_SRtot*K_SRBuffer./((K_SRBuffer+c_SR).^2));

        currtime = toc(StartTimer);     

        % prescribed stimulus - 500ms at 50Hz every 2.5s (applied current
        % I_SL at frequency freq - square pulses of width 1 ms)
        
        I_SL = 0;
        if (freq > 0 && t > 0)
            if continuousStim
                if (mod(t,1/freq) < pulsewidth)
                    I_SL = - ClampCurrent;
                end
            else
                period = 2.5;
                numPeriods = floor(t / period);
                timeInCurrentPeriod = t - numPeriods * period;
                if timeInCurrentPeriod <= 0.5
                    if (mod(timeInCurrentPeriod,1/freq) < pulsewidth)
                        I_SL = - ClampCurrent;
                    end
                end
            end
        end

        %% Rates

        Nf = [1;1000;1;1;100;1000;1000;0.1;1;1;1;1;100000;500;1000;1;1]; %Normalization factor
        if freq == 0 && currtime > 60
            dydt = zeros(17,1);
            fluxes = zeros(1,8);
            currents = zeros(1,13);
            return;
        end
        dydt = [
            J_r6;    % rate for SOCEProb(1)
            (f_SR * KMOLE * (LumpedJ_SERCA - LumpedJ_RyR - J_CaLeak_SR))/vol_SR; %c_SR (2)
            J_r5;    % rate for h_K (3)
            J_r0;    % rate for w_RyR (4)           
            (1000 /device_SL.Capacitance) * (SA_SL * (I_CaLeak_SL + I_Cl + I_DHPR + I_K_DR + I_K_IR + I_NCX_C + I_NCX_N + I_NKX_K + I_NKX_N + I_Na + I_PMCA + I_SOCE) + I_SL);    % rate for Voltage_SL (5)
            KFlux_SL_myo * (J_Na - J_NKX_N + J_NCX_N);    % rate for Na_i (6)
            (J_Cl .* KFlux_SL_myo);    % rate for Cl_i (7)           
            (KFlux_SL_myo * (J_SOCE + J_CaLeak_SL - J_NCX_C + J_DHPR - J_PMCA)) + ((LumpedJ_RyR - LumpedJ_SERCA + J_CaLeak_SR) * KMOLE / vol_myo) - dCP - dCA; % rate for c_i (8)
            J_r4;    % rate for n (9)
            J_r2;    % rate for m (10)
            J_r1;    % rate for h (11)
            J_r3;    % rate for S (12)
            KFlux_SL_myo * ( J_K_IR - J_K_DR + J_NKX_K);    % rate for K_i  (13)
            dCP;     % Rate for Parvalbumin bound Ca2+ (14)
            dMP;     % Rate for Parvalbumin bound Mg2+ (15)
            dCA;     % Rate for ATP bound Ca2+ (16)
            dCT;     % Rate for Trop bound Ca2+ (17)
            ];

        R = abs(dydt ./ Nf); 
        if all(R < 0.00001) && freq==0
           dydt = zeros(17,1);
           fluxes = zeros(1,8);
           currents = zeros(1,13);
           return
        end

        fluxes = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR];
        % convert all fluxes to uM/s
        fluxes(1:5) = fluxes(1:5) * KFlux_SL_myo;
        fluxes(6:8) = fluxes(6:8) * KMOLE / vol_myo;
        currents = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
        currents(1:end-1) = SA_SL * currents(1:end-1);
      
    end
end

