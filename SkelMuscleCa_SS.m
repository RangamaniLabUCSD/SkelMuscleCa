function yInf = SkelMuscleCa_SS(tSpan,lowATP, p, yinit)

% function [T,Y,yInf] = SkelMuscleCa_SS(tSpan,freq, lowATP, p, varargin)
% if isempty(varargin)
%     solve for ss values
% elseif length(varargin)==1
%   yinit = varargin{1};
% else
%   error('too many inputs')

% input:
%     tSpan is a vector of start and stop times
%     lowATP is true for low ATP conditions
%     p is a vector of selected parameters to test
%     yinit is a vector of initial conditions for the state variables
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yInf is the predicted SS values



% solve for SS values

if size(p,1) == 1
    p = p';
end

lb = zeros(17,1);
lb(6) = -120;
ub = 1e6 * ones(17,1);
%ub(13) = 1;
yInf = lsqnonlin(@(y)f(tSpan,y,p,lowATP), yinit, lb, ub); 
if any(yInf([1:5,7:17]) < 0)
    error('Variable should be non negative')
end
% -------------------------------------------------------
% ode rate
    function dydt = f(tSpan,y, p, lowATP) %f(y, p, lowATP)
	% State Variables
	SOCEProb = y(1);
	c_SR = y(2);
	h_K = y(3);
	w_RyR = y(4);
	w_DHPR = y(5);
	Voltage_PM = y(6);
	Na_i = y(7);
	Cl_i = y(8);
	c_i = y(9);
	n = y(10);
	m = y(11);
	h = y(12);
	S = y(13);
	K_i = y(14);
    CaParv = y(15);
    MgParv = y(16);
    CATP = y(17);
	
    
    % Parameters
    A_a = p(1);
    A_hkinf = p(2);
    A_Sinf = p(3);
    alpha_h0 = p(4);
    alpha_m0 = p(5);
    alpha_n0 = p(6);
    beta_h0 = p(7);
    beta_m0 = p(8);
    beta_n0 = p(9);
    C_PM = p(10);
    c_ref = p(11);
    C_SRM = p(12);
    delta = p(13);
    ClampCurrent = p(14);
    f_DHPR = p(15);
    f_RyR = p(16);
    g_Cl = p(17);
    g_K = p(18);
    G_K = p(19);
    g_Na = p(20);
    g_NCX = p(21);
    J_NaK_NKX = p(22);
    K_alphah = p(23);
    K_alpham = p(24);
    K_alphan = p(25);
    K_betah = p(26);
    K_betam = p(27);
    K_betan = p(28);
    K_DHPR = p(29);
    K_K = p(30);
    K_mk_NKX = p(31); 
    K_mNa_NKX = p(32);
    K_PMCA = p(33);
    K_RyR = p(34);
    K_S = p(35);
    K_SERCA = p(36);
    Kdact_NCX = p(37);
    Kmc_i_NCX = p(38);
    ksat_NCX = p(39);
    nu_NCX = p(40);
    Q10NCX = p(41);
    R_fiber = p(42);
    S_i = p(43);
    tau_SOCEProb = p(44);
    vol_SA_ratio = p(45);
    volFraction_cyto = p(46);
    volFraction_SR = p(47);
    volFraction_TT = p(48);
    

    % Constants
    c_EC_init_uM = 1300.0;
    Cl_EC_init_uM = 128000.0;
    K_EC_init_uM = 4000.0;
    Na_EC_init_uM = 147000.0;

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
    T = 300.0;
    KMOLE = 0.001660538783162726;

    V_a = 70.0;
    V_h = -45.0;
    V_hkinf = -40.0;
    V_m = -46.0;
    V_n = -40.0;
    V_Sinf = -78.0;
    V_tau = 90.0;
    VBar_DHPR = -20;
    VBar_RyR = -20.0;

	
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
	I_K_IR = ((G_K .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_PM .* F) ./ (R .* T))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_PM)) .* (1.0 + TTFrac));
	J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* F)) .* 1.0E9);
	Na_EC = Na_EC_init_uM;
	E_Na = (log((Na_EC ./ Na_i)) .* R .* T ./ F);
	KmNa_i_NCX = (12.29 .* 1000.0);
	volFactor = (vol_cyto ./ (PI .* 0.26));
	nu_SERCA = (4875.0 .* volFactor);
	LumpedJ_SERCA = (602.214179 * nu_SERCA * c_i ./ (K_SERCA + c_i)); 
	I_Na = ((g_Na * ((m ./ mUnit) ^ 3.0) * (h ./ hUnit) * (S ./ SUnit) * (E_Na - Voltage_PM)) * (1.0 + (0.1 .* TTFrac)));
	J_Na = ((I_Na ./ (carrierValence_Na .* F)) .* 1E09);
	fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
	I_NKX_K = ((2.0 .* (fPump_NKX_K .* F .* J_NaK_NKX ./ (((1.0 + (K_mk_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_PMCA = (5.37 ./ SA_PM);
	I_PMCA = ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
	fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
	I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* F .* J_NaK_NKX ./ (((1.0 + (K_mk_NKX ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_SOCE = (0.01 ./ 210.44);
	c_EC = c_EC_init_uM;
	E_Ca = (log((c_EC ./ c_i)) .* (R .* T) ./ (2.0 .* F));
	I_SOCE = ((g_SOCE .* (E_Ca - Voltage_PM) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
	E_K_K_DR = (log((K_EC ./ K_i)) .* R .* T ./ F);
	I_K_DR = ((g_K .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_PM)) .* (1.0 + (0.45 .* TTFrac)));
	g_PMLeak = 5.0 .* 2.0E-6;
	I_CaLeak_PM = ((g_PMLeak .* (E_Ca - Voltage_PM)) .* (1.0 + TTFrac));
	L_DHPR = (1000.0 ./ 0.002);
	openDHPR = ((((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ((1.0 + exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0)))) .* (w_DHPR ./ wUnit_DHPR));
	g0_DHPR = (9.39 ./ 100.0);
	I_DHPR = ((g0_DHPR .* openDHPR .* (E_Ca - Voltage_PM)) .* TTFrac);
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
	I_Cl = 0; % ((g_Cl .* (A ^ 4.0) .* (E_Cl - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
	s2_NCX_C = (exp(((nu_NCX - 1.0) .* Voltage_PM .* F ./ (R .* T))) .* (Na_EC ^ 3.0) .* c_i);
	Kmc_EC_NCX_C = (1.3 .* 1000.0);
	KmNa_EC_NCX_C = (87.5 .* 1000.0);
	s3_NCX_C = ((Kmc_i_NCX .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX) ^ 3.0))) + ((KmNa_EC_NCX_C ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX))) + (Kmc_EC_NCX_C .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
	Ka_NCX_C = (1.0 ./ (1.0 + ((Kdact_NCX ./ c_i) ^ 2.0)));
	Qcorr_NCX = ((T - 310.0) ./ 10.0);
	s1_NCX_C = (exp((nu_NCX .* Voltage_PM .* F ./ (R .* T))) .* (Na_i ^ 3.0) .* c_EC);
	I_NCX_C = ((2.0 .* (g_NCX .* (Q10NCX ^ Qcorr_NCX) .* Ka_NCX_C .* (s1_NCX_C - s2_NCX_C) ./ s3_NCX_C ./ (1.0 + (ksat_NCX .* exp(((nu_NCX - 1.0) .* Voltage_PM ./ (R .* T ./ F))))))) .* (1.0 + TTFrac));
	L_RyR = (1000.0 ./ 0.002);
	voltProb = (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ((1.0 + exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
	tau_w_r7 = (1.0 ./ (100.0 .* (1.0 + c_i)));
	tau_w_r0 = (1.0 ./ (100.0 .* (1.0 + c_i)));

	alpha_n = (alpha_n0 .* (Voltage_PM - V_n) ./ (1.0 - exp(( - (Voltage_PM - V_n) ./ K_alphan))));
	alpha_m = (alpha_m0 .* (Voltage_PM - V_m) ./ (1.0 - exp(( - (Voltage_PM - V_m) ./ K_alpham))));
	
    j0_RyR = 300.0 * volFactor;
	openProb = voltProb * w_RyR;
	LumpedJ_RyR = 602.214179 * j0_RyR * openProb * (c_SR - c_i);

    nu_leakSR = 0.2*volFactor;  %1.1338;
    J_CaLeak_SR = 602.214179 * nu_leakSR * (c_SR - c_i);

	alpha_h = (alpha_h0 .* exp(( - (Voltage_PM - V_h) ./ K_alphah)));
	J_CaLeak_PM = ((I_CaLeak_PM ./ (carrierValence_CaLeak_PM .* F)) .* 1.0E09);
	wInf_r7 = (1.0 ./ (1.0 + c_i));
	J_r7 = ((wInf_r7 - w_DHPR) ./ tau_w_r7);
	SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
	J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
	tau_hK = (exp(( - (Voltage_PM + 40.0) ./ 25.75)) .* tauUnit);
	h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_hkinf) ./ A_hkinf))));
	J_r5 = ((h_Kinf - h_K) ./ tau_hK);
	beta_n = (beta_n0 .* exp(( - (Voltage_PM - V_n) ./ K_betan)));
	J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
	S_inf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_Sinf) ./ A_Sinf))));
	tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_PM + V_tau) ./ 100.0) ^ 2.0))));
	J_r3 = ((S_inf - S) ./ tau_S);
	beta_m = (beta_m0 .* exp(( - (Voltage_PM - V_m) ./ K_betam)));
	J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
	beta_h = (beta_h0 ./ (1.0 + exp(( - (Voltage_PM - V_h) ./ K_betah))));
	J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
	wInf_r0 = (1.0 ./ (1.0 + c_i));
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
	device_PM.Capacitance = (C_PM .* SA_PM);
	J_Cl = ((I_Cl ./ (carrierValence_Cl .* F)) .* 1E09);
    

    %% Emmet - added code to VCell function
    % calcium buffering in the cytosol and SR
   
    k_onATP = 0.01364*1000;
    k_offATP = 30*1000;
    k_onParv = 0.0417*1000;
    k_offParv = 0.0005*1000;
    Parv_itot = 1500; 
    ATP_itot = 8000;
    % CaParv = 0.258 * Parv_itot; %Fraction of Parv sites bound with Ca2+ at resting
    % CATP = 0*ATP_itot;
    % MgParv = 0.680 * Parv_itot %Fraction of Parv sites bound with Mg2+ at resting
    % CATP = 0*ATP_itot;
    ATP = ATP_itot - CATP; 
    Parv = Parv_itot - CaParv - MgParv;
    Mg = 0.001; %1mM constant concentration
    
    dCP = k_onParv*c_i*Parv - k_offParv*CaParv;
    dMP = k_onParv*Mg * Parv - k_offParv*MgParv;
    dCA = k_onATP*c_i*ATP - k_offATP*CATP;
    %K_iATP = k_offATP/k_onATP; % ATP binding affinity
    %K_iParv = k_offParv/k_onParv; % Parvalbumin binding afinity
    %f_i = 1/(1 + (K_iATP*ATP_itot/((K_iATP+c_i).^2)) + (K_iParv*Parv_itot/((K_iParv+c_i).^2))); % Assuming rapid buffering


    B_SRtot = 31000;
    K_SRBuffer = 500;
    f_SR = 1/(1 + B_SRtot*K_SRBuffer./((K_SRBuffer+c_SR).^2));
    
    % prescribed stimulus (applied current I_PM at frequency freq - square
    % pulses of width 1 ms)
    I_PM = - ClampCurrent;
 
    % if t>0.1 && (mod(t,1/freq) < .0005 || mod(t,1/freq) > ((1/freq)-.0005))
    %     I_PM = - ClampCurrent;
    % else
    %     I_PM = 0;
    % end
   

    % switch lowATP condition on to test predictions in energy deficit (all
    % pumps go to half activity)
    if lowATP
        LumpedJ_SERCA = LumpedJ_SERCA*0.5;
        J_PMCA = J_PMCA*0.5;
        I_PMCA = I_PMCA*0.5;
        J_NKX_K = J_NKX_K*0.5;
        J_NKX_N = J_NKX_N*0.5;
        I_NKX_K = I_NKX_K*0.5;
        I_NKX_N = I_NKX_N*0.5;
    end
    
	%% Rates
	dydt = [
		J_r6;    % rate for SOCEProb
		(f_SR * KMOLE * (LumpedJ_SERCA - LumpedJ_RyR - J_CaLeak_SR))/vol_SR; %c_SR
        J_r5;    % rate for h_K
		J_r0;    % rate for w_RyR
		J_r7;    % rate for w_DHPR
		(1000 * SA_PM /device_PM.Capacitance) * (I_CaLeak_PM + I_Cl +...
        I_DHPR + I_K_DR + I_K_IR + I_NCX_C + I_NCX_N + I_NKX_K + I_NKX_N ...
        + I_Na + I_PMCA + I_SOCE + I_PM);    % rate for Voltage_PM
		KFlux_PM_cyto * (J_Na - J_NKX_N + J_NCX_N);    % rate for Na_i
		0;%(J_Cl .* KFlux_PM_cyto);    % rate for Cl_i
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
    dydt = dydt ./ yinit;

    end
end 

