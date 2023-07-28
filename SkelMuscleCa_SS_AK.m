function [T,Y,SSParam] = SkelMuscleCa_AK(argTimeSpan, freq, lowATP)
% [T,Y,yinit,param] = SkelMuscleCa_AK(argTimeSpan,argYinit,argParam)
%
% input:
%     argTimeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
%     argYinit is a vector of initial conditions for the state variables (optional)
%     argParam is a vector of values for the parameters (optional)
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yinit is the initial conditions that were used
%     param is the parameter vector that was used
%     allNames is the output solution variable names
%     allValues is the output solution variable values corresponding to the names
%
%     example of running this file: [T,Y,yinit,param,allNames,allValues] = myMatlabFunc; <-(your main function name)
%

%
% Default time span
%
timeSpan = [0.0 1.0];

if nargin >= 1
	if length(argTimeSpan) > 0
		%
		% TimeSpan overridden by function arguments
		%
		timeSpan = argTimeSpan;
	end
end
%
% Default Initial Conditions
%
yinit = [
	0.0122; %4.743E-6;		% yinit(1) is the initial condition for 'SOCEProb'
	1500.0;		% yinit(2) is the initial condition for 'c_SR'
	0.9983;		% yinit(3) is the initial condition for 'h_K'
	0.9091;		% yinit(4) is the initial condition for 'w_RyR'
	0.9091;		% yinit(5) is the initial condition for 'w_DHPR'
	-88.0;		% yinit(6) is the initial condition for 'Voltage_PM'
	14700.0;		% yinit(7) is the initial condition for 'Na_i'
	5830.0;		% yinit(8) is the initial condition for 'Cl_i'
	0.1;		% yinit(9) is the initial condition for 'c_i'
	0.003;		% yinit(10) is the initial condition for 'n'
	0.0128;		% yinit(11) is the initial condition for 'm'
	0.8051;		% yinit(12) is the initial condition for 'h'
	0.8487;		% yinit(13) is the initial condition for 'S'
	154500.0;		% yinit(14) is the initial condition for 'K_i'
];

%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%


%% Parameters 
param = [

    288.0; % param(1) is 'alpha_m0'
    0.17; % param(2) is 'K_PMCA'
    1.0; % param(3) is 'wUnit_RyR'
    0.4; % param(4) is 'delta'
    13000.0; % param(5) is 'K_mNa_NKX'
    0.2; % param(6) is 'f_RyR'
    3.5420E-5; % param(7) is 'J_NaK_NKX'
    0.01; % param(8) is 'tau_SOCEProb'
    7.0; % param(9) is 'K_alphan'
    10.0; % param(10) is 'K_alpham'
    14.7; % param(11) is 'K_alphah'
    500.0; % param(12) is 'c_ref'
    1.0; % param(13) is 'SUnit'
    7.5; % param(14) is 'A_hkinf'
    1.0; % param(15) is 'SOCEProbUnit'
    0.1965; % param(16) is 'g_Cl'
    4380.0; % param(17) is 'beta_h0'
    1.0; % param(18) is 'hUnit'
    0.129; % param(19) is 'g_NCX'
    5.8; % param(20) is 'A_Sinf'
    8.04; % param(21) is 'g_Na'
    6.59; % param(22) is 'Kmc_i_NCX'
    8.1; % param(23) is 'alpha_h0'
    0.01; % param(24) is 'C_PM'
    1000.0; % param(25) is 'K_mK_NKX_N'
    1.0; % param(26) is 'h_KUnit'
    1.0; % param(27) is 'nUnit'
    1.0; % param(28) is 'C_SRM'
    1.57; % param(29) is 'Q10NCX'
    0.648; % param(30) is 'g_K'
    4.5; % param(31) is 'K_RyR'
    1.0; % param(32) is 'wUnit_DHPR'
    67.0; % param(33) is 'beta_n0'
    0.05; % param(34) is 'Kdact_NCX'
    10000.0; % param(35) is 'S_i'
    0.27; % param(36) is 'nu_NCX'
    40.0; % param(37) is 'K_betan'
    18.0; % param(38) is 'K_betam'
    9.0; % param(39) is 'K_betah'
    1380.0; % param(40) is 'beta_m0'
    1.0; % param(41) is 'mUnit'
    0.2; % param(42) is 'f_DHPR'
    4.5; % param(43) is 'K_DHPR'
    1000000.0; % param(44) is 'K_S'
    -20000.0; % param(45) is 'device_totalCurrClampElectrode.F'
    9.5E8; % param(46) is 'K_K'
    13.1; % param(47) is 'alpha_n0'
    0.32; % param(48) is 'ksat_NCX'
    150.0; % param(49) is 'A_a'
    0.111; % param(50) is 'G_K'
    1.0; % param(51) is 'K_SERCA'
    1.0; % param(52) is 'tauUnit'

	
];
%
% invoke the integrator
%
options = odeset('RelTol',1e-8,'MaxStep',.001,'NonNegative',[1:5,7:14]);%,'OutputFcn',@odeplot);
% solve for SS Param
[~,SSParam] = f(0,yinit,param,yinit,freq,lowATP);

[T,Y] = ode15s(@f,timeSpan,yinit,options,param,yinit,freq,lowATP); %pass extra arguments at the end


end

% -------------------------------------------------------
% ode rate
function [dydt,SSParam] = f(t,y,p,y0,freq,lowATP)
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
	
    % Constants    
    Size_cyto = 119381.0;
    carrierValence_K_IR = 1.0;
    Size_SRM = 37698.0;
    VBar_RyR = -20.0;
    V_tau = 90.0;
    Voltage_SRM = 0.0;
    h_K_init_molecules_um_2 = 0.9983;
    Na_i_init_uM = 14700.0;
    carrierValence_Cl = -1.0;
    V_n = -40.0;
    V_m = -46.0;
    carrierValence_NCX_N = 1.0;
    V_h = -45.0;
    V_a = 70.0;
    carrierValence_NCX_C = 2.0;
    Cl_EC_init_uM = 128000.0;
    c_i_init_uM = 0.1;
    h_init_molecules_um_2 = 0.8051;
    Na_EC_init_uM = 147000.0;
    carrierValence_Na = 1.0;
    c_EC_init_uM = 1300.0;
    Size_SR = 6283.0;
    netValence_r7 = 1.0;
    carrierValence_SOCE = 2.0;
    netValence_r6 = 1.0;
    netValence_r5 = 1.0;
    carrierValence_K_DR = 1.0;
    netValence_r4 = 1.0;
    netValence_r3 = 1.0;
    netValence_r2 = 1.0;
    netValence_r1 = 1.0;
    netValence_r0 = 1.0;
    K_i_init_uM = 154500.0;
    Voltage_SRM_init = 0.0;
    Size_EC = 1000000.0;
    Fnmol_ = 9.64853321E-5;
    T = 300.0;
    K_millivolts_per_volt = 1000.0;
    carrierValence_CaLeak_SR = 1.0;
    carrierValence_SERCA = 1.0;
    Cl_i_init_uM = 5830.0;
    Size_PM = 12566.0;
    carrierValence_DHPR = 2.0;
    c_SR_init_uM = 1500.0;
    VBar_DHPR = -20;
    PI = 3.141592653589793;
    V_Sinf = -78.0;
    F = 96485.3321;
    w_DHPR_init_molecules_um_2 = 0.0;
    R = 8314.46261815;
    S_init_molecules_um_2 = 0.8487;
    mlabfix_K_GHK_ = 1.0E-9;
    carrierValence_RyR = 1.0;
    V_hkinf = -40.0;
    carrierValence_PMCA = 2.0;
    K_EC_init_uM = 4000.0;
    N_pmol = 6.02214179E11;
    carrierValence_CaLeak_PM = 2.0 %1.0;
    n_init_molecules_um_2 = 0.003;
    m_init_molecules_um_2 = 0.0128;
    Voltage_PM_init = -88.0;
    TTFrac = 3.75;
    UnitFactor_uM_um3_molecules_neg_1 = 0.001660538783162726;
    w_RyR_init_molecules_um_2 = 0.9091;
    SOCEProb_init_molecules_um_2 = 4.743E-6;
    carrierValence_NKX = 1.0;


    % Parameters
    alpha_m0 = p(1);
    K_PMCA = p(2);
    wUnit_RyR = p(3);
    delta = p(4);
    K_mNa_NKX = p(5);
    f_RyR = p(6);
    J_NaK_NKX = p(7);
    tau_SOCEProb = p(8);
    K_alphan = p(9);
    K_alpham = p(10);
    K_alphah = p(11);
    c_ref = p(12);
    SUnit = p(13);
    A_hkinf = p(14);
    SOCEProbUnit = p(15);
    g_Cl = p(16);
    beta_h0 = p(17);
    hUnit = p(18);
    g_NCX = p(19);
    A_Sinf = p(20);
    g_Na = p(21);
    Kmc_i_NCX = p(22);
    alpha_h0 = p(23);
    C_PM = p(24);
    K_mK_NKX_N = p(25);
    h_KUnit = p(26);
    nUnit = p(27);
    C_SRM = p(28);
    Q10NCX = p(29);
    g_K = p(30);
    K_RyR = p(31);
    wUnit_DHPR = p(32);
    beta_n0 = p(33);
    Kdact_NCX = p(34);
    S_i = p(35);
    nu_NCX = p(36);
    K_betan = p(37);
    K_betam = p(38);
    K_betah = p(39);
    beta_m0 = p(40);
    mUnit = p(41);
    f_DHPR = p(42);
    K_DHPR = p(43);
    K_S = p(44);
    device_totalCurrClampElectrode.F = p(45);
    K_K = p(46);
    alpha_n0 = p(47);
    ksat_NCX = p(48);
    A_a = p(49);
    G_K = p(50);
    K_SERCA = p(51);
    tauUnit = p(52);
	
	% Functions
	
	K_EC = K_EC_init_uM;
	E_K_K_IR = (log((K_EC ./ K_i)) .* R .* T ./ F);
	K_R = (K_EC .* exp(( - delta .* E_K_K_IR .* F ./ (R .* T))));
	I_K_IR = ((G_K .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_PM .* F) ./ (R .* T))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_PM)) .* (1.0 + TTFrac));
	%unitFactor_K_IR = (1.0E9 ./ 1.0);
	J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* F)) .* 1.0E9);
	Na_EC = Na_EC_init_uM;
	E_Na = (log((Na_EC ./ Na_i)) .* R .* T ./ F);
	%unitFactor_CaLeak_PM = (1.0E9 ./ 1.0);
	KmNa_i_NCX = (12.29 .* 1000.0);
	% KmNa_i_NCX = (12.29 .* 1000.0);
	volFactor = (Size_cyto ./ (PI .* 0.26));
	nu_SERCA = (4875.0 .* volFactor ./ Size_SRM);
	LumpedJ_SERCA = (602.214179 .* Size_SRM .* nu_SERCA .* c_i ./ (K_SERCA + c_i));
	I_Na = ((g_Na .* ((m ./ mUnit) ^ 3.0) .* (h ./ hUnit) .* (S ./ SUnit) .* (E_Na - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
	% unitFactor_Na = (1.0E9 ./ 1.0);
	J_Na = ((I_Na ./ (carrierValence_Na .* F)) .* 1E09);
	unitFactor_SOCE = (1.0E9 ./ 1.0);
	unitFactor_K_DR = (1.0E9 ./ 1.0);
	fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
	I_NKX_K = ((2.0 .* (fPump_NKX_K .* F .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX_K ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_PMCA = (5.37 ./ Size_PM);
	I_PMCA = ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
	fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* F ./ (R .* T)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* F ./ (R .* T))))) ^  - 1.0);
	I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* F .* J_NaK_NKX ./ (((1.0 + (K_mK_NKX_N ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_SOCE = (0.01 ./ 210.44);
	c_EC = c_EC_init_uM;
	E_Ca = (log((c_EC ./ c_i)) .* (R .* T) ./ (2.0 .* F));
	I_SOCE = ((g_SOCE .* (E_Ca - Voltage_PM) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
	E_K_K_DR = (log((K_EC ./ K_i)) .* R .* T ./ F);
	I_K_DR = ((g_K .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_PM)) .* (1.0 + (0.45 .* TTFrac)));
	g_PMLeak = 3.542E-4; %3.0455E-5; %2.1362E-6 ; %(5.0 .* 2.0E-6);
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
	A = (1.0 ./ (1.0 + exp(((Voltage_PM - V_a) ./ A_a))));
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
	F_PM = ( - (I_DHPR .* Size_PM) - (I_PMCA .* Size_PM) - (I_K_DR .* Size_PM) - (I_K_IR .* Size_PM) - (I_NKX_K .* Size_PM) - (I_Na .* Size_PM) - (I_Cl .* Size_PM) - (I_NCX_C .* Size_PM) - (I_NKX_N .* Size_PM) - (I_NCX_N .* Size_PM) - (I_SOCE .* Size_PM) - (I_CaLeak_PM .* Size_PM));
	L_RyR = (1000.0 ./ 0.002);
	voltProb = (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ((1.0 + exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
	tau_w_r7 = (1.0 ./ (100.0 .* (1.0 + c_i)));
	tau_w_r0 = (1.0 ./ (100.0 .* (1.0 + c_i)));
	j0_RyR = (300.0 .* volFactor ./ Size_SRM);
	openProb = (voltProb .* (w_RyR ./ wUnit_RyR));
	LumpedJ_RyR = (602.214179 .* Size_SRM .* j0_RyR .* openProb .* (c_SR - c_i));
	alpha_n = (alpha_n0 .* (Voltage_PM - V_n) ./ (1.0 - exp(( - (Voltage_PM - V_n) ./ K_alphan))));
	alpha_m = (alpha_m0 .* (Voltage_PM - V_m) ./ (1.0 - exp(( - (Voltage_PM - V_m) ./ K_alpham))));
	nu_leakSR = 1.1338; %(0.02 .* volFactor .* 10.0 ./ Size_SRM);
	LumpedJ_CaLeak_SR = (602.214179 .* Size_SRM .* nu_leakSR .* (c_SR - c_i));
	alpha_h = (alpha_h0 .* exp(( - (Voltage_PM - V_h) ./ K_alphah)));
	J_CaLeak_PM = ((I_CaLeak_PM ./ (carrierValence_CaLeak_PM .* F)) .* 1E9);
	wInf_r7 = (1.0 ./ (1.0 + c_i));
	J_r7 = ((wInf_r7 - w_DHPR) ./ tau_w_r7);
	% unitFactor_NKX = (1.0E9 ./ 1.0);
	SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
	J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
	tau_hK = (exp(( - (Voltage_PM + 40.0) ./ 25.75)) .* tauUnit);
	h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_hkinf) ./ A_hkinf))));
	J_r5 = ((h_Kinf - h_K) ./ tau_hK);
	beta_n = (beta_n0 .* exp(( - (Voltage_PM - V_n) ./ K_betan)));
	J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
	% unitFactor_NKX = (1.0E9 ./ 1.0);
	S_inf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_Sinf) ./ A_Sinf))));
	tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_PM + V_tau) ./ 100.0) ^ 2.0))));
	J_r3 = ((S_inf - S) ./ tau_S);
	beta_m = (beta_m0 .* exp(( - (Voltage_PM - V_m) ./ K_betam)));
	J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
	beta_h = (beta_h0 ./ (1.0 + exp(( - (Voltage_PM - V_h) ./ K_betah))));
	J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
	wInf_r0 = (1.0 ./ (1.0 + c_i));
	J_r0 = ((wInf_r0 - w_RyR) ./ tau_w_r0);
	% unitFactor_DHPR = (1.0E9 ./ 1.0);
	J_NKX_N =  - ((I_NKX_N ./ (carrierValence_NKX .* F)) .* 1E09);
	J_SOCE = ((I_SOCE ./ (carrierValence_SOCE .* F)) .* 1E09);
	J_NKX_K = ((I_NKX_K ./ (carrierValence_NKX .* F)) .* 1E09);
	J_K_DR =  - ((I_K_DR ./ (carrierValence_K_DR .* F)) .* 1E09);
	%unitFactor_RyR = (1.0 ./ N_pmol);
	LumpedI_RyR =  - ((carrierValence_RyR .* F .* LumpedJ_RyR) ./ N_pmol);
	% unitFactor_PMCA = (1.0E9 ./ 1.0);
	% unitFactor_NCX_N = (1.0E9 ./ 1.0);
	% unitFactor_NCX_C = (1.0E9 ./ 1.0);
	KFlux_SRM_SR = (Size_SRM ./ Size_SR);
	KFlux_PM_cyto = (Size_PM ./ Size_cyto);
	J_DHPR = ((I_DHPR ./ (carrierValence_DHPR .* F)) .* 1E09);
	J_NCX_N = ((I_NCX_N ./ (carrierValence_NCX_N .* F)) .* 1E09);
	UnitFactor_mV_pF_s_neg_1_pA_neg_1 = (1000.0 ./ 1.0);
	J_NCX_C =  - ((I_NCX_C ./ (carrierValence_NCX_C .* F)) .* 1E09);
	KFlux_SRM_cyto = (Size_SRM ./ Size_cyto);
	I_PM =  - device_totalCurrClampElectrode.F;
	J_PMCA =  - ((I_PMCA ./ (carrierValence_PMCA .* F)) .* 1E09);
	device_PM.Capacitance = (C_PM .* Size_PM);
	J_Cl = ((I_Cl ./ (carrierValence_Cl .* F)) .* 1E09);

    % Emmet - added code to VCell function
    % calcium buffering in the cytosol and SR
    B_itot = 300;
    K_iBuffer = 0.5;
    B_SRtot = 31000;
    K_SRBuffer = 500;
    f_i = 1/(1 + B_itot*K_iBuffer./((K_iBuffer+c_i).^2));
    f_SR = 1/(1 + B_SRtot*K_SRBuffer./((K_SRBuffer+c_SR).^2));
    % prescribed stimulus (applied current I_PM at frequency freq - square
    % pulses of width 1 ms)
    if t>0.1 && (mod(t,1/freq) < .0005 || mod(t,1/freq) > ((1/freq)-.0005))
        I_PM = I_PM;
    else
        I_PM = 0;
    end

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

    % Steady State rates
    %KFluxDeficit = (J_K_DR - J_K_IR);
    J_NaK_SS = J_NaK_NKX * ((J_K_DR - J_K_IR) / J_NKX_K);
    J_NKX_K_SS = J_NKX_K * J_NaK_SS / J_NaK_NKX;
    g_NCX_SS = g_NCX*((J_NKX_N * (J_NaK_SS / J_NaK_NKX)) - J_Na)/ J_NCX_N;
    CaFluxDeficit = (J_NCX_C*g_NCX_SS/g_NCX) + J_PMCA - J_SOCE - J_DHPR;
    g_PMleak_SS = g_PMLeak*(CaFluxDeficit/J_CaLeak_PM);
    nu_leakSR_SS = nu_leakSR*(LumpedJ_SERCA - LumpedJ_RyR)/LumpedJ_CaLeak_SR;


    SSParam = [J_NaK_SS,g_NCX_SS,g_PMleak_SS,nu_leakSR_SS]; % J_NaK, etc.
    
	% Rates
	dydt = [
		J_r6;    % rate for SOCEProb
		f_SR*( - (KFlux_SRM_SR .* LumpedJ_RyR .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM) + (KFlux_SRM_SR .* LumpedJ_SERCA .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM) - (LumpedJ_CaLeak_SR .* KFlux_SRM_SR .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM));    % rate for c_SR
		J_r5;    % rate for h_K
		J_r0;    % rate for w_RyR
		J_r7;    % rate for w_DHPR
		((UnitFactor_mV_pF_s_neg_1_pA_neg_1 .* ((I_CaLeak_PM .* Size_PM) + (I_Cl .* Size_PM) + (I_DHPR .* Size_PM) + (I_K_DR .* Size_PM) + (I_K_IR .* Size_PM) + (I_NCX_C .* Size_PM) + (I_NCX_N .* Size_PM) + (I_NKX_K .* Size_PM) + (I_NKX_N .* Size_PM) + (I_Na .* Size_PM) + (I_PMCA .* Size_PM) + (Size_PM .* I_SOCE) + I_PM)) ./ device_PM.Capacitance);    % rate for Voltage_PM
		((KFlux_PM_cyto .* J_Na) - (KFlux_PM_cyto .* J_NKX_N) + (KFlux_PM_cyto .* J_NCX_N));    % rate for Na_i
		0;%(J_Cl .* KFlux_PM_cyto);    % rate for Cl_i
		f_i*((KFlux_SRM_cyto .* LumpedJ_RyR .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM) + (J_DHPR .* KFlux_PM_cyto) - (KFlux_PM_cyto .* J_PMCA) - (KFlux_SRM_cyto .* LumpedJ_SERCA .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM) + (LumpedJ_CaLeak_SR .* KFlux_SRM_cyto .* UnitFactor_uM_um3_molecules_neg_1 ./ Size_SRM) - (KFlux_PM_cyto .* J_NCX_C)  + (KFlux_PM_cyto .* J_SOCE) + (J_CaLeak_PM .* KFlux_PM_cyto));    % rate for c_i
		J_r4;    % rate for n
		J_r2;    % rate for m
		J_r1;    % rate for h
		J_r3;    % rate for S
		( - (KFlux_PM_cyto .* J_K_DR) + (KFlux_PM_cyto .* J_K_IR) + (KFlux_PM_cyto .* J_NKX_K));    % rate for K_i
	];

end