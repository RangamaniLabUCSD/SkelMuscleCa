function [T,Y] = SkelMuscleCa(argTimeSpan,freq, lowATP)
% [T,Y,yinit,param] = SkelMuscleCa(argTimeSpan,argYinit,argParam)
%
% input:
%     argTimeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
% additional inputs from Emmet's alterations:
%     freq: frequency of applied current stimulus in Hz
%     lowATP: flag (logical) for low energy condition
% output:
%     T is the vector of times (in seconds)
%     Y is the vector of state variables
%     yinit is the initial conditions that were used
%     param is the parameter vector that was used
%     allNames is the output solution variable names
%     allValues is the output solution variable values corresponding to the names
%
%     example of running this file: [T,Y,yinit,param,allNames,allValues] = myMatlabFunc; <-(your main function name)
%
%   Variables in Y:
%   Y(1) is 'SOCEProb' - store-operated channel open probability (dimensionless)
%   Y(2) is 'c_SR' - SR calcium conc (uM)
%   Y(3) is 'h_K' - var for potassium delayed rectifier opening (dimensionless)
%   Y(4) is 'w_RyR' - var for ryanodine receptor opening (dimensionless)
%   Y(5) is 'Voltage_PM' - membrane voltage (mV)
%   Y(6) is 'Na_i' - intracellular sodium conc (uM)
%   Y(7) is 'Cl_i' - intracellular chloride conc (uM)
%   Y(8) is 'c_i' - intracellular calcium conc (uM)
%   Y(9) is 'n' - var for potassium delayed rectifier opening (dimensionless)
%   Y(10) is 'm' - var for Na channel opening (dimensionless)
%   Y(11) is 'h' - var for Na channel opening (dimensionless)
%   Y(12) is 'S' - var for Na channel opening (dimensionless)
%   Y(13) is 'K_i' - intracellular potassium conc (uM)

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
	4.743E-6;		% yinit(1) is the initial condition for 'SOCEProb'
	1500.0;		% yinit(2) is the initial condition for 'c_SR'
	0.9983;		% yinit(3) is the initial condition for 'h_K'
	0.9091;		% yinit(4) is the initial condition for 'w_RyR'
	-88.0;		% yinit(5) is the initial condition for 'Voltage_PM'
	14700.0;		% yinit(6) is the initial condition for 'Na_i'
	5830.0;		% yinit(7) is the initial condition for 'Cl_i'
	0.1;		% yinit(8) is the initial condition for 'c_i'
	0.003;		% yinit(9) is the initial condition for 'n'
	0.0128;		% yinit(10) is the initial condition for 'm'
	0.8051;		% yinit(11) is the initial condition for 'h'
	0.8487;		% yinit(12) is the initial condition for 'S'
	154500.0;		% yinit(13) is the initial condition for 'K_i'
];

%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%
param = [
	288.0;		% param(1) is 'alpha_m0'
	1.0;		% param(2) is 'netValence_SERCA'
	1.0;		% param(3) is 'carrierValence_K_IR'
	37698.0;		% param(4) is 'Size_SRM'
	0.17;		% param(5) is 'K_PMCA'
	1.0;		% param(6) is 'wUnit'
	-20.0;		% param(7) is 'VBar_RyR'
	0.4;		% param(8) is 'delta'
	13000.0;		% param(9) is 'K_mNa_NKX_N'
	0.2;		% param(10) is 'f_RyR'
	13000.0;		% param(11) is 'K_mNa_NKX_K'
	90.0;		% param(12) is 'V_tau'
	6.21E-6;		% param(13) is 'J_NaK_NKX_N'
	0.0;		% param(14) is 'Voltage_SRM'
	6.21E-6;		% param(15) is 'J_NaK_NKX_K'
	0.01;		% param(16) is 'tau_SOCEProb'
	0.9983;		% param(17) is 'h_K_init_molecules_um_2'
	14700.0;		% param(18) is 'Na_i_init_uM'
	7.0;		% param(19) is 'K_alphan'
	10.0;		% param(20) is 'K_alpham'
	14.7;		% param(21) is 'K_alphah'
	500.0;		% param(22) is 'c_ref'
	-1.0;		% param(23) is 'carrierValence_Cl'
	1.0;		% param(24) is 'SUnit'
	7.5;		% param(25) is 'A_hkinf'
	1.0;		% param(26) is 'SOCEProbUnit'
	-40.0;		% param(27) is 'V_n'
	-46.0;		% param(28) is 'V_m'
	0.1965;		% param(29) is 'g_Cl'
	4380.0;		% param(30) is 'beta_h0'
	1.0;		% param(31) is 'carrierValence_NCX_N'
	-45.0;		% param(32) is 'V_h'
	119381.0;		% param(33) is 'Size_cyto'
	70.0;		% param(34) is 'V_a'
	2.0;		% param(35) is 'carrierValence_NCX_C'
	128000.0;		% param(36) is 'Cl_EC_init_uM'
	1.0;		% param(37) is 'hUnit'
	0.129;		% param(38) is 'g_NCX_NCX_N'
	1.0;		% param(39) is 'netValence_RyR'
	0.1;		% param(40) is 'c_i_init_uM'
	0.129;		% param(41) is 'g_NCX_NCX_C'
	0.8051;		% param(42) is 'h_init_molecules_um_2'
	147000.0;		% param(43) is 'Na_EC_init_uM'
	1.0;		% param(44) is 'carrierValence_Na'
	5.8;		% param(45) is 'A_Sinf'
	1300.0;		% param(46) is 'c_EC_init_uM'
	8.04;		% param(47) is 'g_Na'
	6283.0;		% param(48) is 'Size_SR'
	6.59;		% param(49) is 'Kmc_i_NCX_N'
	8.1;		% param(50) is 'alpha_h0'
	2.0;		% param(51) is 'carrierValence_SOCE'
	1.0;		% param(52) is 'netValence_r6'
	1.0;		% param(53) is 'netValence_r5'
	1.0;		% param(54) is 'carrierValence_K_DR'
	1.0;		% param(55) is 'netValence_r4'
	6.59;		% param(56) is 'Kmc_i_NCX_C'
	1.0;		% param(57) is 'netValence_r3'
	1.0;		% param(58) is 'netValence_r2'
	1.0;		% param(59) is 'netValence_r1'
	1.0;		% param(60) is 'netValence_r0'
	154500.0;		% param(61) is 'K_i_init_uM'
	0.01;		% param(62) is 'C_PM'
	1.0;		% param(63) is 'netValence_CaLeak_SR'
	0.0;		% param(64) is 'Voltage_SRM_init'
	1000000.0;		% param(65) is 'Size_EC'
	9.64853321E-5;		% param(66) is 'mlabfix_F_nmol_'
	1000.0;		% param(67) is 'K_mK_NKX_N'
	300.0;		% param(68) is 'mlabfix_T_'
	1000.0;		% param(69) is 'K_millivolts_per_volt'
	1000.0;		% param(70) is 'K_mK_NKX_K'
	1.0;		% param(71) is 'h_KUnit'
	1.0;		% param(72) is 'nUnit'
	1.0;		% param(73) is 'C_SRM'
	1.57;		% param(74) is 'Q10NCX_NCX_N'
	5830.0;		% param(75) is 'Cl_i_init_uM'
	0.0;		% param(76) is 'g0_DHPR' - currently zero flux through the DHPR
	12566.0;		% param(77) is 'Size_PM'
	1.57;		% param(78) is 'Q10NCX_NCX_C'
	2.0;		% param(79) is 'carrierValence_DHPR'
	1500.0;		% param(80) is 'c_SR_init_uM'
	0.648;		% param(81) is 'g_K'
	4.5;		% param(82) is 'K_RyR'
	-20.0;		% param(83) is 'VBar_DHPR'
	3.141592653589793;		% param(84) is 'mlabfix_PI_'
	-78.0;		% param(85) is 'V_Sinf'
	96485.3321;		% param(86) is 'mlabfix_F_'
	8314.46261815;		% param(87) is 'mlabfix_R_'
	0.8487;		% param(88) is 'S_init_molecules_um_2'
	67.0;		% param(89) is 'beta_n0'
	0.05;		% param(90) is 'Kdact_NCX_N'
	10000.0;		% param(91) is 'S_i'
	0.05;		% param(92) is 'Kdact_NCX_C'
	1.0E-9;		% param(93) is 'mlabfix_K_GHK_'
	0.27;		% param(94) is 'nu_NCX_N'
	-40.0;		% param(95) is 'V_hkinf'
	2.0;		% param(96) is 'carrierValence_PMCA'
	4000.0;		% param(97) is 'K_EC_init_uM'
	0.27;		% param(98) is 'nu_NCX_C'
	40.0;		% param(99) is 'K_betan'
	18.0;		% param(100) is 'K_betam'
	9.0;		% param(101) is 'K_betah'
	1380.0;		% param(102) is 'beta_m0'
	6.02214179E11;		% param(103) is 'mlabfix_N_pmol_'
	1.0;		% param(104) is 'carrierValence_CaLeak_PM'
	1.0;		% param(105) is 'mUnit'
	0.003;		% param(106) is 'n_init_molecules_um_2'
	0.2;		% param(107) is 'f_DHPR'
	4.5;		% param(108) is 'K_DHPR'
	1000000.0;		% param(109) is 'K_S'
	-20000.0;		% param(110) is 'device_totalCurrClampElectrode.F'
	9.5E8;		% param(111) is 'K_K'
	0.0128;		% param(112) is 'm_init_molecules_um_2'
	13.1;		% param(113) is 'alpha_n0'
	0.32;		% param(114) is 'ksat_NCX_N'
	-88.0;		% param(115) is 'Voltage_PM_init'
	3.75;		% param(116) is 'TTFrac' (Area fraction of T-tubules area rel to sarcolemma area)
	150.0;		% param(117) is 'A_a'
	0.111;		% param(118) is 'G_K'
	0.9091;		% param(119) is 'w_RyR_init_molecules_um_2'
	4.743E-6;		% param(120) is 'SOCEProb_init_molecules_um_2'
	1.0;		% param(121) is 'K_SERCA'
	0.001660538783162726;		% param(122) is 'KMOLE'
	0.32;		% param(123) is 'ksat_NCX_C'
	1.0;		% param(124) is 'carrierValence_NKX_N'
	1.0;		% param(125) is 'tauUnit'
	1.0;		% param(126) is 'carrierValence_NKX_K'
];
%
% invoke the integrator
%
options = odeset('RelTol',1e-8,'MaxStep',.0001,'NonNegative',[1:4,6:13]);%,'OutputFcn',@odeplot);
[T,Y] = ode15s(@f,timeSpan,yinit,options,param,yinit,freq,lowATP); %pass extra arguments at the end
end

% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0,freq,lowATP)
	% State Variables
	SOCEProb = y(1);
	c_SR = y(2);
	h_K = y(3);
	w_RyR = y(4);
	Voltage_PM = y(5);
	Na_i = y(6);
	Cl_i = y(7);
	c_i = y(8);
	n = y(9);
	m = y(10);
	h = y(11);
	S = y(12);
	K_i = y(13);
	% Constants
	alpha_m0 = p(1);
	netValence_SERCA = p(2);
	carrierValence_K_IR = p(3);
	Size_SRM = p(4);
	K_PMCA = p(5);
	wUnit = p(6);
	VBar_RyR = p(7);
	delta = p(8);
	K_mNa_NKX_N = p(9);
	f_RyR = p(10);
	K_mNa_NKX_K = p(11);
	V_tau = p(12);
	J_NaK_NKX_N = p(13);
	Voltage_SRM = p(14);
	J_NaK_NKX_K = p(15);
	tau_SOCEProb = p(16);
	h_K_init_molecules_um_2 = p(17);
	Na_i_init_uM = p(18);
	K_alphan = p(19);
	K_alpham = p(20);
	K_alphah = p(21);
	c_ref = p(22);
	carrierValence_Cl = p(23);
	SUnit = p(24);
	A_hkinf = p(25);
	SOCEProbUnit = p(26);
	V_n = p(27);
	V_m = p(28);
	g_Cl = p(29);
	beta_h0 = p(30);
	carrierValence_NCX_N = p(31);
	V_h = p(32);
	Size_cyto = p(33);
	V_a = p(34);
	carrierValence_NCX_C = p(35);
	Cl_EC_init_uM = p(36);
	hUnit = p(37);
	g_NCX_NCX_N = p(38);
	netValence_RyR = p(39);
	c_i_init_uM = p(40);
	g_NCX_NCX_C = p(41);
	h_init_molecules_um_2 = p(42);
	Na_EC_init_uM = p(43);
	carrierValence_Na = p(44);
	A_Sinf = p(45);
	c_EC_init_uM = p(46);
	g_Na = p(47);
	Size_SR = p(48);
	Kmc_i_NCX_N = p(49);
	alpha_h0 = p(50);
	carrierValence_SOCE = p(51);
	netValence_r6 = p(52);
	netValence_r5 = p(53);
	carrierValence_K_DR = p(54);
	netValence_r4 = p(55);
	Kmc_i_NCX_C = p(56);
	netValence_r3 = p(57);
	netValence_r2 = p(58);
	netValence_r1 = p(59);
	netValence_r0 = p(60);
	K_i_init_uM = p(61);
	C_PM = p(62);
	netValence_CaLeak_SR = p(63);
	Voltage_SRM_init = p(64);
	Size_EC = p(65);
	mlabfix_F_nmol_ = p(66);
	K_mK_NKX_N = p(67);
	mlabfix_T_ = p(68);
	K_millivolts_per_volt = p(69);
	K_mK_NKX_K = p(70);
	h_KUnit = p(71);
	nUnit = p(72);
	C_SRM = p(73);
	Q10NCX_NCX_N = p(74);
	Cl_i_init_uM = p(75);
	g0_DHPR = p(76);
	Size_PM = p(77);
	Q10NCX_NCX_C = p(78);
	carrierValence_DHPR = p(79);
	c_SR_init_uM = p(80);
	g_K = p(81);
	K_RyR = p(82);
	VBar_DHPR = p(83);
	mlabfix_PI_ = p(84);
	V_Sinf = p(85);
	mlabfix_F_ = p(86);
	mlabfix_R_ = p(87);
	S_init_molecules_um_2 = p(88);
	beta_n0 = p(89);
	Kdact_NCX_N = p(90);
	S_i = p(91);
	Kdact_NCX_C = p(92);
	mlabfix_K_GHK_ = p(93);
	nu_NCX_N = p(94);
	V_hkinf = p(95);
	carrierValence_PMCA = p(96);
	K_EC_init_uM = p(97);
	nu_NCX_C = p(98);
	K_betan = p(99);
	K_betam = p(100);
	K_betah = p(101);
	beta_m0 = p(102);
	mlabfix_N_pmol_ = p(103);
	carrierValence_CaLeak_PM = p(104);
	mUnit = p(105);
	n_init_molecules_um_2 = p(106);
	f_DHPR = p(107);
	K_DHPR = p(108);
	K_S = p(109);
	device_totalCurrClampElectrode.F = p(110);
	K_K = p(111);
	m_init_molecules_um_2 = p(112);
	alpha_n0 = p(113);
	ksat_NCX_N = p(114);
	Voltage_PM_init = p(115);
	TTFrac = p(116);
	A_a = p(117);
	G_K = p(118);
	w_RyR_init_molecules_um_2 = p(119);
	SOCEProb_init_molecules_um_2 = p(120);
	K_SERCA = p(121);
	KMOLE = p(122);
	ksat_NCX_C = p(123);
	carrierValence_NKX_N = p(124);
	tauUnit = p(125);
	carrierValence_NKX_K = p(126);
	% Functions
	volFactor = (Size_cyto ./ 0.26);
	nu_leakSR = (0.02 .* volFactor .* 10.0 ./ Size_SRM);
	J_CaLeak_SR = (nu_leakSR .* (c_SR - c_i));
	K_EC = K_EC_init_uM;
	E_K_K_IR = (log((K_EC ./ K_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	K_R = (K_EC .* exp(( - delta .* E_K_K_IR .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))));
	I_K_IR = ((G_K .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_PM .* mlabfix_F_) ./ (mlabfix_R_ .* mlabfix_T_))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_PM)) .* (1.0 + TTFrac));
	unitFactor_K_IR = (8.388606999999999E15 ./ 8388607.0);
	J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* mlabfix_F_)) .* unitFactor_K_IR);
	Na_EC = Na_EC_init_uM;
	E_Na = (log((Na_EC ./ Na_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	unitFactor_CaLeak_PM = (8.388606999999999E15 ./ 8388607.0);
	KmNa_i_NCX_N = (12.29 .* 1000.0);
	nu_SERCA = (4875.0 .* volFactor .* 0.682 ./ Size_SRM);
	J_SERCA = (nu_SERCA .* c_i ./ (K_SERCA + c_i));
	KmNa_i_NCX_C = (12.29 .* 1000.0);
	I_Na = ((g_Na .* ((m ./ mUnit) ^ 3.0) .* (h ./ hUnit) .* (S ./ SUnit) .* (E_Na - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
	unitFactor_Na = (8.388606999999999E15 ./ 8388607.0);
	J_Na = ((I_Na ./ (carrierValence_Na .* mlabfix_F_)) .* unitFactor_Na);
	unitFactor_SOCE = (8.388606999999999E15 ./ 8388607.0);
	unitFactor_K_DR = (8.388606999999999E15 ./ 8388607.0);
	fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67.3)) - 1.0) .* exp(( - Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))))) ^  - 1.0);
	I_NKX_K = ((2.0 .* (fPump_NKX_K .* mlabfix_F_ .* J_NaK_NKX_K ./ (((1.0 + (K_mK_NKX_K ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX_K ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_PMCA = (5.37 ./ Size_PM);
	I_PMCA = ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
	fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67.3)) - 1.0) .* exp(( - Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))))) ^  - 1.0);
	I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* mlabfix_F_ .* J_NaK_NKX_N ./ (((1.0 + (K_mK_NKX_N ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX_N ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_SOCE = (0.01 ./ 210.44);
	c_EC = c_EC_init_uM;
	E_Ca = (log((c_EC ./ c_i)) .* (mlabfix_R_ .* mlabfix_T_) ./ (2.0 .* mlabfix_F_));
	I_SOCE = ((g_SOCE .* (E_Ca - Voltage_PM) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
	E_K_K_DR = (log((K_EC ./ K_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	I_K_DR = ((g_K .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_PM)) .* (1.0 + (0.45 .* TTFrac)));
	g_PMLeak = (5.0 .* 2.0E-6);
	I_CaLeak_PM = ((g_PMLeak .* (E_Ca - Voltage_PM)) .* (1.0 + TTFrac));
	L_DHPR = (1000.0 ./ 0.002);
	openDHPR = (((exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0) .* ((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 3.0)) + (L_DHPR .* exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* ((1.0 + exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 3.0))) ./ (((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ((1.0 + exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0))));
	I_DHPR = ((g0_DHPR .* openDHPR .* (E_Ca - Voltage_PM)) .* TTFrac);
	s1_NCX_N = (exp((nu_NCX_N .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))) .* (Na_i ^ 3.0) .* c_EC);
	Ka_NCX_N = (1.0 ./ (1.0 + ((Kdact_NCX_N ./ c_i) ^ 2.0)));
	Qcorr_NCX_N = ((mlabfix_T_ - 310.0) ./ 10.0);
	s2_NCX_N = (exp(((nu_NCX_N - 1.0) .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))) .* (Na_EC ^ 3.0) .* c_i);
	KmNa_EC_NCX_N = (87.5 .* 1000.0);
	Kmc_EC_NCX_N = (1.3 .* 1000.0);
	s3_NCX_N = ((Kmc_i_NCX_N .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX_N) ^ 3.0))) + ((KmNa_EC_NCX_N ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX_N))) + (Kmc_EC_NCX_N .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
	I_NCX_N = ( - (3.0 .* (g_NCX_NCX_N .* (Q10NCX_NCX_N ^ Qcorr_NCX_N) .* Ka_NCX_N .* (s1_NCX_N - s2_NCX_N) ./ s3_NCX_N ./ (1.0 + (ksat_NCX_N .* exp(((nu_NCX_N - 1.0) .* Voltage_PM ./ (mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_))))))) .* (1.0 + TTFrac));
	A = (1.0 ./ (1.0 + exp(((Voltage_PM - V_a) ./ A_a))));
	Cl_EC = Cl_EC_init_uM;
	E_Cl =  - (log((Cl_EC ./ Cl_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	I_Cl = ((g_Cl .* (A ^ 4.0) .* (E_Cl - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
	s2_NCX_C = (exp(((nu_NCX_C - 1.0) .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))) .* (Na_EC ^ 3.0) .* c_i);
	Kmc_EC_NCX_C = (1.3 .* 1000.0);
	KmNa_EC_NCX_C = (87.5 .* 1000.0);
	s3_NCX_C = ((Kmc_i_NCX_C .* (Na_EC ^ 3.0) .* (1.0 + ((Na_i ./ KmNa_i_NCX_C) ^ 3.0))) + ((KmNa_EC_NCX_C ^ 3.0) .* c_i .* (1.0 + (c_i ./ Kmc_i_NCX_C))) + (Kmc_EC_NCX_C .* (Na_i ^ 3.0)) + ((Na_i ^ 3.0) .* c_EC) + ((Na_EC ^ 3.0) .* c_i));
	Ka_NCX_C = (1.0 ./ (1.0 + ((Kdact_NCX_C ./ c_i) ^ 2.0)));
	Qcorr_NCX_C = ((mlabfix_T_ - 310.0) ./ 10.0);
	s1_NCX_C = (exp((nu_NCX_C .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))) .* (Na_i ^ 3.0) .* c_EC);
	I_NCX_C = ((2.0 .* (g_NCX_NCX_C .* (Q10NCX_NCX_C ^ Qcorr_NCX_C) .* Ka_NCX_C .* (s1_NCX_C - s2_NCX_C) ./ s3_NCX_C ./ (1.0 + (ksat_NCX_C .* exp(((nu_NCX_C - 1.0) .* Voltage_PM ./ (mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_))))))) .* (1.0 + TTFrac));
	F_PM = ( - (I_DHPR .* Size_PM) - (I_PMCA .* Size_PM) - (I_K_DR .* Size_PM) - (I_K_IR .* Size_PM) - (I_NKX_K .* Size_PM) - (I_Na .* Size_PM) - (I_Cl .* Size_PM) - (I_NCX_C .* Size_PM) - (I_NKX_N .* Size_PM) - (I_NCX_N .* Size_PM) - (I_SOCE .* Size_PM) - (I_CaLeak_PM .* Size_PM));
	L_RyR = (1000.0 ./ 0.002);
	voltProb = (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR))) .* (f_RyR ^  - 2.0))) ^ 4.0) + (L_RyR .* ((1.0 + exp(((Voltage_PM - VBar_RyR) ./ (4.0 .* K_RyR)))) ^ 4.0))));
	alpha_n = (alpha_n0 .* (Voltage_PM - V_n) ./ (1.0 - exp(( - (Voltage_PM - V_n) ./ K_alphan))));
	alpha_m = (alpha_m0 .* (Voltage_PM - V_m) ./ (1.0 - exp(( - (Voltage_PM - V_m) ./ K_alpham))));
	alpha_h = (alpha_h0 .* exp(( - (Voltage_PM - V_h) ./ K_alphah)));
	J_CaLeak_PM = ((I_CaLeak_PM ./ (carrierValence_CaLeak_PM .* mlabfix_F_)) .* unitFactor_CaLeak_PM);
	unitFactor_NKX_N = (8.388606999999999E15 ./ 8388607.0);
	SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
	J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
	tau_hK = (exp(( - (Voltage_PM + 40.0) ./ 25.75)) .* tauUnit);
	h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_hkinf) ./ A_hkinf))));
	J_r5 = ((h_Kinf - h_K) ./ tau_hK);
	beta_n = (beta_n0 .* exp(( - (Voltage_PM - V_n) ./ K_betan)));
	J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
	unitFactor_NKX_K = (8.388606999999999E15 ./ 8388607.0);
	S_inf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_Sinf) ./ A_Sinf))));
	tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_PM + V_tau) ./ 100.0) ^ 2.0))));
	J_r3 = ((S_inf - S) ./ tau_S);
	beta_m = (beta_m0 .* exp(( - (Voltage_PM - V_m) ./ K_betam)));
	J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
	beta_h = (beta_h0 ./ (1.0 + exp(( - (Voltage_PM - V_h) ./ K_betah))));
	J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
	tau_w = (1.0 ./ (100.0 .* (1.0 + c_i)));
	wInf = (1.0 ./ (1.0 + c_i));
	J_r0 = ((wInf - w_RyR) ./ tau_w);
	unitFactor_DHPR = (8.388606999999999E15 ./ 8388607.0);
	J_NKX_N =  - ((I_NKX_N ./ (carrierValence_NKX_N .* mlabfix_F_)) .* unitFactor_NKX_N);
	J_SOCE = ((I_SOCE ./ (carrierValence_SOCE .* mlabfix_F_)) .* unitFactor_SOCE);
	J_NKX_K = ((I_NKX_K ./ (carrierValence_NKX_K .* mlabfix_F_)) .* unitFactor_NKX_K);
	J_K_DR =  - ((I_K_DR ./ (carrierValence_K_DR .* mlabfix_F_)) .* unitFactor_K_DR);
	unitFactor_PMCA = (8.388606999999999E15 ./ 8388607.0);
	unitFactor_NCX_N = (8.388606999999999E15 ./ 8388607.0);
	unitFactor_NCX_C = (8.388606999999999E15 ./ 8388607.0);
	j0_RyR = (300.0 .* volFactor ./ (3.0 .* Size_SRM));
	KFlux_SRM_SR = (Size_SRM ./ Size_SR);
	KFlux_PM_cyto = (Size_PM ./ Size_cyto);
	J_DHPR = ((I_DHPR ./ (carrierValence_DHPR .* mlabfix_F_)) .* unitFactor_DHPR);
	J_NCX_N = ((I_NCX_N ./ (carrierValence_NCX_N .* mlabfix_F_)) .* unitFactor_NCX_N);
	UnitFactor_mV_pF_s_neg_1_pA_neg_1 = (1000.0 ./ 1.0);
	J_NCX_C =  - ((I_NCX_C ./ (carrierValence_NCX_C .* mlabfix_F_)) .* unitFactor_NCX_C);
	KFlux_SRM_cyto = (Size_SRM ./ Size_cyto);
	I_PM =  - device_totalCurrClampElectrode.F;
	J_PMCA =  - ((I_PMCA ./ (carrierValence_PMCA .* mlabfix_F_)) .* unitFactor_PMCA);
	device_PM.Capacitance = (C_PM .* Size_PM);
	unitFactor_Cl = (8.388606999999999E15 ./ 8388607.0);
	device_totalCurrClampElectrode.I = device_totalCurrClampElectrode.F;
	openProb = (voltProb .* (w_RyR ./ wUnit));
	J_RyR = (j0_RyR .* openProb .* (c_SR - c_i));
	J_Cl = ((I_Cl ./ (carrierValence_Cl .* mlabfix_F_)) .* unitFactor_Cl);

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
        J_SERCA = J_SERCA*0.5;
        J_PMCA = J_PMCA*0.5;
        I_PMCA = I_PMCA*0.5;
        J_NKX_K = J_NKX_K*0.5;
        J_NKX_N = J_NKX_N*0.5;
        I_NKX_K = I_NKX_K*0.5;
        I_NKX_N = I_NKX_N*0.5;
    end


	% Rates
    % altered to include f_i and f_SR for calcium buffering (see above)
    % currently set change in K_i, Na_i, and Cl_i to zero - need to fix
    % this so they settle to a reasonable steady state
	dydt = [
		J_r6;    % rate for SOCEProb
		f_SR*( - (KFlux_SRM_SR .* J_RyR) + (KFlux_SRM_SR .* J_SERCA) - (KFlux_SRM_SR .* J_CaLeak_SR));    % rate for c_SR
		J_r5;    % rate for h_K
		J_r0;    % rate for w_RyR
		((UnitFactor_mV_pF_s_neg_1_pA_neg_1 ./ device_PM.Capacitance) .* (I_PM - ( - (I_DHPR .* Size_PM) - (I_PMCA .* Size_PM) - (I_K_DR .* Size_PM) - (I_K_IR .* Size_PM) - (I_NKX_K .* Size_PM) - (I_Na .* Size_PM) - (I_Cl .* Size_PM) - (I_NCX_C .* Size_PM) - (I_NKX_N .* Size_PM) - (I_NCX_N .* Size_PM) - (I_SOCE .* Size_PM) - (I_CaLeak_PM .* Size_PM))));    % rate for Voltage_PM
		0;%((KFlux_PM_cyto .* J_Na) - (KFlux_PM_cyto .* J_NKX_N) + (KFlux_PM_cyto .* J_NCX_N));    % rate for Na_i
		0;%(KFlux_PM_cyto .* J_Cl);    % rate for Cl_i
		f_i*((KFlux_SRM_cyto .* J_RyR) + (KFlux_PM_cyto .* J_DHPR) - (KFlux_PM_cyto .* J_PMCA) - (KFlux_SRM_cyto .* J_SERCA) + (KFlux_SRM_cyto .* J_CaLeak_SR) - (KFlux_PM_cyto .* J_NCX_C) + (KFlux_PM_cyto .* J_SOCE) + (KFlux_PM_cyto .* J_CaLeak_PM));    % rate for c_i
		J_r4;    % rate for n
		J_r2;    % rate for m
		J_r1;    % rate for h
		J_r3;    % rate for S
		0;%( - (KFlux_PM_cyto .* J_K_DR) + (KFlux_PM_cyto .* J_K_IR) + (KFlux_PM_cyto .* J_NKX_K));    % rate for K_i
	];
end