function [T,Y,SSParam] = SkelMuscleCa(argTimeSpan, freq, lowATP, selParam)
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
param = [
	288.0;		% param(1) is 'alpha_m0'
	1.0;		% param(2) is 'carrierValence_K_IR'
	37698.0;		% param(3) is 'Size_SRM'
	0.17;		% param(4) is 'K_PMCA'
	1.0;		% param(5) is 'wUnit_RyR'
	-20.0;		% param(6) is 'VBar_RyR'
	0.4;		% param(7) is 'delta'
	13000.0;		% param(8) is 'K_mNa_NKX_N'
	0.2;		% param(9) is 'f_RyR'
	13000.0;		% param(10) is 'K_mNa_NKX_K'
	90.0;		% param(11) is 'V_tau'
	.1e-6;%6.21E-6;		% param(12) is 'J_NaK_NKX_N'
	0.0;		% param(13) is 'Voltage_SRM'
	.1e-6;%6.21E-6;		% param(14) is 'J_NaK_NKX_K'
	0.01;		% param(15) is 'tau_SOCEProb'
	0.9983;		% param(16) is 'h_K_init_molecules_um_2'
	14700.0;		% param(17) is 'Na_i_init_uM'
	7.0;		% param(18) is 'K_alphan'
	10.0;		% param(19) is 'K_alpham'
	14.7;		% param(20) is 'K_alphah'
	500.0;		% param(21) is 'c_ref'
	-1.0;		% param(22) is 'carrierValence_Cl'
	1.0;		% param(23) is 'SUnit'
	7.5;		% param(24) is 'A_hkinf'
	1.0;		% param(25) is 'SOCEProbUnit'
	-40.0;		% param(26) is 'V_n'
	-46.0;		% param(27) is 'V_m'
	0.1965;		% param(28) is 'g_Cl'
	4380.0;		% param(29) is 'beta_h0'
	1.0;		% param(30) is 'carrierValence_NCX_N'
	-45.0;		% param(31) is 'V_h'
	119381.0;		% param(32) is 'Size_cyto'
	70.0;		% param(33) is 'V_a'
	2.0;		% param(34) is 'carrierValence_NCX_C'
	128000.0;		% param(35) is 'Cl_EC_init_uM'
	1.0;		% param(36) is 'hUnit'
	0.129;		% param(37) is 'g_NCX_NCX_N'
	0.1;		% param(38) is 'c_i_init_uM'
	0.129;		% param(39) is 'g_NCX_NCX_C'
	0.8051;		% param(40) is 'h_init_molecules_um_2'
	147000.0;		% param(41) is 'Na_EC_init_uM'
	1.0;		% param(42) is 'carrierValence_Na'
	5.8;		% param(43) is 'A_Sinf'
	1300.0;		% param(44) is 'c_EC_init_uM'
	8.04;		% param(45) is 'g_Na'
	6283.0;		% param(46) is 'Size_SR'
	6.59;		% param(47) is 'Kmc_i_NCX_N'
	8.1;		% param(48) is 'alpha_h0'
	1.0;		% param(49) is 'netValence_r7'
	2.0;		% param(50) is 'carrierValence_SOCE'
	1.0;		% param(51) is 'netValence_r6'
	1.0;		% param(52) is 'netValence_r5'
	1.0;		% param(53) is 'carrierValence_K_DR'
	1.0;		% param(54) is 'netValence_r4'
	6.59;		% param(55) is 'Kmc_i_NCX_C'
	1.0;		% param(56) is 'netValence_r3'
	1.0;		% param(57) is 'netValence_r2'
	1.0;		% param(58) is 'netValence_r1'
	1.0;		% param(59) is 'netValence_r0'
	154500.0;		% param(60) is 'K_i_init_uM'
	0.01;		% param(61) is 'C_PM'
	0.0;		% param(62) is 'Voltage_SRM_init'
	1000000.0;		% param(63) is 'Size_EC'
	9.64853321E-5;		% param(64) is 'mlabfix_F_nmol_'
	1000.0;		% param(65) is 'K_mK_NKX_N'
	300.0;		% param(66) is 'mlabfix_T_'
	1000.0;		% param(67) is 'K_millivolts_per_volt'
	1000.0;		% param(68) is 'K_mK_NKX_K'
	1.0;		% param(69) is 'h_KUnit'
	1.0;		% param(70) is 'nUnit'
	1.0;		% param(71) is 'C_SRM'
	1.0;		% param(72) is 'carrierValence_CaLeak_SR'
	1.57;		% param(73) is 'Q10NCX_NCX_N'
	1.0;		% param(74) is 'carrierValence_SERCA'
	5830.0;		% param(75) is 'Cl_i_init_uM'
	12566.0;		% param(76) is 'Size_PM'
	1.57;		% param(77) is 'Q10NCX_NCX_C'
	2.0;		% param(78) is 'carrierValence_DHPR'
	1500.0;		% param(79) is 'c_SR_init_uM'
	0.648;		% param(80) is 'g_K'
	4.5;		% param(81) is 'K_RyR'
	-20;		% param(82) is 'VBar_DHPR'
	3.141592653589793;		% param(83) is 'mlabfix_PI_'
	-78.0;		% param(84) is 'V_Sinf'
	96485.3321;		% param(85) is 'mlabfix_F_'
	0.0;		% param(86) is 'w_DHPR_init_molecules_um_2'
	1.0;		% param(87) is 'wUnit_DHPR'
	8314.46261815;		% param(88) is 'mlabfix_R_'
	0.8487;		% param(89) is 'S_init_molecules_um_2'
	67.0;		% param(90) is 'beta_n0'
	0.05;		% param(91) is 'Kdact_NCX_N'
	10000.0;		% param(92) is 'S_i'
	0.05;		% param(93) is 'Kdact_NCX_C'
	1.0E-9;		% param(94) is 'mlabfix_K_GHK_'
	0.27;		% param(95) is 'nu_NCX_N'
	1.0;		% param(96) is 'carrierValence_RyR'
	-40.0;		% param(97) is 'V_hkinf'
	2.0;		% param(98) is 'carrierValence_PMCA'
	4000.0;		% param(99) is 'K_EC_init_uM'
	0.27;		% param(100) is 'nu_NCX_C'
	40.0;		% param(101) is 'K_betan'
	18.0;		% param(102) is 'K_betam'
	9.0;		% param(103) is 'K_betah'
	1380.0;		% param(104) is 'beta_m0'
	6.02214179E11;		% param(105) is 'mlabfix_N_pmol_'
	2.0;		% param(106) is 'carrierValence_CaLeak_PM'
	1.0;		% param(107) is 'mUnit'
	0.003;		% param(108) is 'n_init_molecules_um_2'
	0.2;		% param(109) is 'f_DHPR'
	4.5;		% param(110) is 'K_DHPR'
	1000000.0;		% param(111) is 'K_S'
	-20000.0;		% param(112) is 'device_totalCurrClampElectrode.F'
	9.5E8;		% param(113) is 'K_K'
	0.0128;		% param(114) is 'm_init_molecules_um_2'
	13.1;		% param(115) is 'alpha_n0'
	0.32;		% param(116) is 'ksat_NCX_N'
	-88.0;		% param(117) is 'Voltage_PM_init'
	3.75;		% param(118) is 'TTFrac'
	0.001660538783162726;		% param(119) is 'KMOLE'
	150.0;		% param(120) is 'A_a'
	0.111;		% param(121) is 'G_K'
	0.9091;		% param(122) is 'w_RyR_init_molecules_um_2'
	4.743E-6;		% param(123) is 'SOCEProb_init_molecules_um_2'
	1.0;		% param(124) is 'K_SERCA'
	0.32;		% param(125) is 'ksat_NCX_C'
	1.0;		% param(126) is 'carrierValence_NKX_N'
	1.0;		% param(127) is 'tauUnit'
	1.0;		% param(128) is 'carrierValence_NKX_K'
];

param(45) = selParam*param(45);


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
	alpha_m0 = p(1);
	carrierValence_K_IR = p(2);
	Size_SRM = p(3);
	K_PMCA = p(4);
	wUnit_RyR = p(5);
	VBar_RyR = p(6);
	delta = p(7);
	K_mNa_NKX_N = p(8);
	f_RyR = p(9);
	K_mNa_NKX_K = p(10);
	V_tau = p(11);
	J_NaK_NKX_N = p(12);
	Voltage_SRM = p(13);
	J_NaK_NKX_K = p(14);
	tau_SOCEProb = p(15);
	h_K_init_molecules_um_2 = p(16);
	Na_i_init_uM = p(17);
	K_alphan = p(18);
	K_alpham = p(19);
	K_alphah = p(20);
	c_ref = p(21);
	carrierValence_Cl = p(22);
	SUnit = p(23);
	A_hkinf = p(24);
	SOCEProbUnit = p(25);
	V_n = p(26);
	V_m = p(27);
	g_Cl = p(28);
	beta_h0 = p(29);
	carrierValence_NCX_N = p(30);
	V_h = p(31);
	Size_cyto = p(32);
	V_a = p(33);
	carrierValence_NCX_C = p(34);
	Cl_EC_init_uM = p(35);
	hUnit = p(36);
	g_NCX_NCX_N = p(37);
	c_i_init_uM = p(38);
	g_NCX_NCX_C = p(39);
	h_init_molecules_um_2 = p(40);
	Na_EC_init_uM = p(41);
	carrierValence_Na = p(42);
	A_Sinf = p(43);
	c_EC_init_uM = p(44);
	g_Na = p(45);
	Size_SR = p(46);
	Kmc_i_NCX_N = p(47);
	alpha_h0 = p(48);
	netValence_r7 = p(49);
	carrierValence_SOCE = p(50);
	netValence_r6 = p(51);
	netValence_r5 = p(52);
	carrierValence_K_DR = p(53);
	netValence_r4 = p(54);
	Kmc_i_NCX_C = p(55);
	netValence_r3 = p(56);
	netValence_r2 = p(57);
	netValence_r1 = p(58);
	netValence_r0 = p(59);
	K_i_init_uM = p(60);
	C_PM = p(61);
	Voltage_SRM_init = p(62);
	Size_EC = p(63);
	mlabfix_F_nmol_ = p(64);
	K_mK_NKX_N = p(65);
	mlabfix_T_ = p(66);
	K_millivolts_per_volt = p(67);
	K_mK_NKX_K = p(68);
	h_KUnit = p(69);
	nUnit = p(70);
	C_SRM = p(71);
	carrierValence_CaLeak_SR = p(72);
	Q10NCX_NCX_N = p(73);
	carrierValence_SERCA = p(74);
	Cl_i_init_uM = p(75);
	Size_PM = p(76);
	Q10NCX_NCX_C = p(77);
	carrierValence_DHPR = p(78);
	c_SR_init_uM = p(79);
	g_K = p(80);
	K_RyR = p(81);
	VBar_DHPR = p(82);
	mlabfix_PI_ = p(83);
	V_Sinf = p(84);
	mlabfix_F_ = p(85);
	w_DHPR_init_molecules_um_2 = p(86);
	wUnit_DHPR = p(87);
	mlabfix_R_ = p(88);
	S_init_molecules_um_2 = p(89);
	beta_n0 = p(90);
	Kdact_NCX_N = p(91);
	S_i = p(92);
	Kdact_NCX_C = p(93);
	mlabfix_K_GHK_ = p(94);
	nu_NCX_N = p(95);
	carrierValence_RyR = p(96);
	V_hkinf = p(97);
	carrierValence_PMCA = p(98);
	K_EC_init_uM = p(99);
	nu_NCX_C = p(100);
	K_betan = p(101);
	K_betam = p(102);
	K_betah = p(103);
	beta_m0 = p(104);
	mlabfix_N_pmol_ = p(105);
	carrierValence_CaLeak_PM = p(106);
	mUnit = p(107);
	n_init_molecules_um_2 = p(108);
	f_DHPR = p(109);
	K_DHPR = p(110);
	K_S = p(111);
	device_totalCurrClampElectrode.F = p(112);
	K_K = p(113);
	m_init_molecules_um_2 = p(114);
	alpha_n0 = p(115);
	ksat_NCX_N = p(116);
	Voltage_PM_init = p(117);
	TTFrac = p(118);
	KMOLE = p(119);
	A_a = p(120);
	G_K = p(121);
	w_RyR_init_molecules_um_2 = p(122);
	SOCEProb_init_molecules_um_2 = p(123);
	K_SERCA = p(124);
	ksat_NCX_C = p(125);
	carrierValence_NKX_N = p(126);
	tauUnit = p(127);
	carrierValence_NKX_K = p(128);
	% Functions
	UnitFactor_uM_um3_molecules_neg_1 = (1.0 .* (KMOLE ^ 1.0));
	K_EC = K_EC_init_uM;
	E_K_K_IR = (log((K_EC ./ K_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	K_R = (K_EC .* exp(( - delta .* E_K_K_IR .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))));
	I_K_IR = ((G_K .* ((K_R ^ 2.0) ./ (K_K + (K_R ^ 2.0))) .* (1.0 - ((1.0 + ((K_S ./ ((S_i ^ 2.0) .* exp(((2.0 .* (1.0 - delta) .* Voltage_PM .* mlabfix_F_) ./ (mlabfix_R_ .* mlabfix_T_))))) .* (1.0 + ((K_R ^ 2.0) ./ K_K)))) ^  - 1.0)) .* (E_K_K_IR - Voltage_PM)) .* (1.0 + TTFrac));
	unitFactor_K_IR = (1.0E9 ./ 1.0);
	J_K_IR = ((I_K_IR ./ (carrierValence_K_IR .* mlabfix_F_)) .* unitFactor_K_IR);
	Na_EC = Na_EC_init_uM;
	E_Na = (log((Na_EC ./ Na_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	unitFactor_CaLeak_PM = (1.0E9 ./ 1.0);
	KmNa_i_NCX_N = (12.29 .* 1000.0);
	KmNa_i_NCX_C = (12.29 .* 1000.0);
	volFactor = (Size_cyto ./ (3.1416 .* 0.26));
	nu_SERCA = (4875.0 .* volFactor ./ Size_SRM);
	LumpedJ_SERCA = (602.214179 .* Size_SRM .* nu_SERCA .* c_i ./ (K_SERCA + c_i));
	I_Na = ((g_Na .* ((m ./ mUnit) ^ 3.0) .* (h ./ hUnit) .* (S ./ SUnit) .* (E_Na - Voltage_PM)) .* (1.0 + (0.1 .* TTFrac)));
	unitFactor_Na = (1.0E9 ./ 1.0);
	J_Na = ((I_Na ./ (carrierValence_Na .* mlabfix_F_)) .* unitFactor_Na);
	unitFactor_SOCE = (1.0E9 ./ 1.0);
	unitFactor_K_DR = (1.0E9 ./ 1.0);
	fPump_NKX_K = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))))) ^  - 1.0);
	I_NKX_K = ((2.0 .* (fPump_NKX_K .* mlabfix_F_ .* J_NaK_NKX_K ./ (((1.0 + (K_mK_NKX_K ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX_K ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_PMCA = (5.37 ./ Size_PM);
	I_PMCA = ( - (g_PMCA .* c_i ./ (K_PMCA + c_i)) .* (1.0 + TTFrac));
	fPump_NKX_N = ((1.0 + (0.12 .* exp(( - 0.1 .* Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_)))) + ((0.04 ./ 7.0) .* (exp((Na_EC ./ 67300.0)) - 1.0) .* exp(( - Voltage_PM .* mlabfix_F_ ./ (mlabfix_R_ .* mlabfix_T_))))) ^  - 1.0);
	I_NKX_N = ( - (3.0 .* (fPump_NKX_N .* mlabfix_F_ .* J_NaK_NKX_N ./ (((1.0 + (K_mK_NKX_N ./ K_EC)) ^ 2.0) .* ((1.0 + (K_mNa_NKX_N ./ Na_i)) ^ 3.0)))) .* (1.0 + (0.1 .* TTFrac)));
	g_SOCE = (0.01 ./ 210.44);
	c_EC = c_EC_init_uM;
	E_Ca = (log((c_EC ./ c_i)) .* (mlabfix_R_ .* mlabfix_T_) ./ (2.0 .* mlabfix_F_));
	I_SOCE = ((g_SOCE .* (E_Ca - Voltage_PM) .* (SOCEProb ./ SOCEProbUnit)) .* (1.0 + TTFrac));
	E_K_K_DR = (log((K_EC ./ K_i)) .* mlabfix_R_ .* mlabfix_T_ ./ mlabfix_F_);
	I_K_DR = ((g_K .* ((n ./ nUnit) ^ 4.0) .* (h_K ./ h_KUnit) .* (E_K_K_DR - Voltage_PM)) .* (1.0 + (0.45 .* TTFrac)));
	g_PMLeak = 2.1362E-6 ; %(5.0 .* 2.0E-6);
	I_CaLeak_PM = ((g_PMLeak .* (E_Ca - Voltage_PM)) .* (1.0 + TTFrac));
	L_DHPR = (1000.0 ./ 0.002);
	openDHPR = ((((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) ./ (((1.0 + (exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR))) .* (f_DHPR ^  - 2.0))) ^ 4.0) + (L_DHPR .* ((1.0 + exp(((Voltage_PM - VBar_DHPR) ./ (4.0 .* K_DHPR)))) ^ 4.0)))) .* (w_DHPR ./ wUnit_DHPR));
	g0_DHPR = (9.39 ./ 100.0);
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
	J_CaLeak_PM = ((I_CaLeak_PM ./ (carrierValence_CaLeak_PM .* mlabfix_F_)) .* unitFactor_CaLeak_PM);
	wInf_r7 = (1.0 ./ (1.0 + c_i));
	J_r7 = ((wInf_r7 - w_DHPR) ./ tau_w_r7);
	unitFactor_NKX_N = (1.0E9 ./ 1.0);
	SOCEProb_inf = (1.0 ./ (1.0 + ((c_SR ./ c_ref) ^ 4.0)));
	J_r6 = ((SOCEProb_inf - SOCEProb) ./ tau_SOCEProb);
	tau_hK = (exp(( - (Voltage_PM + 40.0) ./ 25.75)) .* tauUnit);
	h_Kinf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_hkinf) ./ A_hkinf))));
	J_r5 = ((h_Kinf - h_K) ./ tau_hK);
	beta_n = (beta_n0 .* exp(( - (Voltage_PM - V_n) ./ K_betan)));
	J_r4 = ((alpha_n .* (1.0 - n)) - (beta_n .* n));
	unitFactor_NKX_K = (1.0E9 ./ 1.0);
	S_inf = (1.0 ./ (1.0 + exp(((Voltage_PM - V_Sinf) ./ A_Sinf))));
	tau_S = (60.0 ./ (0.2 + (5.65 .* (((Voltage_PM + V_tau) ./ 100.0) ^ 2.0))));
	J_r3 = ((S_inf - S) ./ tau_S);
	beta_m = (beta_m0 .* exp(( - (Voltage_PM - V_m) ./ K_betam)));
	J_r2 = ((alpha_m .* (1.0 - m)) - (beta_m .* m));
	beta_h = (beta_h0 ./ (1.0 + exp(( - (Voltage_PM - V_h) ./ K_betah))));
	J_r1 = ((alpha_h .* (1.0 - h)) - (beta_h .* h));
	wInf_r0 = (1.0 ./ (1.0 + c_i));
	J_r0 = ((wInf_r0 - w_RyR) ./ tau_w_r0);
	unitFactor_DHPR = (1.0E9 ./ 1.0);
	J_NKX_N =  - ((I_NKX_N ./ (carrierValence_NKX_N .* mlabfix_F_)) .* unitFactor_NKX_N);
	J_SOCE = ((I_SOCE ./ (carrierValence_SOCE .* mlabfix_F_)) .* unitFactor_SOCE);
	J_NKX_K = ((I_NKX_K ./ (carrierValence_NKX_K .* mlabfix_F_)) .* unitFactor_NKX_K);
	J_K_DR =  - ((I_K_DR ./ (carrierValence_K_DR .* mlabfix_F_)) .* unitFactor_K_DR);
	unitFactor_RyR = (1.0 ./ 6.02214179E11);
	LumpedI_RyR =  - ((carrierValence_RyR .* mlabfix_F_ .* LumpedJ_RyR) .* unitFactor_RyR);
	unitFactor_PMCA = (1.0E9 ./ 1.0);
	unitFactor_NCX_N = (1.0E9 ./ 1.0);
	unitFactor_NCX_C = (1.0E9 ./ 1.0);
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
	unitFactor_Cl = (1.0E9 ./ 1.0);
	device_totalCurrClampElectrode.I = device_totalCurrClampElectrode.F;
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
    J_NaK_SS = J_NaK_NKX_K * ((J_K_DR - J_K_IR) / J_NKX_K);
    J_NKX_K_SS = J_NKX_K * J_NaK_SS / J_NaK_NKX_K;
    g_NCX_SS = g_NCX_NCX_N*((J_NKX_N * (J_NaK_SS / J_NaK_NKX_N)) - J_Na)/ J_NCX_N;
    CaFluxDeficit = (J_NCX_C*g_NCX_SS/g_NCX_NCX_C) + J_PMCA - J_SOCE - J_DHPR;
    g_PMleak_SS = g_PMLeak*(CaFluxDeficit/J_CaLeak_PM);
    nu_leakSR_SS = nu_leakSR*(LumpedJ_SERCA - LumpedJ_RyR)/LumpedJ_CaLeak_SR;


    SSParam = [J_NaK_SS,g_NCX_SS,g_PMleak_SS,nu_leakSR_SS]; % J_NaK, etc.

    I_Cl = 0;

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