%% Sweep parameters for quick testing of objective function
yinit = [
    0.0122; 	% yinit(1) is the initial condition for 'SOCEProb'
    1500.0;		% yinit(2) is the initial condition for 'c_SR'
    0.9983;		% yinit(3) is the initial condition for 'h_K'
    0.9091;		% yinit(4) is the initial condition for 'w_RyR'
    -88.0;		% yinit(5) is the initial condition for 'Voltage_PM'
    14700.0;	% yinit(6) is the initial condition for 'Na_i'
    5830.0;		% yinit(7) is the initial condition for 'Cl_i'
    0.1;		% yinit(8) is the initial condition for 'c_i'
    0.003;		% yinit(9) is the initial condition for 'n'
    0.0128;		% yinit(10) is the initial condition for 'm'
    0.8051;		% yinit(11) is the initial condition for 'h'
    0.8487;		% yinit(12) is the initial condition for 'S'
    154500.0;	% yinit(13) is the initial condition for 'K_i'
    387;        % yinit(14) is the initial condition for 'CaParv'
    1020;       % yinit(15) is the initial condition for 'MgParv'
    0.3632;     % yinit(16) is the initial consition for 'CATP'
    10.004;     % yinit(17) is the initial condition for 'CaTrop'
    0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    0;	    	% yinit(19) is the initial condition for 'D_2'
    0;	    	% yinit(20) is the initial condition for 'Pre_Pow'
    0;	    	% yinit(21) is the initial condition for 'Post_Pow'
    0;	    	% yinit(22) is the initial condition for 'MgATP'
    8000;       % yinit(23) is the initial condition for 'ATP'
    3000        % yinit(24) is the initial condition for 'p_i_SR'
    0           % yinit(25) is the initial condition for 'PiCa'
    3000        % yinit(26) is the initial condition for 'Pi_Myo'
    ];


param = importdata('InputParam1.xlsx'); % load default parameters
p0 =  param.data;
% load PSO_25-Apr-2024.mat pSol % load best fit parameters from PSO - particle swarm optimization
% pSol(12) = pSol(12)*0.2;
% pSol(46:96) = 1;
% pPSO = pSol .* p0';
pSol = [1.33684210526316	0.847368421052632	1.05263157894737	1.66842105263158	1.36842105263158	0.547368421052632	0.705263157894737	1.33684210526316	0.547368421052632	1.49473684210526	1.33684210526316	1.98421052631579	0.500000000000000	0.563157894736842	1.02105263157895	1.19473684210526	1.02105263157895	0.578947368421053	0.642105263157895	1.43157894736842	0.610526315789474	1.51052631578947	1.87368421052632	1.06842105263158	1.51052631578947	0.500000000000000	1.43157894736842	0.500000000000000	1.58947368421053	1.30526315789474	0.594736842105263	1.36842105263158	1.00526315789474	1.51052631578947	0.736842105263158	0.689473684210526	1.08421052631579	0.784210526315790	1.16315789473684	0.563157894736842	1.08421052631579	1.87368421052632	1.85789473684211	0.894736842105263	1.74736842105263	1.55789473684211	0.926315789473684	0.926315789473684	0.768421052631579	1.27368421052632	1.10000000000000	1.68421052631579	1.14736842105263	1.28947368421053	1.81052631578947	0.878947368421053	1.08421052631579	1.54210526315789	1.22631578947368	1.25789473684211	1.87368421052632	1.79473684210526	1.76315789473684	1.33684210526316	1.10000000000000	1.40000000000000	1.70000000000000	1.58947368421053	0.515789473684211	1.33684210526316	1.88947368421053	1.05263157894737	1.17894736842105	0.926315789473684	0.878947368421053	1.36842105263158	1.25789473684211	1.22631578947368	1.70000000000000	0.784210526315790	0.942105263157895	1.62105263157895	1.76315789473684	0.642105263157895	0.563157894736842	1.51052631578947	1.66842105263158	1.77894736842105	1.87368421052632	0.610526315789474	0.815789473684211	0.831578947368421	1.58947368421053	1.52631578947368	1.74736842105263	1.87368421052632];
pSol(95:96) = 1;
pPSO = pSol' .* p0;


[TimeSS,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pPSO, tic, 2); % compute steady state solution
tSol = 0:.0001:10;
freq = 60; 
% ySS(end,26) = 0;
% ySS(end,28) = 2500;
[Time,Y] = SkelMuscleCa_dydt(tSol, freq, 0, ySS(end,:), pPSO, tic, 1); % compute time-dependent solution
% plot calcium (the 8th variable)


% figure
% scatter(Time, max(Y(:,8)))
% title('Time vs max Ca2+')
% xlabel('Time (seconds)');
% ylabel('[Ca2+] (µM)'); 
 
figure
plot(Time, Y(:,8))
title('Time vs Ca2+')
xlabel('Time (seconds)');
ylabel('[Ca2+] (µM)'); 

figure
plot(Time, Y(:,2))
title('Time vs Ca2+ SR')
xlabel('Time (seconds)');
ylabel('[Ca2+]_{SR} (µM)'); 
% % % 
figure
plot(Time, Y(:,23))
title('Time vs ATP, freq=60 Hz')
xlabel('Time (seconds)');
ylabel('[ATP] (µM)'); 
% 
figure
plot(Time, Y(:,21))
title('Time vs Force, freq=60 Hz')
xlabel('Time (seconds)');
ylabel('Post Power Stroke'); 

figure
plot(Time, Y(:,5))
title('Time vs SL Voltage ')
xlabel('Time (seconds)');
ylabel('SL Voltage'); 
% % 
% figure
% plot(Time, Y(:,6))
% title('Time vs Sodium, k_{onParvMg} = 0.033*(1/1.5)')
% xlabel('Time (seconds)');
% ylabel('Sodium'); 
% figure
% plot(Time, Y(:,7))
% title('Time vs Chlorine, k_{onParvMg} = 0.033*(1/1.5)')
% xlabel('Time (seconds)');
% ylabel('SL Chlorine'); 
% figure
% plot(Time, Y(:,13))
% title('Time vs Potassium, k_{onParvMg} = 0.033*(1/1.5)')
% xlabel('Time (seconds)');
% ylabel('Potassium'); 
% 
% 
% figure
% plot(Time, Y(:,26))
% title('26s')
% % legend('Calcium ion Conc');  
% xlabel('Time (seconds)');
% ylabel('[Ca2+] (µM)'); 
% 
% figure
% plot(Time, Y(:,27))
% title('time vs 27, freq=100 Hz')
% legend('CaTrop Conc');  
% xlabel('Time (seconds)');
% ylabel('[CaTrop] (µM)'); 
% % 
% figure
% plot(Time, Y(:,28))
figure
subplot(3,1,1)
plot(Time, Y(:,24))
title('time vs p_i_{SR}, freq=60 Hz')
legend('p_i_{SR}');  
xlabel('Time (seconds)');
ylabel('[p_i_{SR}] (µM)'); 
subplot(3,1,2)
plot(Time, Y(:,25))
title('time vs PiCa, freq=60 Hz')
legend('PiCa Conc');  
xlabel('Time (seconds)');
ylabel('[PiCa] (µM)'); 
subplot(3,1,3)
plot(Time, Y(:,26))
title('time vs Pi_{Myo}, freq=60 Hz')
legend('Pi_{Myo} Conc');  
xlabel('Time (seconds)');
ylabel('[Pi_{Myo}] (µM)'); 

% 
% figure
% plot(Time, Y(:,25))
% title('time vs ATP, k_{onParvMg} = 0.033*(1/1.5)')
% legend('ATP Conc');  
% xlabel('Time (seconds)');
% ylabel('[ATP] (µM)'); 
% 
% figure
% plot(Time, Y(:,23))
% title('time vs A2, k_{onParvMg} = 0.033*(1/1.5)')
% legend('ATP Conc');  
% xlabel('Time (seconds)');
% ylabel('[A2] (µM)'); 
%% 
% figure
% subplot(3,1,1)
% plot(Time, Y(:,8))
% title('time vs Ca2+, freq=100 Hz')
% legend('Calcium ion Conc');  
% xlabel('Time (seconds)');
% ylabel('[Ca2+] (µM)'); 
% subplot(3,1,2)
% plot(Time, Y(:,17))
% title('time vs CaTrop, freq=100 Hz')
% legend('CaTrop Conc');  
% xlabel('Time (seconds)');
% ylabel('[CaTrop] (µM)'); 
% subplot(3,1,3)
% plot(Time, Y(:,18))
% title('time vs CaCaTrop, freq=100 Hz')
% legend('CaTrop Conc');  
% xlabel('Time (seconds)');
% ylabel('[CaTrop] (µM)'); 

% % 
% figure
% subplot(2,3,1)
% plot(Time, Y(:,17))
% title('time vs dCT')
% subplot(2,3,2)
% plot(Time, Y(:,18))
% title('time vs dCCT')
% subplot(2,3,3)
% plot(Time, Y(:,19))
% title('time vs dD0')
% subplot(2,3,4)
% plot(Time, Y(:,19))
% title('time vs dD1')
% subplot(2,3,5)
% plot(Time, Y(:,22))
% title('time vs dPrePower')
% subplot(2,3,6)
% plot(Time, Y(:,23))
% title('time vs dPostPower')

%% 
% 
% figure
% plot(Time, Y(:,8))

%Adding Senneff Corssbridge Cycling
Sol= {Time*1000, Y(:,8)}; %to convert Anusha's time (s) to (ms)
y0=[0 ; 0; 0; 0 ; 0 ];
kTon = 0.0885; % /0.04425 
kToff = 0.115; 
k0on = 0;
k0off =  0.15; 
kCaOn = 0.15;
kCaOff = 0.05; 
h0 = 0.24; % 0.08
hP =  0.18; % 0.06 
f0 = 1.5; %0.5
fP =  15; % /5 
g0 = 0.12; % /0.04
Ca = 1e60; 
%
kCATPon= 0.15;
kCATPoff= 0.3;
kMATPon = 0.0015;
kMATPoff= 0.15;
kH = 1000;
kHYD = 0.1;
%
TropTotal = 140;
paramVec = [kTon, kToff, k0on, k0off, kCaOn, kCaOff, h0, hP, f0, fP, g0, Ca, TropTotal, kCATPon, kCATPoff, kMATPon, kMATPoff, kH, kHYD];

%Adding Senneff ATP Hydrolysis

tspan= [0 1000]; 
odeoptions = odeset();
[t,y]= ode15s(@crossbridge_model, tspan, y0, odeoptions, paramVec, Sol);

figure
plot(t,y(:,5))
legend('A2');  
xlabel('Time');
ylabel('[A2]'); 
figure
plot(t,y)
legend('CaT', 'CaCaT',  'D2', 'A1', 'A2');  


function dydt=crossbridge_model(t,y,paramVec,Sol)
    dydt= zeros(5,1);  
    kTon = paramVec(1); % /0.04425 
    kToff = paramVec(2); 
    k0on = paramVec(3);
    k0off =  paramVec(4); 
    kCaOn = paramVec(5);
    kCaOff = paramVec(6); 
    h0 = paramVec(7);  
    hP =  paramVec(8);  
    f0 = paramVec(9); 
    fP =  paramVec(10);  
    g0 = paramVec(11);  
    Ca = paramVec(12);
    T0 = paramVec(13) - y(2)  -y(3) - y(4) - y(5) - y(1); 
    % Sol= paramVec(14); 
    % kCATPon = paramVec(14);
    % kCATPoff = paramVec(15);
    % kMATPon = paramVec(16); &&&&&&&&&& 
    % kMATPoff = paramVec(17);
    kH = paramVec(18);
    kHYD = paramVec(19); 
    %JhydMyo =  kHYD *( ATP_myo / (kH + ATP_myo) );

    % TimeFinal= Sol{1}; 
    % CaSol = Sol{2};
    % Ca= interp1(TimeFinal, CaSol,t); 

    dydt(1)= (kTon * Ca * T0) - ( kToff* y(1)) - (kTon*Ca *y(1)) + (kToff* y(2));% - (k0on * y(1));% + (k0off * y(4));
    dydt(2)= (kTon * Ca * y(1)) - ( kToff* y(2)) - (kCaOn * y(2)) + (kCaOff* y(5));
    % dydt(3)=  -(kTon * Ca * y(3)) + ( kToff* y(4)) + (k0on * T0) - (k0off * y(3));
    % dydt(4)=  (kTon * Ca  * y(3)) - ( kToff* y(4)) + (k0on * y(1)) - (k0off * y(4)) - (kTon *Ca*y(4)) + (kToff* y(5));
    % dydt(5)=  (kTon * Ca * y(4)) - ( kToff* y(5)) + (kCaOn * y(2)) - (kCaOff * y(5)) - (f0 * y(5)) + (fP * y(6)) + (g0* y(7)/700);
    dydt(3)=  - ( kToff* y(3)) + (kCaOn * y(2)) - (kCaOff * y(3))  + (fP * y(4)) + (g0* y(5)/700);
    dydt(4)=   (f0* y(3)) - ( fP * y(4)) + (hP * y(5)) - (h0 * y(4));
    dydt(5)=  -(hP * y(5)) +  (h0 * y(4)) - (g0 * y(5)/700);
    %atp: 
    % dydt(8)= -JhydM - (kCATPon* Ca * ATPm - kCATPoff*CaATP) - (kMATPon* Mg * ATPm - kMATPoff*MgATPm ) + tATP*(ATPtm - ATPm)/Vm;
    % 
end
  