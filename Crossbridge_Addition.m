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
    0.3632;     % yinit(16) is the inital consition for 'CATP'
    10.004;     % yinit(17) is the initial condition for 'CaTrop'
    0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    0;	    	% yinit(19) is the initial condition for 'D_0'
    0;	    	% yinit(20) is the initial condition for 'D_1'
    0;	    	% yinit(21) is the initial condition for 'D_2'
    0;	    	% yinit(22) is the initial condition for 'Pre_Pow'
    0;	    	% yinit(23) is the initial condition for 'Post_Pow'
    0;	    	% yinit(24) is the initial condition for 'MgATP'
    8000;       % yinit(25) is the initial condition for 'ATP'
    3000          % yinit(26) is the initial condition for 'p_i_SR' 
    0       % yinit(27) is the initial condition for 'PiCa'
    3000       % yinit(28) is the initial condition for 'Pi_Myo'
    ];

param = importdata('InputParam1.xlsx'); % load default parameters
p0 =  param.data;
load PSO_25-Apr-2024.mat pSol % load best fit parameters from PSO - particle swarm optimization
% pSol(12) = pSol(12)*0.2;
pSol(46:75) = 1;
pPSO = pSol.*p0';
[TimeSS,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pPSO, tic, 2); % compute steady state solution
tSol = 0:.0001:10;
freq = 60; 
[Time,Y] = SkelMuscleCa_dydt(tSol, freq, 0, ySS(end,:), pPSO, tic, 1); % compute time-dependent solution
% plot calcium (the 8th variable)

figure
plot(Time, Y(:,8))
title('time vs Ca2+, freq=100 Hz Jprod=100*kHyd ; atp depend on pumps')
% legend('Calcium ion Conc');  
% xlabel('Time (seconds)');
% ylabel('[Ca2+] (µM)'); 

figure
plot(Time, Y(:,25))
title('Time vs ATP, freq=100 Hz Jprod=100*kHyd; atp depend on pumps')
% legend('Calcium ion Conc');  
% xlabel('Time (seconds)');
% ylabel('[Ca2+] (µM)'); 
% 
% figure
% plot(Time, Y(:,26))
% title('26s')
% % % legend('Calcium ion Conc');  
% % xlabel('Time (seconds)');
% % ylabel('[Ca2+] (µM)'); 
% 
% figure
% plot(Time, Y(:,27))
% title('time vs 27, freq=100 Hz')
% legend('CaTrop Conc');  
% xlabel('Time (seconds)');
% ylabel('[CaTrop] (µM)'); 
% 
% figure
% plot(Time, Y(:,28))
figure
subplot(3,1,1)
plot(Time, Y(:,26))
title('time vs p_i_SR, freq=100 Hz')
legend('p_i_SR');  
xlabel('Time (seconds)');
ylabel('[p_i_SR] (µM)'); 
subplot(3,1,2)
plot(Time, Y(:,27))
title('time vs PiCa, freq=100 Hz')
legend('PiCa Conc');  
xlabel('Time (seconds)');
ylabel('[PiCa] (µM)'); 
subplot(3,1,3)
plot(Time, Y(:,28))
title('time vs Pi_Myo, freq=100 Hz')
legend('Pi_Myo Conc');  
xlabel('Time (seconds)');
ylabel('[Pi_Myo] (µM)'); 



% figure
% plot(Time, Y(:,25))
% title('time vs ATP')
% legend('ATP Conc');  
% xlabel('Time (seconds)');
% ylabel('[ATP] (µM)'); 
% 
% figure
% plot(Time, Y(:,23))
% title('time vs A2')
% legend('ATP Conc');  
% xlabel('Time (seconds)');
% ylabel('[A2] (µM)'); 
% % 
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

% 
figure
subplot(2,3,1)
plot(Time, Y(:,17))
title('time vs dCT')
subplot(2,3,2)
plot(Time, Y(:,18))
title('time vs dCCT')
subplot(2,3,3)
plot(Time, Y(:,19))
title('time vs dD0')
subplot(2,3,4)
plot(Time, Y(:,19))
title('time vs dD1')
subplot(2,3,5)
plot(Time, Y(:,22))
title('time vs dPrePower')
subplot(2,3,6)
plot(Time, Y(:,23))
title('time vs dPostPower')

%% 
 
figure
plot(Time, Y(:,8))

%Adding Senneff Corssbridge Cycling
Sol= {Time*1000, Y(:,8)}; %to convert Anusha's time (s) to (ms)
y0=[0 ; 0; 0; 0 ; 0 ; 0 ; 0 ];
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
Ca = 100 ; 
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
plot(t,y(:,7))
legend('A2');  
xlabel('Time');
ylabel('[A2]'); 
figure
plot(t,y)
legend('CaT', 'CaCaT', 'D0', 'D1', 'D2', 'A1', 'A2', 'MgATP', 'ATP');  


function dydt=crossbridge_model(t,y,paramVec,Sol)
    dydt= zeros(7,1);  
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
    T0 = paramVec(13) - y(2) - y(3) - y(4) -y(5) - y(6) - y(7) - y(1); 
    % Sol= paramVec(14); 
    % kCATPon = paramVec(14);
    % kCATPoff = paramVec(15);
    % kMATPon = paramVec(16); &&&&&&&&&& 
    % kMATPoff = paramVec(17);
    kH = paramVec(18);
    kHYD = paramVec(19); 
    %JhydMyo =  kHYD *( ATP_myo / (kH + ATP_myo) );

    TimeFinal= Sol{1}; 
    CaSol = Sol{2};
    Ca= interp1(TimeFinal, CaSol,t); 

    dydt(1)= (kTon * Ca * T0) - ( kToff* y(1)) - (kTon*Ca *y(1)) + (kToff* y(2)) - (k0on * y(1)) + (k0off * y(4));
    dydt(2)= (kTon * Ca * y(1)) - ( kToff* y(2)) - (kCaOn * y(2)) + (kCaOff* y(5));
    dydt(3)=  -(kTon * Ca * y(3)) + ( kToff* y(4)) + (k0on * T0) - (k0off * y(3));
    dydt(4)=  (kTon * Ca  * y(3)) - ( kToff* y(4)) + (k0on * y(1)) - (k0off * y(4)) - (kTon *Ca*y(4)) + (kToff* y(5));
    dydt(5)=  (kTon * Ca * y(4)) - ( kToff* y(5)) + (kCaOn * y(2)) - (kCaOff * y(5)) - (f0 * y(5)) + (fP * y(6)) + (g0* y(7));
    dydt(6)=  (f0* y(5)) - ( fP * y(6)) + (hP * y(7)) - (h0 * y(6));
    dydt(7)=  -(hP * y(7)) +  (h0 * y(6)) - (g0 * y(7));
    %atp: 
    dydt(8)= -JhydM - (kCATPon* Ca * ATPm - kCATPoff*CaATP) - (kMATPon* Mg * ATPm - kMATPoff*MgATPm ) + tATP*(ATPtm - ATPm)/Vm;
    
end
  