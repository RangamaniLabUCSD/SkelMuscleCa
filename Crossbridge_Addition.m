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
    0;%387;        % yinit(14) is the initial condition for 'CaParv'
    0;%1020;       % yinit(15) is the initial condition for 'MgParv'
    0.3632;     % yinit(16) is the initial consition for 'CATP'
    0;%10.004;     % yinit(17) is the initial condition for 'CaTrop'
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

load PSO_17-Sep-2024.mat pSol
highSensIdx = [1,5,15,16,19,22,24,30,32,34,35,43,74,83,91,92];
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).* pPSO(highSensIdx);
% number 69 is the rate of phosphate transport into SR
% pPSO(69) = pPSO(69)*100;

% compute steady state solution without SOCE or phosphate accumulation
[TimeSS,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pPSO, tic, 2, false); 
pPSO(12) = pPSO(12)*ySS(end,2); % set c_ref according to SS c_SR
% tSol = 0:.0001:10;
tSol = [0, 420]; % 420 is max time for HIIT
freq = 10; 
% ySS(end,26) = 0;
% ySS(end,28) = 2500;
% expt 1: const stim, with SOCE, expt 2: const stim, no SOCE
% expt 3: fig 3.4 stim, with SOCE, expt 4: fig 3.4 stim, no SOCE
% expt 1: HIIT stim, with SOCE, expt 2: HIIT stim, no SOCE
expt = 6;
phosphateAccum = false;
[Time,Y,~,fluxes,currents] = SkelMuscleCa_dydt(tSol, freq, 0, ySS(end,:), pPSO, tic, expt, phosphateAccum); % compute time-dependent solution

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

%%
figure
plot(Time, Y(:,2))
title('Time vs Ca2+ SR')
xlabel('Time (seconds)');
ylabel('[Ca2+]_{SR} (µM)'); 

figure 
plot(Time, Y(:,23))
% % % % 
% figure
% plot(Time, Y(:,23))
% title('Time vs ATP, freq=100 Hz')
% xlabel('Time (seconds)');
% ylabel('[ATP] (µM)'); 
% % % 
figure
plot(Time, Y(:,21))
title('Time vs Force, freq=100 Hz')
xlabel('Time (seconds)');
ylabel('Post Power Stroke'); 

figure
plot(Time, Y(:,5))
title('Time vs SL Voltage ')
xlabel('Time (seconds)');
ylabel('SL Voltage'); 

 
% Input Stimulus
plot1 = figure;
axes1 = axes('Parent', plot1);
plot(Time,currents(:,end),'LineWidth',1);                          
xlabel('Time (s)', 'Fontsize',16)
ylabel('I_{SL} (pA)','Fontsize',16,'FontSmoothing','on')
set(axes1,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot1,"Renderer","painters");

figure 
yyaxis left
plot(Time,currents(:,end) / 1e9)
ylim([-1e-4 1e-4])
yyaxis right
plot(Time,Y(:,5))
xlim([0 0.05])

% 
% figure 
% yyaxis left
% plot(Time,Y(:,8))
% ylabel('[Ca2+]_{SR} (µM)'); 
% ylim([0 20])
% yyaxis right
% plot(Time,Y(:,21))
% ylabel('Post Power Stroke'); 
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
title('time vs p_i_{SR}, freq=100 Hz')
legend('p_i_{SR}');  
xlabel('Time (seconds)');
ylabel('[p_i_{SR}] (µM)'); 
subplot(3,1,2)
plot(Time, Y(:,25))
title('time vs PiCa, freq=100 Hz')
legend('PiCa Conc');  
xlabel('Time (seconds)');
ylabel('[PiCa] (µM)'); 
subplot(3,1,3)
plot(Time, Y(:,26))
title('time vs Pi_{Myo}, freq=100 Hz')
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


%% Adding Senneff Corssbridge Cycling
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
  