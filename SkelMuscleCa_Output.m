clc 

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
    0.3632;    % yinit(16) is the inital consition for 'CATP'
    ];
 
% % Parameter values
param = importdata('InputParam1.xlsx');
p =  param.data; 

load PSO_15-Mar-2024_NEW.mat pSol
pSol(44) = 1;
% load PSO_14-Mar-2024_2.mat pSol
% pSol(42:44) = 1;
p = p(:) .* pSol(:);


%% Steady State and Dynamics calculation
expt = 1;
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt); % Steady State values of variable 
[Time,Y, ~, fluxes] = SkelMuscleCa_dydt([0 60], 50, 0, ySS(end,:), p, tic, expt); % Dynamics computation

% Steady State and Dynamics calculation without any SOCE.
expt = 2;
[~,ySS_noSOCE] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt); % Steady State values of variable 
[Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE] = SkelMuscleCa_dydt([0 60], 50, 0, ySS_noSOCE(end,:), p, tic, expt); % Dynamics computation

%% Figure 3
Fig3a(Time,Y(:,5),Y(:,8))
Fig3b(Time,Y(:,5),Y(:,8)) % V_SL vs Ca2+_Myo Zoomedin
Fig3c(Time,Y(:,1),Y(:,2)) % Density of activated Orai1 channel vs SR [Ca^{2+}] 
Fig3d(Time,[fluxes(:,4), fluxes(:,3)]) % DHPR and NCX Myo Fluxes
Fig3e(Time,[fluxes(:,6), fluxes(:,7), fluxes(:,8)]) % SR fluxes
%% Subplot for flux comparision of different channels in SRM and Sarcolemma
figure
subplot(1,3,1)
plot(Time,fluxes(:,4))
hold on
plot(Time,fluxes(:,3))
hold off
xlim([0 1]);
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (\muM/s)','Fontsize',16)
legend('J_{DHPR}','J_{NCX}','Fontsize',14)

subplot(1,3,2)
plot(Time,fluxes(:,6))
hold on
plot(Time,fluxes(:,7))
hold on
plot(Time,fluxes(:,8))
hold off
xlim([0 1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (\muM/s)','Fontsize',16)
legend('J_{RyR}','J_{SERCA}','J_{SRLeak}','Fontsize',14)

subplot(1,3,3)
plot(Time,fluxes(:,1))
hold on
plot(Time,fluxes(:,5))
hold off
xlim([0 1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (\muM/s)','Fontsize',16)
legend('J_{SOCE}','J_{PMCA}','Fontsize',14)

%% Figure 4
%t = 0 : 2.5 : max(Time);
t = 0:max(Time);

% Peak Ca2+ per stimulus plot
MaxCa = zeros(length(t),1);
MaxCa_noSOCE = zeros(length(t),1);
for j = 1 : (length(t) - 1)
    t_index = (Time > t(j)) & (Time < t(j+1));
    timepts = Time(t_index);
    MaxCa(j) = max(Y(t_index,8));  

    t_index_noSOCE = (Time_noSOCE > t(j)) & (Time_noSOCE < t(j+1));
    timepts_noSOCE = Time_noSOCE(t_index_noSOCE);
    MaxCa_noSOCE(j) = max(Y_noSOCE(t_index_noSOCE,8));  
end

figure
subplot1 = subplot(2,1,1);
scatter(t(2:end),MaxCa(1:end-1),'filled')
hold on
scatter(t(2:end),MaxCa_noSOCE(1:end-1),'filled')
xlabel('Time (s)','FontSize',16);
ylabel('Max [Ca^{2+}_{myo}] (\muM)','FontSize',16);
title ('Max [Ca^{2+}_{myo}] per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold');
legend('Control','Orai1 blocked','FontSize',14,'EdgeColor','none');
set(subplot1,'FontSize',14)
hold off

% Tail integral per stimulus plot
AUC = zeros(length(t),1);
AUC_noSOCE = zeros(length(t),1);
for k = 1 : (length(t) - 1)
    t_index = (Time >= (t(k) + 0.5)) & (Time < t(k+1));
    timepts = Time(t_index);
    AUC(k) = trapz(timepts,Y(t_index,8));   

    t_index_noSOCE = (Time_noSOCE >= (t(k) + 0.5)) & (Time_noSOCE < t(k+1));
    timepts_noSOCE = Time_noSOCE(t_index_noSOCE);
    AUC_noSOCE(k) = trapz(timepts_noSOCE,Y_noSOCE(t_index_noSOCE,8));   
end

subplot2 = subplot(2,1,2);
scatter(t(2:end),AUC(1:end-1),'filled')
hold on
scatter(t(2:end),AUC_noSOCE(1:end-1),'filled')
xlabel('Time (s)','FontSize',16);
ylabel('Tail Integral of [Ca^{2+}]_{myo} (\muM)','FontSize',16);
title ('Tail Integral of [Ca^{2+}]_{myo} per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold')
legend('Control','Orai1 blocked','FontSize',14,'EdgeColor','none')
set(subplot2,'FontSize',14)
hold off


