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
    10.004;    % yinit(17) is the initial condition for 'CaTrop'
    ];

% Parameter values
% param = importdata('InputParam1.xlsx');

p =  param.data;
load PSO_03-Apr-2024.mat pSol
p = p(:) .* pSol(:);

%% Steady State and Dynamics calculation for expt stimulus of 50Hz for 0.5s every 2.5s----------------------------------
expt = 1;
[~,ySS1] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt);                                                                       % Steady State values of variable
[Time1,Y1, ~, fluxes1,currents1] = SkelMuscleCa_dydt([0 150], 50, 0, ySS1(end,:), p, tic, expt);                                        % Dynamics computation

%% Steady State and Dynamics calculation without any SOCE for expt stimulus of 50Hz for 0.5s every 2.5s -----------------
expt = 0;
[~,ySS_noSOCE0] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt);                                                                % Steady State values of variable
[Time_noSOCE0,Y_noSOCE0, ~, fluxes_noSOCE0,currents_noSOCE0] = SkelMuscleCa_dydt([0 150], 50, 0, ySS_noSOCE0(end,:), p, tic, expt);     % Dynamics computation

%% Steady State and Dynamics calculation for Cont Stimulation of 10Hz for 30s
expt = 2;
[~,ySS2] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt);                    % Steady State values of variable
[Time2,Y2, ~, fluxes2,currents2] = SkelMuscleCa_dydt([0 30], 10, 0, ySS2(end,:), p, tic, expt);   % Dynamics computation
%% Steady State and Dynamics calculation without any SOCE for Cont Stimulation of 10Hz for 30s
expt = 4;
[~,ySS_noSOCE4] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, expt);                                                                % Steady State values of variable
[Time_noSOCE4,Y_noSOCE4, ~, fluxes_noSOCE4,currents_noSOCE4] = SkelMuscleCa_dydt([0 30], 10, 0, ySS_noSOCE4(end,:), p, tic, expt);      % Dynamics computation
%% Figure 3 ---------------------------------------------------------------
%Fluxes = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR]
%currents = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
%Only expt 2 is used for fig3. 

Fig3a(Time2,Y2(:,5),Y2(:,8))                            % V_SL and [Ca^2+]_myo over time

%% Stimulus Plot 

plot1 = figure;
axes1  = axes('Parent',plot1);
plot(Time2,currents2(:,13),"LineWidth",2)
ylabel('I_{SL} (pA)','Fontsize',16)
%xlim(axes1,[0 0.1])
set(axes1,'FontSize',16,'Box', 'off','FontSmoothing','on');
set(plot1,"Renderer","painters");

%% Myoplasmic Calcium zoomed in [0 0.1]
plot2 = figure;
axes2 = axes('Parent', plot2);
plot(Time2,Y2(:,8),'LineWidth',2,'Color',[0.635294117647059 0.0784313725490196 0.184313725490196]);                           % Ca2+_Myo Zoomedin
xlim(axes2,[0 0.1])
%xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{myo} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes2,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot2,"Renderer","painters");

%% Fluxes plots zoomed in [0 0.1]
plot3 = figure;
axes3 = axes('Parent', plot3);
plot(Time2,fluxes2(:,1),'LineWidth',2,'Color',[0.72,0.27,1.00])
hold on
plot(Time2,fluxes2(:,2),'LineWidth',2,'Color',[0.93,0.69,0.13])
hold on
plot(Time2,fluxes2(:,3),'LineWidth',2,'Color',[0.77,0.11,0.23])
hold on
plot(Time2,fluxes2(:,5),'LineWidth',2,'Color',[0.49,0.18,0.56])
hold off
xlim([0 0.1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (μM/s)','Fontsize',16)
legend('SOCE','Leak_{myo}','NCX','PMCA','Fontsize',14,'Edgecolor','none', 'Color','none')
set(axes3,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot3,"Renderer","painters");
%-------------------------------------------------------------------------------
plot4 = figure;
axes4 = axes('Parent', plot4);
plot(Time2,fluxes2(:,4),'LineWidth',2);
hold on
plot(Time2,fluxes2(:,6),'LineWidth',2)
hold off
xlim([0 0.1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (μM/s)','Fontsize',16)
legend('DHPR','RyR','Fontsize',14,'Edgecolor','none', 'Color','none')
set(axes4,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot4,"Renderer","painters");
%-------------------------------------------------------------------------------
plot5 = figure;
axes5 = axes('Parent', plot5);
plot(Time2,fluxes2(:,7),'LineWidth',2)
hold on
plot(Time2,fluxes2(:,8),'LineWidth',2)
hold off
xlim([0 0.1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (μM/s)','Fontsize',16)
legend('SERCA','Leak_{SR}','Fontsize',14,'Edgecolor','none', 'Color','none')
set(axes5,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot5,"Renderer","painters");
%% SR Calcium over time
plot6 = figure;
axes6 = axes('Parent', plot6);
plot(Time2,Y2(:,2),'LineWidth',2);                           
xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{SR} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes6,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot6,"Renderer","painters");

%%
Fig3c(Time2,Y2(:,1),Y2(:,2))                           % Density of activated Orai1 channel vs SR [Ca^{2+}]
Fig3d(Time2,[fluxes2(:,4), fluxes2(:,3)])              % DHPR and NCX Myo Fluxes
Fig3e(Time2,[fluxes2(:,6), fluxes2(:,7), fluxes2(:,8)]) % SR fluxes

%% Subplot for flux comparision of different channels in SRM and Sarcolemma
% figure
% subplot(1,3,1)
% plot(Time,fluxes(:,4))
% hold on
% plot(Time,fluxes(:,3))
% hold off
% xlim([0 0.01])
% xlabel('Time (s)', 'Fontsize',16)
% ylabel('Flux (\muM/s)','Fontsize',16)
% legend('J_{DHPR}','J_{NCX}','Fontsize',14)

% subplot(1,3,2)
% semilogy(Time,fluxes(:,6))
% hold on
% semilogy(Time,fluxes(:,7))
% hold on
% semilogy(Time,fluxes(:,8))
% hold off
% xlim([0 1])
% xlabel('Time (s)', 'Fontsize',16)
% ylabel('Flux (\muM/s)','Fontsize',16)
% legend('J_{RyR}','J_{SERCA}','J_{SRLeak}','Fontsize',14)

% subplot(1,3,3)
% plot(Time,fluxes(:,1))
% hold on
% plot(Time,fluxes(:,5))
% hold off
% xlim([0 1])
% xlabel('Time (s)', 'Fontsize',16)
% ylabel('Flux (\muM/s)','Fontsize',16)
% legend('J_{SOCE}','J_{PMCA}','Fontsize',14)

figure
semilogy(Time,fluxes(:,6))
hold on
plot(Time,fluxes(:,7))
hold on
plot(Time,fluxes(:,8))
hold on
plot(Time_noSOCE,fluxes_noSOCE(:,6))
hold on
plot(Time_noSOCE,fluxes_noSOCE(:,7))
hold on
plot(Time_noSOCE,fluxes_noSOCE(:,8))
hold off
xlim([0 0.01])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (\muM/s)','Fontsize',16)
legend('Control_{RyR}','Control_{SERCA}','Control_{SRLeak}','Orai1 Blocked_{RyR}','Orai1 Blocked_{SERCA}','Orai1 Blocked_{SRLeak}','Fontsize',14,'Edgecolor','none','Color','none')

figure
% subplot(1,3,1)
plot(Time,fluxes(:,4))
hold on
plot(Time,fluxes(:,3))
hold on
plot(Time_noSOCE,fluxes_noSOCE(:,4))
hold on
plot(Time_noSOCE,fluxes_noSOCE(:,3))
hold off
xlim([0 0.01])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (\muM/s)','Fontsize',16)
legend('Control_{DHPR}','Control_{NCX}','Orai1 Blocked_{DHPR}','Orai1 Blocked_{NCX}','Fontsize',14,'Edgecolor','none','Color','none')

%% Figure 4 ---------------------------------------------------------------
t = 0 : 2.5 : max(Time);
%t = 0:max(Time);

% Peak Ca2+ per stimulus plot----------------------------------------------
MaxCa = zeros(length(t),1);
MaxCa_noSOCE = zeros(length(t),1);
for j = 1 : (length(t) - 1)
    t_index = (Time > t(j)) & (Time < t(j+1));
    timepts = Time(t_index);
    MaxCa(j) = max(Y(t_index,8));
    % When Orai1 channel is blocked
    t_index_noSOCE = (Time_noSOCE > t(j)) & (Time_noSOCE < t(j+1));
    timepts_noSOCE = Time_noSOCE(t_index_noSOCE);
    MaxCa_noSOCE(j) = max(Y_noSOCE(t_index_noSOCE,8));
end

subplot1 = figure; %subplot(2,1,1);
scatter(1:length(MaxCa)-1,MaxCa(1:end-1),'filled') %t(2:end)
hold on
scatter(1:length(MaxCa_noSOCE)-1,MaxCa_noSOCE(1:end-1),'filled')
xlabel('Consecutive stimuli','FontSize',16); %Time (s)
ylabel('Max [Ca^{2+}_{myo}] (\muM)','FontSize',16);
title ('Max [Ca^{2+}_{myo}] per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold');
legend('Control','Orai1 blocked','FontSize',14,'EdgeColor','none');
set(subplot1,'FontSize',14)
hold off

%% Tail integral per stimulus plot -----------------------------------------
AUC = zeros(length(t),1);
AUC_noSOCE = zeros(length(t),1);
for k = 1 : (length(t) - 1)
    t_index = (Time >= (t(k) + 0.5)) & (Time < t(k+1));
    timepts = Time(t_index);
    AUC(k) = trapz(timepts,Y(t_index,8));
    % When Orai1 channel is blocked
    t_index_noSOCE = (Time_noSOCE >= (t(k) + 0.5)) & (Time_noSOCE < t(k+1));
    timepts_noSOCE = Time_noSOCE(t_index_noSOCE);
    AUC_noSOCE(k) = trapz(timepts_noSOCE,Y_noSOCE(t_index_noSOCE,8));
end

subplot2 = figure; %subplot(2,1,2);

scatter(t(2:end),AUC(1:end-1),'filled')
hold on
scatter(t(2:end),AUC_noSOCE(1:end-1),'filled')
xlabel('Time (s)','FontSize',16);
ylabel('Tail Integral of [Ca^{2+}]_{myo} (\muM)','FontSize',16);
title ('Tail Integral of [Ca^{2+}]_{myo} per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold')
legend('Control','Orai1 blocked','FontSize',14,'EdgeColor','none')
set(subplot2,'FontSize',14)
hold off
% -------------------------------------------------------------------------

