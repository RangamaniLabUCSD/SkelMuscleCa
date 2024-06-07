%% Color code
Delta = [0.89,0.47,0.97]; % - Pink shade for delta
Ca_control = [0.64,0.08,0.18]; % - Maroon shade for control Calcium
SOCE_Blocked = [0.10,0.85,0.83]; % - Cyan shade for SOCE Blocked calcium
V = [0.49,0.18,0.56] ; %-Dark Purple for Voltage


% Flux# = [J_SOCE, J_CaLeak_SL , J_NCX_C, J_DHPR, J_PMCA, LumpedJ_RyR, LumpedJ_SERCA, J_CaLeak_SR]
% Current# = [I_CaLeak_SL, I_Cl, I_DHPR, I_K_DR, I_K_IR, I_NCX_C, I_NCX_N, I_NKX_K, I_NKX_N, I_Na, I_PMCA, I_SOCE, I_SL];
% Expt # : 
% 1,2 - Cont Stim (Control & no SOCE)
% 3,4 - Wei-Lapierre (Control & noSOCE)
% 5,6 - Resistance (Control & noSOCE)
% 7,8 - HIIT (Control & noSOCE)

% Load data for plotting 

% Define variables below
%Current = ;
%Y = ;
%Flux = ;


%% Input Stimulus
plot1 = figure;
axes1 = axes('Parent', plot1);
plot(Time,Current(:,end),'LineWidth',2);                          
xlabel('Time (s)', 'Fontsize',16)
ylabel('I_{SL (pA)','Fontsize',16,'FontSmoothing','on')
set(axes1,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot1,"Renderer","painters");

%% Myoplasmic Calcium 
plot2 = figure;
axes2 = axes('Parent', plot2);
plot(Time,Y(:,8),'LineWidth',2,'Color',Ca_control);                           
xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{myo} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes2,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot2,"Renderer","painters");

%% Membrane Potential
plot3 = figure;
axes3 = axes('Parent', plot3);
plot(Time,Y(:,5),'LineWidth',2,'Color',V)
xlabel('Time (s)', 'Fontsize',16)
ylabel('V_{SL} (mV)','Fontsize',16)
set(axes3,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot3,"Renderer","painters");
%% DHPR - RyR Flux --------------------------------------------------------
plot4 = figure;
axes4 = axes('Parent', plot4);
semilogy(Time2,Flux(:,4),'LineWidth',2);
hold on
semilogy(Time2,Flux(:,6),'LineWidth',2)
hold off
xlim([0 0.1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (μM/s)','Fontsize',16)
legend('DHPR','RyR','Fontsize',14,'Edgecolor','none', 'Color','none')
set(axes4,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot4,"Renderer","painters");

% SERCA and SR Leak -------------------------------------------------------
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

% Myo Fluxes --------------------------------------------------------------
plot6 = figure;
axes6 = axes('Parent',plot6);
plot(Time, Flux(:,1:3),'LineWidth',2)
hold on
plot(Time, Flux(:,5),'LineWidth',2)
hold off
xlim([0 0.1])
xlabel('Time (s)', 'Fontsize',16)
ylabel('Flux (μM/s)','Fontsize',16)
legend('SOCE','Leak_{SL}','NCX','PMCA','Fontsize',14,'Edgecolor','none', 'Color','none')
set(axes5,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot5,"Renderer","painters");

%% SR Calcium 
plot7 = figure;
axes7 = axes('Parent', plot7);
plot(Time2,Y(:,2),'LineWidth',2);                           
xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{SR} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes7,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot7,"Renderer","painters");

%% Comparison of Calcium in presence and absence of SOCE
% Cont stimulus (expt 1 and 2)
plot8_SR = figure;
axes8_SR = axes('Parent', plot8_SR);
plot(Time1,Y1(:,2),'LineWidth',2);
hold on 
plot(Time_noSOCE2,Y_noSOCE2(:,2),'LineWidth',2);
xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{SR} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes8_SR,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot8_SR,"Renderer","painters");
% -------------------------------------------------------------------------
plot8_myo = figure;
axes8_myo = axes('Parent', plot8_myo);
plot(Time1,Y1(:,2),'LineWidth',2,'Color',Ca_control);
hold on 
plot(Time_noSOCE2,Y_noSOCE2(:,2),'LineWidth',2,'Color',SOCE_Blocked);
xlabel('Time (s)', 'Fontsize',16)
ylabel('[Ca^{2+}]_{myo} (μM)','Fontsize',16,'FontSmoothing','on')
set(axes8_myo,'FontSize',16,'Box','off','FontSmoothing','on')
set(plot8_myo,"Renderer","painters");

%% Wei La-Pierre expt comparison-------------------------------------------
t = 0 : 2.5 : max(Time);

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

figure
plot9_control = figure;
axes9_control = axes('Parent', plot9_control);
scatter(t(2:end),MaxCa(1:end-1),'filled','MarkerFaceColor',Ca_control)
hold on
scatter(t(2:end),MaxCa_noSOCE(1:end-1),'filled','MarkerFaceColor',SOCE_Blocked)
xlabel('Time (s)','FontSize',16);
ylabel('Max [Ca^{2+}_{myo}] (μM)','FontSize',16);
title ('Max [Ca^{2+}_{myo}] per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold');
legend('Control','SOCE blocked','FontSize',14,'EdgeColor','none');
set(axes9_control,'FontSize',14)
hold off
set(plot9_control,"Renderer","painters");

% Tail integral per stimulus plot -----------------------------------------
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

plot9_SOCEblocked = figure;
axes9_SOCEblocked = axes('Parent', plot9_SOCEblocked);
subplot2 = subplot(2,1,2);
scatter(t(2:end),AUC(1:end-1),'filled','MarkerFaceColor',Ca_control)
hold on
scatter(t(2:end),AUC_noSOCE(1:end-1),'filled','MarkerFaceColor',SOCE_Blocked)
xlabel('Time (s)','FontSize',16);
ylabel('Tail Integral of [Ca^{2+}]_{myo} (μM)','FontSize',16);
title ('Tail Integral of [Ca^{2+}]_{myo} per stimulus period of 2.5s ','FontSize',16,'FontWeight','bold')
legend('Control','Orai1 blocked','FontSize',14,'EdgeColor','none')
set(axes9_SOCEblocked,'FontSize',14)
hold off
set(plot9_SOCEblocked,"Renderer","painters");
% -------------------------------------------------------------------------

%% Load exercise data Resistance.mat HIIT.mat

%% Resistance Plots

% DeltaMax Resistance ----------------------------------------------------
plot10 = figure;
axes10 = axes('Parent', plot10);
scatter(f_Resistance,deltaMax_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes10,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot10,"Renderer","painters");

% DeltaAuc Resistannce ---------------------------------------------------
plot11 = figure;
axes11 = axes('Parent', plot11);
scatter(f_Resistance,deltaAUC_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes11,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot11,"Renderer","painters");

%  Avg Ca ----------------------------------------------------------------
plot12 = figure;
axes12 = axes('Parent', plot12);
scatter(f_Resistance,AUC5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,AUC_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes12,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot12,"Renderer","painters");

% Max Ca -----------------------------------------------------------------
plot13 = figure;
axes13 = axes('Parent', plot13);
scatter(f_Resistance,MaxCa5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,MaxCa_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes13,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot13,"Renderer","painters");

% Force vs freq -----------------------------------------------------------
plot14 = figure;
axes14 = axes('Parent', plot14);
scatter(f_Resistance,Force5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,Force_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Force (%)','Fontsize',18)
title('Average force during Resistance exercise')
set(axes14,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot14,"Renderer","painters");

% Delta F vs freq ---------------------------------------------------------
plot15 = figure;
axes15 = axes('Parent', plot15);
scatter(f_Resistance,DeltaForce_Resistance,'filled','MarkerFaceColor',[1,0,1]); 
%ylim([0 100]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('% Change in Force','Fontsize',18)
title('Relative force during Resistance exercise')
set(axes15,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot15,"Renderer","painters");

%% HIIT Plots 

% DeltaAvg ----------------------------------------------------------------
plot16 = figure;
axes16 = axes('Parent', plot16);
scatter(f_HIIT,deltaAUC_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes16,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot16,"Renderer","painters");

% Avg Ca ------------------------------------------------------------------
plot17 = figure;
axes17 = axes('Parent', plot17);
scatter(f_HIIT,AUC7,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_HIIT,AUC_noSOCE8,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes17,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot17,"Renderer","painters");

% Max Ca ------------------------------------------------------------------
plot18 = figure;
axes18 = axes('Parent', plot18);
scatter(f_HIIT,MaxCa7,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_HIIT,MaxCa_noSOCE8,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes18,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot18,"Renderer","painters");

% Delta MaxCa -------------------------------------------------------------
plot19 = figure;
axes19 = axes('Parent', plot19);
scatter(f_HIIT,deltaMax_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes19,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot19,"Renderer","painters");

% Force vs freq -----------------------------------------------------------
plot20 = figure;
axes20 = axes('Parent', plot20);
scatter(f_HIIT,Force7,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_HIIT,Force_noSOCE8,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Normalized Force (%)','Fontsize',18)
title('Average force during HIIT exercise')
set(axes20,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot20,"Renderer","painters");

% Delta F vs freq ---------------------------------------------------------
plot21 = figure;
axes21 = axes('Parent', plot21);
scatter(f_HIIT,DeltaForce_HIIT,'filled','MarkerFaceColor',[1,0,1]); 
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('% Change in Force','Fontsize',18)
title('Change in force during HIIT exercise')
set(axes21,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot21,"Renderer","painters");

%% Ca vs Time
plot22 = figure;
axes22= axes('Parent', plot22);
plot(Time5{5}(indices1),Ca5{5}(indices1),'Color',[0.64,0.08,0.18]); 
hold on
plot(Time_noSOCE6{5}(indices2),Ca_noSOCE6{5}(indices2),'Color',[0.10,0.85,0.83]); 
xlabel('Time (s)', 'Fontsize',18)
ylabel('[Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes22,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot22,"Renderer","painters");
 