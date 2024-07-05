clc
clear
%delete(gcp('nocreate'))
load PSO_25-Apr-2024.mat p pSol 
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
    0.01        % yinit(26) is the initial condition for 'p_i_SR' 
    0.001       % yinit(27) is the initial condition for 'PiCa'
    0.001       % yinit(28) is the initial condition for 'Pi_Myo'
    ];

p = p(:) .* pSol(:);
% Expt : 
% 1,2 - Cont Stim (Control & no SOCE)
% 3,4 - Wei-Lapierre (Control & noSOCE)
% 5,6 - Resistance (Control & noSOCE)
% 7,8 - HIIT (Control & noSOCE)


f_HIIT = 60 : 10 : 180;
l_HIIT = length(f_HIIT);
t_HIIT = [0 60];

Time7 = cell(l_HIIT,1);
Ca7 = cell(l_HIIT,1);
Flux7 = cell(l_HIIT,1);
Current7 = cell(l_HIIT,1);
MaxCa7 = zeros(l_HIIT,1);
AUC7 = zeros(l_HIIT,1); 
Force7 = cell(l_HIIT,1);
Max_F7 = zeros(l_HIIT,1);
AUC_F7 = zeros(l_HIIT,1);

Time_noSOCE8 = cell(l_HIIT,1);
Ca_noSOCE8 = cell(l_HIIT,1);
Flux_noSOCE8 = cell(l_HIIT,1);
Current_noSOCE8 = cell(l_HIIT,1);
MaxCa_noSOCE8 = zeros(l_HIIT,1);
AUC_noSOCE8 = zeros(l_HIIT,1);
Force_noSOCE8 = cell(l_HIIT,1);
MaxF_noSOCE8 = zeros(l_HIIT,1);
AUC_F_noSOCE8 = zeros(l_HIIT,1);

f_Resistance = 10:2:50;
l_Resistance = length(f_Resistance);
t_Resistance = [0 420];

Time5 = cell(l_Resistance,1);
Ca5 = cell(l_Resistance,1);
Flux5 = cell(l_Resistance,1);
Current5 = cell(l_Resistance,1);
MaxCa5 = zeros(l_Resistance,1);
AUC5 = zeros(l_Resistance,1);
Force5  = cell(l_Resistance,1);
Max_F5 = zeros(l_Resistance,1);
AUC_F5 = zeros(l_Resistance,1);

Time_noSOCE6 = cell(l_Resistance,1);
Ca_noSOCE6 = cell(l_Resistance,1);
Flux_noSOCE6 = cell(l_Resistance,1);
Current_noSOCE6 = cell(l_Resistance,1);
MaxCa_noSOCE6 = zeros(l_Resistance,1);
AUC_noSOCE6 = zeros(l_Resistance,1);
Force_noSOCE6 = cell(l_Resistance,1);
MaxF_noSOCE6 = zeros(l_Resistance,1);
AUC_F_noSOCE6 = zeros(l_Resistance,1);

% Hill equation variables
n = 3.5947; % Hill coeff
Ca50 = 1.1694; % Calcium conc at 50% max force.
c0 = 1.03167e-04 ; % Vertical shift 
%% Steady State calculation -----------------------------------------------
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 1);                    % Steady State values of variable
yinf = ySS(end,:);
[~,ySS_noSOCE] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 2);  
yinf_noSOCE = ySS_noSOCE(end,:);

% %% HIIT -------------------------------------------------------------------
% parpool(13)
parfor i = 1:l_HIIT
    freq = f_HIIT(i);

    expt = 7;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt(t_HIIT, freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time7{i} = Time;
    Ca7{i} = Y(:,8);
    Flux7{i} = fluxes;
    Current7{i} = currents;
    Force7{i} = Y(:,23);

    expt = 8;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt(t_HIIT, freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE8{i} = Time_noSOCE;
    Ca_noSOCE8{i} = Y_noSOCE(:,8);
    Flux_noSOCE8{i} = fluxes_noSOCE;
    Current_noSOCE8{i} = currents_noSOCE;
    Force_noSOCE8{i} = Y_noSOCE(:,23);

end

% Max and Avg

for i = 1:l_HIIT
    time_start7 = find(Time7{i} < 1);
    start_index7 = length(time_start7);
    MaxCa7(i) = max(Ca7{i}(start_index7:end));
    AUC7(i) = trapz(Time7{i},Ca7{i})/ max(Time7{i}); 
    Max_F7(i) = max(Force7{i}(start_index7:end));
    AUC_F7(i) = trapz(Time7{i},Force7{i} / max(Time7{i}));

    time_start8 = find(Time_noSOCE8{i} < 1);
    start_index8 = length(time_start8);
    MaxCa_noSOCE8(i) = max(Ca_noSOCE8{i}(start_index8:end));
    AUC_noSOCE8(i) =  trapz(Time_noSOCE8{i},Ca_noSOCE8{i})/max(Time_noSOCE8{i}); 
    MaxF_noSOCE8(i) = max(Force_noSOCE8{i}(start_index8:end));
    AUC_F_noSOCE8(i) =  trapz(Time_noSOCE8{i},Force_noSOCE8{i})/max(Time_noSOCE8{i}); 

end

deltaMax_HIIT = MaxCa7 - MaxCa_noSOCE8;
deltaAUC_HIIT = (AUC7 - AUC_noSOCE8);

deltaMax_F_HIIT = Max_F7 - MaxF_noSOCE8;
deltaAUC_F_HIIT = (AUC_F7 - AUC_F_noSOCE8);

 
% %Force
% for i = 1 : length(f_HIIT)
%     Force7(i) = c0 + (100 * ( AUC7(i)^n )/ ( (Ca50 ^ n)  + (AUC7(i) ^ n)) );
%     Force_noSOCE8(i) = c0 + (100 * (AUC_noSOCE8(i)^n )/ ( (Ca50 ^ n)  + (AUC_noSOCE8(i) ^ n)) );  
% end
% 
% DeltaForce_HIIT = Force7 - Force_noSOCE8;


%% Resistance -------------------------------------------------------------
% delete(gcp('nocreate'))
% parpool(21)
parfor i = 1:l_Resistance
    freq = f_Resistance(i);

    expt = 5;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt(t_Resistance, freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time5{i} = Time;
    Ca5{i} = Y(:,8);
    Flux5{i} = fluxes;
    Current5{i} = currents;
    Force5{i} = Y(:,23);

    expt = 6;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt(t_Resistance, freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE6{i} = Time_noSOCE;
    Ca_noSOCE6{i} = Y_noSOCE(:,8);
    Flux_noSOCE6{i} = fluxes_noSOCE;
    Current_noSOCE6{i} = currents_noSOCE;
    Force_noSOCE6{i} = Y_noSOCE(:,23);

end

%% Max and Avg

for i = 1:l_Resistance

    indices5 = (Time5{i} >=8 & Time5{i} <= 60)| (Time5{i} >= 188 & Time5{i} <= 240) | (Time5{i} >= 368 & Time5{i} <= 420);
    indices6 = (Time_noSOCE6{i} >= 8 & Time_noSOCE6{i} <= 60)| (Time_noSOCE6{i}  >= 188 & Time_noSOCE6{i}  <= 240) | (Time_noSOCE6{i}  >= 368 & Time_noSOCE6{i}  <= 420);
    MaxCa5(i) = max(Ca5{i}(indices5));
    MaxCa_noSOCE6(i) = max(Ca_noSOCE6{i}(indices6));
    Max_F5(i) = max(Force5{i}(indices5));
    MaxF_noSOCE6(i) = max(Force_noSOCE6{i}(indices6));


    AUC_indices5 = (Time5{i} >=0 & Time5{i} <= 60)| (Time5{i} >= 180 & Time5{i} <= 240) | (Time5{i} >= 360 & Time5{i} <= 420);
    AUC_indices6 = (Time_noSOCE6{i} >= 0 & Time_noSOCE6{i} <= 60)| (Time_noSOCE6{i}  >= 180 & Time_noSOCE6{i}  <= 240) | (Time_noSOCE6{i}  >= 360 & Time_noSOCE6{i}  <= 420);
    StimTime = 180; % Only considering the period of stimulation 
    AUC5(i) = trapz(Time5{i}(AUC_indices5),Ca5{i}(AUC_indices5)) / StimTime; 
    AUC_noSOCE6(i) = trapz(Time_noSOCE6{i}(AUC_indices6),Ca_noSOCE6{i}(AUC_indices6)) / StimTime; 
    AUC_F5(i) = trapz(Time5{i}(AUC_indices5),Force5{i}(AUC_indices5)) / StimTime; 
    AUC_F_noSOCE6(i) = trapz(Time_noSOCE6{i}(AUC_indices6),Force_noSOCE6{i}(AUC_indices6)) / StimTime; 


end
deltaMax_Resistance = MaxCa5 - MaxCa_noSOCE6;
deltaAUC_Resistance = AUC5 - AUC_noSOCE6; 

deltaMax_F_Resistance = Max_F5 - MaxF_noSOCE6;
deltaAUC_F_Resistance = AUC_F5 - AUC_F_noSOCE6; 
 

% % Force  
% for i = 1 : length(f_Resistance)
%     Force5(i) = c0 + (100 * ( AUC5(i)^n )/ ( (Ca50 ^ n)  + (AUC5(i) ^ n)) );
%     Force_noSOCE6(i) = c0 + (100 * (AUC_noSOCE6(i)^n )/ ( (Ca50 ^ n)  + (AUC_noSOCE6(i) ^ n)) );
% end
% DeltaForce_Resistance =  (Force5 - Force_noSOCE6);


%% Resistance Plots

% DeltaMax Resistance ----------------------------------------------------
plot10 = figure;
axes10 = axes('Parent', plot10);
scatter(f_Resistance,deltaMax_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes10,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot10,"Renderer","painters");

% DeltaAuc Resistance ---------------------------------------------------
plot11 = figure;
axes11 = axes('Parent', plot11);
scatter(f_Resistance,deltaAUC_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes11,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot11,"Renderer","painters");

%  Avg Ca vs freq ----------------------------------------------------------------
plot12 = figure;
axes12 = axes('Parent', plot12);
scatter(f_Resistance,AUC5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,AUC_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes12,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot12,"Renderer","painters");

% Max Ca vs freq -----------------------------------------------------------------
plot13 = figure;
axes13 = axes('Parent', plot13);
scatter(f_Resistance,MaxCa5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,MaxCa_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes13,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot13,"Renderer","painters");

% Avg Force vs freq -----------------------------------------------------------
plot14 = figure;
axes14 = axes('Parent', plot14);
scatter(f_Resistance,AUC_F5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,AUC_F_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Average Force','Fontsize',18)
title('Average force during Resistance exercise')
set(axes14,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot14,"Renderer","painters");

% Max Force vs freq -----------------------------------------------------------
plot14 = figure;
axes14 = axes('Parent', plot14);
scatter(f_Resistance,Max_F5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,MaxF_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Max Force','Fontsize',18)
title('Max force during Resistance exercise')
set(axes14,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot14,"Renderer","painters");

% DeltaAUC F vs freq ---------------------------------------------------------
plot15 = figure;
axes15 = axes('Parent', plot15);
scatter(f_Resistance,deltaAUC_F_Resistance,'filled','MarkerFaceColor',[1,0,1]); 
%ylim([0 100]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average Force','Fontsize',18)
title('Average force during Resistance exercise')
set(axes15,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot15,"Renderer","painters");

% DeltaMax Force Resistance ----------------------------------------------------
plot10 = figure;
axes10 = axes('Parent', plot10);
scatter(f_Resistance,deltaMax_F_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max Force','Fontsize',18)
set(axes10,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot10,"Renderer","painters");


%% HIIT Plots 
% DeltaAvg ----------------------------------------------------------------
plot16 = figure;
axes16 = axes('Parent', plot16);
scatter(f_HIIT,deltaAUC_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes16,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot16,"Renderer","painters");

% Delta MaxCa -------------------------------------------------------------
plot19 = figure;
axes19 = axes('Parent', plot19);
scatter(f_HIIT,deltaMax_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes19,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot19,"Renderer","painters");

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

% Avg Force vs freq -----------------------------------------------------------
plot20 = figure;
axes20 = axes('Parent', plot20);
scatter(f_HIIT,AUC_F7,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_HIIT,AUC_F_noSOCE8,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Normalized (AUC) Force (%)','Fontsize',18)
title('Average force during HIIT exercise')
set(axes20,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot20,"Renderer","painters");

% Max Force vs freq -----------------------------------------------------------
plot20 = figure;
axes20 = axes('Parent', plot20);
scatter(f_HIIT,Max_F7,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_HIIT,MaxF_noSOCE8,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
hold off
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Normalized (Max) Force (%)','Fontsize',18)
title('Average force during HIIT exercise')
set(axes20,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot20,"Renderer","painters");

% DeltaAUC F vs freq ---------------------------------------------------------
plot21 = figure;
axes21 = axes('Parent', plot21);
scatter(f_HIIT,deltaAUC_F_HIIT,'filled','MarkerFaceColor',[1,0,1]); 
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta AUC Force','Fontsize',18)
title('Change in AUC Force during HIIT exercise')
set(axes21,'FontSize',18,'Box','off','FontSmoothing','on')
set(plot21,"Renderer","painters");

% DeltaMax F vs freq ---------------------------------------------------------
plot21 = figure;
axes21 = axes('Parent', plot21);
scatter(f_HIIT,deltaMax_F_HIIT,'filled','MarkerFaceColor',[1,0,1]); 
xlim([60 180]);
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta AUC Force','Fontsize',18)
title('Change in Max Force during HIIT exercise')
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
 
save('HIIT_Resistance.mat')
