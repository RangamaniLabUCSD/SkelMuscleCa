clc
clear
load HIIT7_5-9.mat
load HIIT_noSOCE8_5-9.mat
load Resistance5_5-9.mat
load Resistance_noSOCE6_5-9.mat

%% Expt : 
%1,2 - Cont Stim (Control & no SOCE)
%3,4 - Wei-Lapierre (Control & noSOCE)
%5,6 - Resistance (Control & noSOCE)
%%7,8 - HIIT (Control & noSOCE)

%% Color code
%[0.89,0.47,0.97] - Pink shade for delta
%[0.64,0.08,0.18] - Maroon shade for control Calcium
%[0.10,0.85,0.83] - Cyan shade for SOCE Blocked calcium
%[0.49,0.18,0.56] -Dark Purple for Voltage

MaxCa5 = zeros(l_Resistance,1);
AUC5 = zeros(l_Resistance,1);

MaxCa_noSOCE6 = zeros(l_Resistance,1);
AUC_noSOCE6 = zeros(l_Resistance,1);

MaxCa7 = zeros(l_HIIT,1);
AUC7 = zeros(l_HIIT,1); 

MaxCa_noSOCE8 = zeros(l_HIIT,1);
AUC_noSOCE8 = zeros(l_HIIT,1);

%% Max and AUC -------------------------------------------------------------
for i = 1:l_HIIT
    time_start7 = find(Time7{i} < 1);
    start_index7 = length(time_start7);
    MaxCa7(i) = max(Ca7{i}(start_index7:end));
    AUC7(i) = trapz(Time7{i},Ca7{i}) / 60;

    time_start8 = find(Time_noSOCE8{i} < 1);
    start_index8 = length(time_start8);
    MaxCa_noSOCE8(i) = max(Ca7{i}(start_index7:end));
    AUC_noSOCE8(i) = trapz(Time7{i},Ca7{i}) / 60;

end

deltaMax_HIIT = MaxCa7 - MaxCa_noSOCE8;
deltaAUC_HIIT = AUC7 - AUC_noSOCE8;
%%
for i = 1:l_Resistance

    indices5 = (Time5{i} >=1 & Time5{i} <= 60) | (Time5{i} >= 181 & Time5{i} <= 240) | (Time5{i} >= 361 & Time5{i} <= 420);
    indices6 = (Time_noSOCE6{i} >= 1 & Time_noSOCE6{i} <= 60) | (Time_noSOCE6{i}  >= 181 & Time_noSOCE6{i}  <= 240) | (Time_noSOCE6{i}  >= 361 & Time_noSOCE6{i}  <= 420);
    % time_start5_1 = find((Time5{i} < 1);
    % start_index5_1 = length(time_start5_1);
    % time_start5_2 = find(Time5{i} < 60);
    % start_index5_2 = length(time_start5_2);
    % time_start5_3 = find(Time5{i} < 180 );
    % start_index5_3 = length(time_start5_3);
    % time_start5_4 = find(Time5{i} < 240 );
    % start_index5_4 = length(time_start5_4);
    % time_start5_5 = find(Time5{i} < 360);
    % start_index5_5 = length(time_start5_5);

    %MaxCa5(i) = max(Ca5{i}(start_index5_1:));
    MaxCa5(i) = max(Ca5{i}(indices5));
    AUC5(i) = trapz(Time5{i},Ca5{i}) / 420;

    % time_start6 = find(Time_noSOCE6{i} < 1);
    % start_index6 = length(time_start6);
    %MaxCa_noSOCE6(i) = max(Ca_noSOCE6{i}(start_index6:end));
    MaxCa_noSOCE6(i) = max(Ca_noSOCE6{i}(indices6));
    AUC_noSOCE6(i) = trapz(Time_noSOCE6{i},Ca_noSOCE6{i}) / 420;
 
end

deltaMax_Resistance = MaxCa5 - MaxCa_noSOCE6;
deltaAUC_Resistance = AUC5 - AUC_noSOCE6; 
%% Plots ------------------------------------------------------------------
figure1 = figure;
axes1 = axes('Parent', figure1);
scatter(f_Resistance,deltaMax_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes1,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure1,"Renderer","painters");

%% DeltaAuc Resistannce
figure3 = figure;
axes3 = axes('Parent', figure3);
scatter(f_Resistance,deltaAUC_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes3,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure3,"Renderer","painters");

%%  AVg Ca
figure3 = figure;
axes3 = axes('Parent', figure3);
scatter(f_Resistance,AUC5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,AUC_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes3,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure3,"Renderer","painters");
%% Max Ca
figure3 = figure;
axes3 = axes('Parent', figure3);
scatter(f_Resistance,MaxCa5,'filled','MarkerFaceColor',[0.64,0.08,0.18]); 
hold on
scatter(f_Resistance,MaxCa_noSOCE6,'filled','MarkerFaceColor',[0.10,0.85,0.83]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Max [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes3,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure3,"Renderer","painters");
%% HIIT Plots -------------------------------------------------------------
figure2 = figure;
axes2 = axes('Parent', figure2);
scatter(f_HIIT,deltaAUC_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes2,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure2,"Renderer","painters");%% 
figure4 = figure;
axes4 = axes('Parent', figure4);
scatter(f_HIIT,deltaAUC_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes4,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure4,"Renderer","painters");

%% Ca
figure3 = figure;
axes3 = axes('Parent', figure3);
plot(Time5{5}(indices1),Ca5{5}(indices1),'Color',[0.64,0.08,0.18]); 
hold on
plot(Time_noSOCE6{5}(indices2),Ca_noSOCE6{5}(indices2),'Color',[0.10,0.85,0.83]); 
xlabel('Time (s)', 'Fontsize',18)
ylabel('[Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes3,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure3,"Renderer","painters");
