clc
load HIIT_5-9.mat
load HIIT_noSOCE_5-9.mat
load Resistance_5-9.mat
load Resistance_noSOCE_5-9.mat

%% Expt : 
%1,2 - Cont Stim (Control & no SOCE)
%3,4 - Wei-Lapierre (Control & noSOCE)
%5,6 - Resistance (Control & noSOCE)
%%7,8 - HIIT (COntrol & noSOCE)

%% Color code
%[0.89,0.47,0.97] - Pink shade for delta
%[0.64,0.08,0.18] - Maroon shade for control Calcium
%[0.10,0.85,0.83] - Cyan shade for SOCE Blocked calcium
%[0.49,0.18,0.56] -Dark Purple for Voltage

MaxCa5 = zeros(l_Resistance,1);
MaxCa_noSOCE6 = zeros(l_Resistance,1);
MaxCa7 = zeros(l_HIIT,1);
MaxCa_noSOCE8 = zeros(l_HIIT,1);

deltaMaxCa5 = zeros(l_Resistance,1);
deltaMaxCa_noSOCE6 = zeros(l_Resistance,1);
deltaMaxCa7 = zeros(l_HIIT,1);
deltaMaxCa_noSOCE8 = zeros(l_HIIT,1);

AUC5 = zeros(l_Resistance,1);
AUC_noSOCE6 = zeros(l_Resistance,1);
AUC7 = zeros(l_HIIT,1);
AUC_noSOCE8 = zeros(l_HIIT,1);

deltaAUC5 = zeros(l_Resistance,1);
deltaAUC_noSOCE6 = zeros(l_Resistance,1);
deltaAUC7 = zeros(l_HIIT,1);
deltaAUC_noSOCE8 = zeros(l_HIIT,1);

% Max and AUC -------------------------------------------------------------
for i = 1:l_HIIT
    time_start7 = find(Time7{i} < 1);
    start_index7 = length(time_start7);
    MaxCa7(i) = max(Ca7{i}(start_index7:end));
    deltaMaxCa7(i) = MaxCa7(i) - MaxCa7(1);

    AUC7(i) = trapz(Time7{i},Ca7{i}) / 60;
    deltaAUC7(i) = AUC7(i) - AUC(1);

    time_start8 = find(Time_noSOCE8{i} < 1);
    start_index8 = length(time_start8);
    MaxCa_noSOCE8(i) = max(Ca7{i}(start_index7:end));
    deltaMaxCa_noSOCE8(i) = MaxCa_noSOCE8(i) - MaxCa_noSOCE8(1);

    AUC_noSOCE8(i) = trapz(Time7{i},Ca7{i}) / 60;
    deltaAUC_noSOCE8(i) = AUC_noSOCE8(i) - AUC_noSOCE8(1);
end

for i = 1:l_Resistance
    time_start5 = find(Time5{i} < 1);
    start_index5 = length(time_start5);
    MaxCa5(i) = max(Ca5{i}(start_index5:end));
    deltaMaxCa5(i) = MaxCa5(i) - MaxCa5(1);

    AUC7(i) = trapz(Time7{i},Ca7{i}) / 420;
    deltaAUC7(i) = AUC7(i) - AUC(1);

    time_start6 = find(Time_noSOCE6{i} < 1);
    start_index6 = length(time_start6);
    MaxCa_noSOCE6(i) = max(Ca_noSOCE6{i}(start_index6:end));
    deltaMaxCa_noSOCE6(i) = MaxCa_noSOCE6(i) - MaxCa_noSOCE6(1);

    AUC_noSOCE6(i) = trapz(Time7{i},Ca7{i}) / 420;
    deltaAUC_noSOCE6(i) = AUC_noSOCE6(i) - AUC_noSOCE6(1);
end
deltaMax_HIIT = MaxCa7 - MaxCa_noSOCE8;
deltaMax_Resistance = MaxCa5 - MaxCa_noSOCE6;

deltaAUC_HIIT = AUC7 - AUC_noSOCE8;
deltaAUC_Resistance = AUC5 - AUC_noSOCE6; 

figure1 = figure;
axes1 = axes('Parent', figure1);
scatter(f_Resistance,deltaMax_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes1,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure1,"Renderer","painters");

figure2 = figure;
axes2 = axes('Parent', figure2);
scatter(f_Resistance,deltaAUC_HIIT,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes2,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure2,"Renderer","painters");

figure3 = figure;
axes3 = axes('Parent', figure3);
scatter(f_Resistance,deltaMax_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes3,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure3,"Renderer","painters");

figure4 = figure;
axes4 = axes('Parent', figure4);
scatter(f_Resistance,deltaAUC_Resistance,'filled','MarkerFaceColor',[0.89,0.47,0.97]); 
xlabel('Frequency (Hz)', 'Fontsize',18)
ylabel('\Delta Average [Ca^{2+}]_{myo} (μM)','Fontsize',18)
set(axes4,'FontSize',18,'Box','off','FontSmoothing','on')
set(figure4,"Renderer","painters");
