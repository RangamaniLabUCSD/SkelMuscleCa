clc
clear
delete(gcp('nocreate'))

load PSO_25-Apr-2024.mat p pSol yinit
p = p(:) .* pSol(:);
% Expt : 
% 1,2 - Cont Stim (Control & no SOCE)
% 3,4 - Wei-Lapierre (Control & noSOCE)
% 5,6 - Resistance (Control & noSOCE)
% 7,8 - HIIT (Control & noSOCE)


f_HIIT = 60:10:180;
l_HIIT = length(f_HIIT);
t_HIIT = [0 60];

Time7 = cell(l_HIIT,1);
Ca7 = cell(l_HIIT,1);
Flux7 = cell(l_HIIT,1);
Current7 = cell(l_HIIT,1);
MaxCa7 = zeros(l_HIIT,1);
AUC7 = zeros(l_HIIT,1); 
Force7 = zeros(l_HIIT,1);

Time_noSOCE8 = cell(l_HIIT,1);
Ca_noSOCE8 = cell(l_HIIT,1);
Flux_noSOCE8 = cell(l_HIIT,1);
Current_noSOCE8 = cell(l_HIIT,1);
MaxCa_noSOCE8 = zeros(l_HIIT,1);
AUC_noSOCE8 = zeros(l_HIIT,1);
Force_noSOCE8 = zeros(l_HIIT,1);


f_Resistance = 10:2:50;
l_Resistance = length(f_Resistance);
t_Resistance = [0 420];

Time5 = cell(l_Resistance,1);
Ca5 = cell(l_Resistance,1);
Flux5 = cell(l_Resistance,1);
Current5 = cell(l_Resistance,1);
MaxCa5 = zeros(l_Resistance,1);
AUC5 = zeros(l_Resistance,1);
Force5  = zeros(l_Resistance,1);

Time_noSOCE6 = cell(l_Resistance,1);
Ca_noSOCE6 = cell(l_Resistance,1);
Flux_noSOCE6 = cell(l_Resistance,1);
Current_noSOCE6 = cell(l_Resistance,1);
MaxCa_noSOCE6 = zeros(l_Resistance,1);
AUC_noSOCE6 = zeros(l_Resistance,1);
Force_noSOCE6 = zeros(l_Resistance,1);

% Hill equation variables
n = 3.5947; % Hill coeff
Ca50 = 1.1694; % Calcium conc at 50% max force.
c0 = 1.03167e-04 ; % Vertical shift 
%% Steady State calculation -----------------------------------------------
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 1);                    % Steady State values of variable
yinf = ySS(:,end);
[~,ySS_noSOCE] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 2);  
yinf_noSOCE = ySS_noSOCE(end,:);

%% HIIT -------------------------------------------------------------------
parpool(13)
parfor i = 1:l_HIIT
    freq = f_HIIT(i);

    expt = 7;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt(t_HIIT, freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time7{i} = Time;
    Ca7{i} = Y(:,8);
    Flux7{i} = fluxes;
    Current7{i} = currents;

    expt = 8;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt(t_HIIT, freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE8{i} = Time_noSOCE;
    Ca_noSOCE8{i} = Y_noSOCE(:,8);
    Flux_noSOCE8{i} = fluxes_noSOCE;
    Current_noSOCE8{i} = currents_noSOCE;
    
end

% Max and Avg
for i = 1:l_HIIT
    time_start7 = find(Time7{i} < 1);
    start_index7 = length(time_start7);
    MaxCa7(i) = max(Ca7{i}(start_index7:end));
    AUC7(i) = trapz(Time7{i},Ca7{i})/ max(Time7{i}); 

    time_start8 = find(Time_noSOCE8{i} < 1);
    start_index8 = length(time_start8);
    MaxCa_noSOCE8(i) = max(Ca_noSOCE8{i}(start_index8:end));
    AUC_noSOCE8(i) =  trapz(Time_noSOCE8{i},Ca_noSOCE8{i})/max(Time_noSOCE8{i}); 

end

deltaMax_HIIT = MaxCa7 - MaxCa_noSOCE8;
deltaAUC_HIIT = (AUC7 - AUC_noSOCE8);

% Force
for i = 1 : length(f_HIIT)
    Force7(i) = c0 + (100 * ( AUC7(i)^n )/ ( (Ca50 ^ n)  + (AUC7(i) ^ n)) );
    Force_noSOCE8(i) = c0 + (100 * (AUC_noSOCE8(i)^n )/ ( (Ca50 ^ n)  + (AUC_noSOCE8(i) ^ n)) );  
end

DeltaForce_HIIT = Force7 - Force_noSOCE8;
%% Resistance -------------------------------------------------------------
delete(gcp('nocreate'))
parpool(21)
parfor i = 1:l_Resistance
    freq = f_Resistance(i)

    expt = 5;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt(t_Resistance, freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time5{i} = Time;
    Ca5{i} = Y(:,8);
    Flux5{i} = fluxes;
    Current5{i} = currents;

    expt = 6;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt(t_Resistance, freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE6{i} = Time_noSOCE;
    Ca_noSOCE6{i} = Y_noSOCE(:,8);
    Flux_noSOCE6{i} = fluxes_noSOCE;
    Current_noSOCE6{i} = currents_noSOCE;
    
end

% Max and Avg
for i = 1:l_Resistance

    indices5 = (Time5{i} >=8 & Time5{i} <= 60)| (Time5{i} >= 188 & Time5{i} <= 240) | (Time5{i} >= 368 & Time5{i} <= 420);
    indices6 = (Time_noSOCE6{i} >= 8 & Time_noSOCE6{i} <= 60)| (Time_noSOCE6{i}  >= 188 & Time_noSOCE6{i}  <= 240) | (Time_noSOCE6{i}  >= 368 & Time_noSOCE6{i}  <= 420);
    MaxCa5(i) = max(Ca5{i}(indices5));
    AUC5(i) = trapz(Time5{i}(indices5),Ca5{i}(indices5)) / 180;

    MaxCa_noSOCE6(i) = max(Ca_noSOCE6{i}(indices6));
    AUC_noSOCE6(i) = trapz(Time_noSOCE6{i}(indices6),Ca_noSOCE6{i}(indices6)) / 180;

end
deltaMax_Resistance = MaxCa5 - MaxCa_noSOCE6;
deltaAUC_Resistance = AUC5 - AUC_noSOCE6; 

% Force  
for i = 1 : length(f_Resistance)
    Force5(i) = c0 + (100 * ( AUC5(i)^n )/ ( (Ca50 ^ n)  + (AUC5(i) ^ n)) );
    Force_noSOCE6(i) = c0 + (100 * (AUC_noSOCE6(i)^n )/ ( (Ca50 ^ n)  + (AUC_noSOCE6(i) ^ n)) );
end
DeltaForce_Resistance =  (Force5 - Force_noSOCE6);
