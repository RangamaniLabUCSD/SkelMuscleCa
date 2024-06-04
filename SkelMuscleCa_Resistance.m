clc
clear
delete(gcp('nocreate'))
load PSO_25-Apr-2024.mat pSol yinit p

p = p(:) .* pSol(:);

% Expt : 
%  1,2 - Cont Stim (Control & no SOCE)
%  3,4 - Wei-Lapierre (Control & noSOCE)
%  5,6 - Resistance (Control & noSOCE)
%  7,8 - HIIT (COntrol & noSOCE)

f_Resistance = 10:30;

l_Resistance = length(f_Resistance);
Time5 = cell(l_Resistance,1);
Ca5 = cell(l_Resistance,1);
CaSR5 = cell(l_Resistance,1);
SOCEProb5 = cell(l_Resistance,1);
Flux5 = cell(l_Resistance,1);
Current5 = cell(l_Resistance,1);

Time_noSOCE6 = cell(l_Resistance,1);
Ca_noSOCE6 = cell(l_Resistance,1);
CaSR_noSOCE6 = cell(l_Resistance,1);
SOCEProb_noSOCE6 = cell(l_Resistance,1);
Flux_noSOCE6 = cell(l_Resistance,1);
Current_noSOCE6 = cell(l_Resistance,1);

%% Steady State calculation -----------------------------------------------
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 1);                    % Steady State values of variable
yinf = ySS(end,:);

[~,ySS_noSOCE] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 2);  
yinf_noSOCE = ySS_noSOCE(end,:);

%% Resistance Dynamics calculation ---------------------------------------------------
parpool(10)
parfor i = 1:l_Resistance
    freq = f_Resistance(i)

    expt = 5;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt([0 420], freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time5{i} = Time;
    Ca5{i} = Y(:,8);
    CaSR5{i} = Y(:,2);
    SOCEProb5{i} = Y(:,1);
    Flux5{i} = fluxes;
    Current5{i} = currents;

    expt = 6;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt([0 420], freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE6{i} = Time_noSOCE;
    Ca_noSOCE6{i} = Y_noSOCE(:,8);
    CaSR_noSOCE6{i} = Y_noSOCE(:,2);
    SOCEProb_noSOCE6 = Y_noSOCE(:,1);
    Flux_noSOCE6{i} = fluxes_noSOCE;
    Current_noSOCE6{i} = currents_noSOCE;
    
end

save('Resistance5-21.mat');
