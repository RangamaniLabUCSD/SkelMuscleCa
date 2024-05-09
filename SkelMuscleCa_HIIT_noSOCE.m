clc
clear
%delete(gcp('nocreate'))

load PSO_25-Apr-2024.mat pSol yinit p
p = p(:) .* pSol(:);

%% Expt : 
%1,2 - Cont Stim (Control & no SOCE)
%3,4 - Wei-Lapierre (Control & noSOCE)
%5,6 - Resistance (Control & noSOCE)
%%7,8 - HIIT (COntrol & noSOCE)

f_HIIT = 60:10:180;
l_HIIT = length(f_HIIT);

Time_noSOCE8 = cell(l_HIIT,1);
Ca_noSOCE8 = cell(l_HIIT,1);
CaSR_noSOCE8 = cell(l_HIIT,1);
SOCEProb_noSOCE8 = cell(l_HIIT,1);
Flux_noSOCE8 = cell(l_HIIT,1);
Current_noSOCE8 = cell(l_HIIT,1);
%% Steady State calculation -----------------------------------------------
[~,ySS_noSOCE] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 2);  
yinf_noSOCE = ySS_noSOCE(end,:);

parpool(15)
% HIIT Dynamics calculation ---------------------------------------------------
for i = 1:l_HIIT
    freq = f_HIIT(i);    
    expt = 8;
    [Time_noSOCE,Y_noSOCE, ~, fluxes_noSOCE,currents_noSOCE] = SkelMuscleCa_dydt([0 6], freq, 0, yinf_noSOCE, p, tic, expt);   % Dynamics computation
    Time_noSOCE8{i} = Time_noSOCE;
    Ca_noSOCE8{i} = Y_noSOCE(:,8);
    CaSR_noSOCE8{i} = Y_noSOCE(:,2);
    SOCEProb_noSOCE8 = Y_noSOCE(:,1);
    Flux_noSOCE8{i} = fluxes_noSOCE;
    Current_noSOCE8{i} = currents_noSOCE;
end

save('HIIT_noSOCE8_5-9.mat');
