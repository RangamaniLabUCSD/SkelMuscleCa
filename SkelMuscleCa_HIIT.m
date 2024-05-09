clc
clear
delete(gcp('nocreate'))

load PSO_25-Apr-2024.mat pSol yinit p
p = p(:) .* pSol(:);

%% Expt : 
%1,2 - Cont Stim (Control & no SOCE)
%3,4 - Wei-Lapierre (Control & noSOCE)
%5,6 - Resistance (Control & noSOCE)
%%7,8 - HIIT (COntrol & noSOCE)

f_HIIT = 60:10:180;
l_HIIT = length(f_HIIT);
Time7 = cell(l_HIIT,1);
Ca7 = cell(l_HIIT,1);
CaSR7 = cell(l_HIIT,1);
SOCEProb7 = cell(l_HIIT,1);
Flux7 = cell(l_HIIT,1);
Current7 = cell(l_HIIT,1);

%% Steady State calculation -----------------------------------------------
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p, tic, 1);                    % Steady State values of variable
yinf = ySS(:,end);

parpool(15)
% HIIT Dynamics calculation ---------------------------------------------------
parfor i = 1:l_HIIT
    freq = f_HIIT(i);

    expt = 7;
    [Time,Y, ~, fluxes,currents] = SkelMuscleCa_dydt([0 60], freq, 0, yinf, p, tic, expt);   % Dynamics computation
    Time7{i} = Time;
    Ca7{i} = Y(:,8);
    CaSR7{i} = Y(:,2);
    SOCEProb7{i} = Y(:,1);
    Flux7{i} = fluxes;
    Current7{i} = currents;
    
end

save('HIIT7_5-9.mat');


