function QOI = SkelMuscleCa_SAOutput(param)
% Output Function for Senstivity Analysis
% Output:
%   QOI are the quantities of interst 
% Input:
%   param is the list of variable parameters

input = importdata('InputParam1.xlsx'); % Load default parameters
parameters = input.data; % Parameter values
% load PSO_25-Apr-2024.mat yinit 
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
    10          % yinit(26) is the initial condition for 'p_i_SR'
    0.001       % yinit(27) is the initial condition for 'PiCa'
    0.001       % yinit(28) is the initial condition for 'Pi_Myo'
    ];

QOI = zeros(size(param,1),9);
parpool(50);

   parfor (i = 1 : length(param(:,1)))
        p = param(i,:)'.*parameters;
        t = [0 1000];
        tSpan = [0 1];
        freqVec = 100;
        StartTimer = tic;
    
         try
            [tSS,ySS,currtimeSS] = SkelMuscleCa_dydt(t,0, false, yinit, p,StartTimer,2); % Steady State
            yInf = ySS(end,:);               
   
            [Time,Y,~] = SkelMuscleCa_dydt(tSpan, freqVec, false, yInf, p,StartTimer,2);  % Dynamics - Peak Ca2+, Peak Voltage, Ca2+ Area under curve
            for j = 1 : size(Y,1)
                for k = 1:size(Y,2)
                    if ~isreal(Y(j,k))
                        Y(j,k) = 0;
                    end
                end
            end
            MaxCaF = max(Y(:,8));                                            % Maximum [Ca2+] conc for control
            MaxVF = max(Y(:,5));                                             % Maximum Voltage for control
            AreaF = trapz(Time,Y(:,8));                                       % Area under curve for control
            QOI(i,:) = [yInf(2), yInf(5), yInf(6), yInf(7),yInf(8), yInf(13),MaxCaF,MaxVF,AreaF]; % Quantitites of interest         
            fprintf('Session %d of %d , time = %.2f s \n',i,size(param,1),currtimeSS);
    
        catch
            yInf = yinit;
            [Time,Y,~] = SkelMuscleCa_dydt(tSpan, freqVec, false, yInf, p,tic,2);         % Control.
            for j = 1 : size(Y,1)
                for k = 1:size(Y,2)
                    if ~isreal(Y(j,k))
                        Y(j,k) = 0;
                    end
                end
            end
            MaxCaF = max(Y(:,8));                                            % Maximum [Ca2+] conc for control
            MaxVF = max(Y(:,5));                                             % Maximum Voltage for control
            AreaF = trapz(Time,Y(:,8));                                       % Area under curve for control
            QOI(i,:) = [0,0,0,0,0,0,MaxCaF,MaxVF,AreaF];                       % Quantitites of interest
            fprintf('Session %d of %d, with yinit \n',i,size(param,1));

        end
    end
end
