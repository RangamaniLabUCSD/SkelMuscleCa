function QOI = SkelMuscleCa_SAOutput(param)
% Output Function for Senstivity Analysis 
% Output:
%   QOI are the quantities of interst
% Input:
%   param is the list of variable parameters

input = importdata('InputParam1.xlsx'); % Load default parameters
parameters = input.data; % Parameter values
% load PSO_25-Apr-2024.mat yinit


totSize = size(param,1);
QOI = zeros(size(param,1),14*4+1);
% parpool(50);
for (i = 1 : length(param(:,1)))
    pVec = param(i,:)'.*parameters;
    [objVal, qoiList] = SkelMuscleObj(pVec); % objVal is a scalar, qoiList is a vector of length 14*4
    QOI(i,:) = [qoiList, objVal]; % 14*4 + 1 in length

    % t = 0:1000;
    % tSpan = 0:0.0001:1;
    % freqVec = 50;
    % StartTimer = tic;
    % yinit = [
    %     0.0122; 	% yinit(1) is the initial condition for 'SOCEProb'
    %     1500.0;		% yinit(2) is the initial condition for 'c_SR'
    %     0.9983;		% yinit(3) is the initial condition for 'h_K'
    %     0.9091;		% yinit(4) is the initial condition for 'w_RyR'
    %     -88.0;		% yinit(5) is the initial condition for 'Voltage_PM'
    %     14700.0;	% yinit(6) is the initial condition for 'Na_i'
    %     5830.0;		% yinit(7) is the initial condition for 'Cl_i'
    %     0.1;		% yinit(8) is the initial condition for 'c_i'
    %     0.003;		% yinit(9) is the initial condition for 'n'
    %     0.0128;		% yinit(10) is the initial condition for 'm'
    %     0.8051;		% yinit(11) is the initial condition for 'h'
    %     0.8487;		% yinit(12) is the initial condition for 'S'
    %     154500.0;	% yinit(13) is the initial condition for 'K_i'
    %     0;%387;        % yinit(14) is the initial condition for 'CaParv'
    %     0;%1020;       % yinit(15) is the initial condition for 'MgParv'
    %     0.3632;     % yinit(16) is the inital consition for 'CATP'
    %     0;%10.004;     % yinit(17) is the initial condition for 'CaTrop'
    %     0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    %     0;	    	% yinit(19) is the initial condition for 'D_2'
    %     0;	    	% yinit(20) is the initial condition for 'Pre_Pow'
    %     0;	    	% yinit(21) is the initial condition for 'Post_Pow'
    %     0;	    	% yinit(22) is the initial condition for 'MgATP'
    %     8000;       % yinit(23) is the initial condition for 'ATP'
    %     p(74)       % yinit(24) is the initial condition for 'p_i_SR'
    %     0           % yinit(25) is the initial condition for 'PiCa'
    %     p(75)       % yinit(26) is the initial condition for 'Pi_Myo'
    %     ];
    % 
    % try
    %     % [tSS,ySS,currtimeSS] = SkelMuscleCa_dydt(t,0, false, yinit, p,StartTimer,2); % Steady State
    %     [tSS,ySS,currtimeSS] = SkelMuscleCa_dydtEst(t,0, false, yinit, p,StartTimer,2); 
    %     yInf = ySS(end,:);
    %     ssQOI = [yInf(2), yInf(5), yInf(6), yInf(7),yInf(8), yInf(13), yInf(21)];
    % catch
    %     % yInf = yinit;
    %     % ssQOI = zeros(1,7);
    %     % currtimeSS = toc(StartTimer);
    %     fprintf('error in SS computation for session %d of %d  \n',i,totSize);
    %     continue
    % end
    % 
    % StartTimer = tic;
    % try
    %     % [Time,Y,~] = SkelMuscleCa_dydt(tSpan, freqVec, false, yInf, p,StartTimer,2);  % Dynamics - Peak Ca2+, Peak Voltage, Ca2+ Area under curve
    %     [Time,Y,~] = SkelMuscleCa_dydtEst(tSpan, freqVec, false, yInf, p,StartTimer,2);  % Dynamics - Peak Ca2+, Peak Voltage, Ca2+ Area under curve
    %     for j = 1 : size(Y,1)
    %         for k = 1:size(Y,2)
    %             if ~isreal(Y(j,k))
    %                 Y(j,k) = 0;
    %             end
    %         end
    %     end
    %     MaxCaF = max(Y(:,8));                                            % Maximum [Ca2+] conc for control
    %     MaxVF = max(Y(:,5));                                             % Maximum Voltage for control
    %     MaxPost = max(Y(:,21));                                          % Maximum Force for control
    %     AvgF = trapz(Time,Y(:,8)) / (Time(end)-Time(1));                 % Area under curve for Calcium control
    %     AvgPost = trapz(Time,Y(:,21)) / (Time(end)-Time(1));             % Area under curve for Force control
    %     AvgVolt = trapz(Time,Y(:,5)) / (Time(end)-Time(1));              % Area under curve for Voltage control
    %     %
    %     defaultVal = 0;
    %     if all(Y(:,5) > -60)
    %         VoltWidth = defaultVal;
    %     elseif Y(1,5) > -60
    %         VoltWidth = defaultVal;
    %     else
    %         low = find(Y(:,5)>= -60,1,'first');
    %         if all(Y(low:end,5)>-60)
    %             VoltWidth = defaultVal;
    %         else
    %             high = find(Y(low:end,5)<=-60,1,'first');
    %             high = high + low - 1;
    %             VoltWidth = interp1(Y(high-1:high,5),Time(high-1:high),-60) - interp1(Y(low-1:low,5),Time(low-1:low),-60);
    %         end
    %     end
    %     %
    %     QOI(i,:) = [ssQOI, MaxCaF,MaxVF, MaxPost, AvgF, AvgPost, AvgVolt, VoltWidth]; % Quantities of interest
        fprintf('Session %d of %d \n',i,totSize);
    % 
    % catch
    %     QOI(i,:) = [ssQOI, zeros(1,7)];
    %     fprintf('error in dynamics comp for session %d of %d  \n',i,totSize);
    % 
    % end
end
end
% ssA2, maxA2, AUCa2 (or avg), AUCvoltage/avg(--> width of action potentials)