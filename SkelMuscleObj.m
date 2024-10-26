function [objVal, qoiList, simSaved] = SkelMuscleObj(pVec, varargin)
%% Function for calculating the objective value for estimation
% Input:
%       pVec - Vector of particles
%       CreatePlot (optional) - logical variable for whether to plot
%       outputs
% Output:
%       objVal - Minimum error between the experimental and model
%       output
if isempty(varargin)
    Createplot = false;
else
    Createplot = varargin{1};
end

% Createplot = varargin{1};
if length(pVec) < 97 || max(pVec) < 1000
    highSensIdx = [2,4,5,6,8,10,14,15,16,24,25,26,28,29,32,33,35,37,40,42,43,60,68,74,77,78,80,81,83,90,91,92];
    pRef = ones(97,1);
    pRef(highSensIdx) = pVec;
    paramStruct = importdata('InputParam1.xlsx');
    p0 = paramStruct.data;
    pVec = pRef(:) .* p0(:);
end

%Initialize values
T_max = zeros(1,9);
InterpExpt = cell(1,9);
Expt_t = cell(1,9);
InterpComp = cell(1,9);
InterpComp_base = cell(1,9);
CompV = cell(1,5);
CompC = cell(1,5);
StartTimer = tic;
param = pVec;

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
    0;%387;        % yinit(14) is the initial condition for 'CaParv'
    0;%1020;       % yinit(15) is the initial condition for 'MgParv'
    0.3632;     % yinit(16) is the inital consition for 'CATP'
    0;%10.004;     % yinit(17) is the initial condition for 'CaTrop'
    0;	    	% yinit(18) is the initial condition for 'CaCaTrop'
    0;	    	% yinit(19) is the initial condition for 'D_2'
    0;	    	% yinit(20) is the initial condition for 'Pre_Pow'
    0;	    	% yinit(21) is the initial condition for 'Post_Pow'
    0;	    	% yinit(22) is the initial condition for 'MgATP'
    8000;       % yinit(23) is the initial condition for 'ATP'
    param(74)       % yinit(24) is the initial condition for 'p_i_SR'
    0           % yinit(25) is the initial condition for 'PiCa'
    param(75)       % yinit(26) is the initial condition for 'Pi_Myo'
    ];


tSS = 0:1000;
penaltyVal = 100000;

%Experimental Data
%Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C], [K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_T MJ_V]};
load Exptdata.mat Expt
freq = [100, 100, 67, 67,67,60, 60, 60,60];
T_max = [0.03 0.12 0.06 0.045 0.08 0.025 0.01 0.02 0.002];
expt_title = ["Rincon","Calderon et al. (2010)", "Baylor et al. (2007)", "Hollingworth", "Baylor & Hollingworth", "Yonemura","Bibollet et al. (2023)", "Miranda et al.(2020)","Wallinga"];
expt_n = [3 2 7 8]; % Index of the experimental used for estimation

%Interpolating experimental values
for m_index = 1 :length(expt_n) %:9
    m = expt_n(m_index);
    T_max(m) = max(Expt{m}(:,1))/1000;
    Expt_t{m} = Expt{m}(:,1)/1000;
    if m < 6
        Expt{m}(:,2) = Expt{m}(:,2) + 0.1;
    end
    InterpExpt{m} = interp1(Expt_t{m},Expt{m}(:,2),0:0.0001:T_max(m));
end

qoiList = zeros([1,14*4]);
badSRCa = false;
simSaved = cell(8,1);

for n_index = 1 :length(expt_n)
    n = expt_n(n_index);
    try
        param0 = param;
        param0(95) = 0; % set SOCE flux to zero
        [~,ySS] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, param0,StartTimer,n);
        if size(ySS,1) < length(tSS) || any(isnan(ySS))
            cSR0 = yinit(2);
        else
            cSR0 = ySS(end,2);
        end
    catch
        cSR0 = yinit(2);
    end
    StartTimer = tic;
    try
        % error("fake error")
        % param(12) is cratio (default value of 0.25)
        param(12) = param(12) * cSR0;
        [TimeSS,ySS] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, param,StartTimer,n);
        if size(ySS,1) < length(tSS)
            objVal = penaltyVal;
            simSaved{n} = zeros(1,29);
            return
        end
        if any(isnan(ySS))
            objVal = penaltyVal;
            simSaved{n} = zeros(1,29);
            return
        end
        % if ySS(end,2) < 800 || ySS(end,2) > 2000 || param(12) > 0.5*ySS(end,2) %changed from 0.4
        %     badSRCa = true;
        % end
        yinf = ySS(end,:);
        ssQOI = [yinf(2), yinf(5), yinf(6), yinf(7),yinf(8), yinf(13), yinf(21)];
    catch
        yinf = yinit;
        fprintf('error in SS computation \n');
        ssQOI = [yinf(2), yinf(5), yinf(6), yinf(7),yinf(8), yinf(13), yinf(21)];
        % continue
    end


    %% Calculate Dynamics
    t = 0:0.0001:T_max(n);
    StartTimer = tic;
    try
        [~,y] = SkelMuscleCa_dydtEst(t,freq(n), 0, yinf, param, StartTimer,n);
        if size(y,1) < length(t)
            objVal = penaltyVal;
            simSaved{n} = zeros(1,29);
            return
        end

        simSaved{n} = [t(:),y];
        % Baseline model prediction
        % SS computation
        if Createplot
            paramStruct = importdata('InputParam1.xlsx');
            p0 = paramStruct.data;
            param0 = p0;
            param0(95) = 0; % set SOCE flux to zero
            [~,ySS_baseinit] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, param0,StartTimer,n);
            cSR0 = ySS_baseinit(end,2);
            p0(12) = p0(12) * cSR0;
            [TimeSS_base,ySS_base] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, p0,StartTimer,n);
            [~,Y_base] = SkelMuscleCa_dydtEst(t,freq(n), 0, ySS_base(end,:), p0, StartTimer,n);
        end

        if n <= 5
            InterpComp{n} = y(:,8) ; % Calcium Calculations
            CompV{n} = y(:,5);
            if Createplot
                InterpComp_base{n} = Y_base(:,8);
            end
        elseif n > 5
            InterpComp{n} = y(:,5); % Voltage Calculations
            CompC{n} = y(:,8);
            if Createplot
                InterpComp_base{n} = Y_base(:,5);
            end
        end

        Y= y;

        for j = 1 : size(Y,1)
            for k = 1:size(Y,2)
                if ~isreal(Y(j,k))
                    Y(j,k) = 0;
                end
            end
        end
        MaxCaF = max(Y(:,8));                                            % Maximum [Ca2+] conc for control
        MaxVF = max(Y(:,5));                                             % Maximum Voltage for control
        MaxPost = max(Y(:,21));                                          % Maximum Force for control
        AvgF = trapz(t,Y(:,8)) / (t(end)-t(1));                 % Area under curve for Calcium control
        AvgPost = trapz(t,Y(:,21)) / (t(end)-t(1));             % Area under curve for Force control
        AvgVolt = trapz(t,Y(:,5)) / (t(end)-t(1));              % Area under curve for Voltage control
 
        defaultVal = 0;
        if all(Y(:,5) > -60)
            VoltWidth = defaultVal;
        elseif Y(1,5) > -60
            VoltWidth = defaultVal;
        else
            low = find(Y(:,5)>= -60,1,'first');
            if all(Y(low:end,5)>-60)
                VoltWidth = defaultVal;
            else
                high = find(Y(low:end,5)<=-60,1,'first');
                high = high + low - 1;
                VoltWidth = interp1(Y(high-1:high,5),t(high-1:high),-60) - interp1(Y(low-1:low,5),t(low-1:low),-60);
            end
        end

        qoiList((n_index-1)*14+1:(n_index*14)) = [ssQOI, MaxCaF, MaxVF, MaxPost, AvgF, AvgPost, AvgVolt, VoltWidth] ;
        % fprintf('Session %d of %d , time = %.2f s \n',i,totSize,currtimeSS);

    catch
        qoiList((n_index-1)*14+1:(n_index*14)) = [ssQOI, zeros(1,7)];
        fprintf('error in dynamics comp \n');

    end
end

%% Objective Value Calc
% if badSRCa
%     objVal = penaltyVal;
% else
delta = cell(1,9);
sum_delta = zeros(1,9);
Error = cell(1,9);

for j_index = 1 :length(expt_n) %:9
    j = expt_n(j_index);
    weight = length(InterpExpt{j}) ;
    sigma_C = 0.5;
    sigma_V = 5 ;
    delta{j} = InterpComp{j}' - InterpExpt{j};
    if j < 6
        Error{j} = ((delta{j} ./ sigma_C ) .^ 2 )./ weight;
    elseif j > 5
        Error{j} = ((delta{j} ./ sigma_V ) .^ 2 )./ weight;
    end
    sum_delta(j) = sum(Error{j});

end

objVal = sum(sum_delta);


%% Plots
if Createplot

    figure
    for index = 1:2 % Update according to the number of Calcium expts used
        i = expt_n(index);
        subplot(2,2,index)
        x = 0:0.0001:T_max(i);
        Time_Comp = [x, fliplr(x)];
        Ca_Comp = [InterpExpt{i} + 2*sigma_C, fliplr(InterpExpt{i} - 2*sigma_C)];
        plot(0:0.0001:T_max(i),InterpComp{i},'LineWidth',3, 'color',[0.49,0.18,0.56]) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',3, 'Color',[0.10,0.85,0.83])
        hold on
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',3,'Linestyle','--')
        fill(Time_Comp,Ca_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)
        xlabel('Time (s)', 'FontSize',18);
        title('[Ca^{2+}] for '+ expt_title(i), 'FontSize',18);
        ylabel('Concentration (uM)', 'FontSize',18);
        ax = gca;
        ax.FontSize = 18;


    end

    for index = 3:4
        i = expt_n(index); %% Update according to the number of Voltage expts used
        subplot(2,2,index)
        x = 0:0.0001:T_max(i);
        Time_Comp = [x, fliplr(x)];
        V_Comp = [InterpExpt{i} + 2*sigma_V, fliplr(InterpExpt{i} - 2*sigma_V)];
        plot(0:0.0001:T_max(i),InterpComp{i},'LineWidth',3, 'color',[0.49,0.18,0.56] ) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',3, 'Color',[0.10,0.85,0.83])
        hold on
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineStyle', '--' ,'LineWidth',3)
        fill(Time_Comp,V_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)       %'color',[0.00,0.00,1.00]
        xlabel('Time (s)', 'FontSize',18);
        ylabel('Membrane Potential (mV)', 'FontSize',18);
        title('V_{SL} for expt - '+ expt_title(i), 'FontSize',18);
        ax = gca;
        ax.FontSize = 18;

    end
    legend('Calibrated','Pre-Calibrated','Experiment', '95% Confidence Interval', 'FontSize',18)
end
end