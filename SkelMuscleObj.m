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
    saveProgress = false;
elseif length(varargin)==1 %#ok<ISCL>
    Createplot = varargin{1};
    saveProgress = false;
elseif length(varargin)==2
    Createplot = varargin{1};
    progressPath = varargin{2};
    if isfolder(progressPath)
        saveProgress = true;
    end
end

if length(pVec) < 105 || max(pVec) < 1000
    highSensIdx = [2,6,10,14,15,18,20,21,23,24,28,32,33,35,37,40,42,43,45,69,72,77,78,79,80,81,83,86,90,91];
    pRef = ones(105,1);
    pRef(highSensIdx) = pVec;
    pVecStruct = importdata('InputParam1.xlsx');
    p0 = pVecStruct.data;
    pVec = pRef(:) .* p0(:);
end

%Initialize values
InterpExpt = cell(1,9);
Expt_t = cell(1,9);
CompInterp = cell(1,9);
InterpComp_base = cell(1,9);

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
    pVec(74);       % yinit(24) is the initial condition for 'p_i_SR'
    0;          % yinit(25) is the initial condition for 'PiCa'
    pVec(75)       % yinit(26) is the initial condition for 'Pi_Myo'
    1300;        % c_o (µM)
    147000.0;   % Na_o (µM)
    4000.0;      % K_o (µM)
    128000.0;   % Cl_o (µM)
    15000;      % CSQ
    ];
juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
yinit = [yinit(juncLocLogic); yinit(bulkLocLogic)];

% save indices for later
totIdx = sum(juncLocLogic) + sum(bulkLocLogic);
cSRJuncIdx = sum(juncLocLogic(1:2));
cSRBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:2));
ciJuncIdx = sum(juncLocLogic(1:8));
ciBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:8));
forceIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
TTVoltIdx = sum(juncLocLogic(1:5));
SLVoltIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:5));

tSS = 0:1000;
load Exptdata.mat Expt
freq = [100, 100, 67, 67,67,60, 60, 60,60];
T_max = [0.03 0.12 0.06 0.045 0.08 0.025 0.01 0.02 0.002];
expt_title = ["Rincon","Calderon et al. (2010)", "Baylor et al. (2007)",...
    "Hollingworth", "Baylor & Hollingworth", "Yonemura","Bibollet et al. (2023)",...
    "Miranda et al.(2020)","Wallinga"];
expt_n = [3 2 7 8]; % Indices of the experiments used for estimation

%Interpolating experimental values
for m_index = 1 :length(expt_n) %:9
    m = expt_n(m_index);
    T_max(m) = max(Expt{m}(:,1))/1000;
    Expt_t{m} = Expt{m}(:,1)/1000;
    if m < 6 % calcium cases assume baseline 0.1 uM
        Expt{m}(:,2) = Expt{m}(:,2) + 0.1;
    end
    InterpExpt{m} = interp1(Expt_t{m},Expt{m}(:,2),0:0.0001:T_max(m));
end

qoiList = zeros([1,14*4]);
simSaved = cell(8,1);
phosphateAccum = true;

% geom info for averaging
vol_Fiber = pi * (20 ^ 2) * 100 ;
vol_SA_ratio = 0.01;
volFraction_TT = 0.003 ;
vol_myo = 0.95 * vol_Fiber ;
SA_TT = volFraction_TT * vol_Fiber / vol_SA_ratio ;
diffusion_length = pVec(99);
vol_myoJ = SA_TT*diffusion_length;
JFrac = vol_myoJ / vol_myo;
BFrac = 1 - JFrac;

for n_index = 1 :length(expt_n)
    n = -expt_n(n_index); % negative n to indicate this is estimation 
    % (determines EC conc, temperature, and stimulus)
    try
        pVec0 = pVec;
        pVec0(95) = 0; % set SOCE flux to zero
        [~,ySS] = SkelMuscleCa_dydt(tSS,0, yinit, pVec0,tic,n,phosphateAccum);
        if size(ySS,1) < length(tSS) || any(isnan(ySS(:)))
            cSR0 = yinit(cSRBulkIdx);
        else
            cSR0 = ySS(end,cSRBulkIdx);
        end
    catch
        cSR0 = yinit(cSRBulkIdx);
    end
    try
        % pVec(12) is cratio (default value of 0.25)
        pVecCur = pVec;
        pVecCur(12) = pVecCur(12) * cSR0;
        [~,ySS] = SkelMuscleCa_dydt(tSS, 0, yinit, pVecCur, tic, n, phosphateAccum);
        if size(ySS,1) < length(tSS) || any(isnan(ySS(:)))
            yinf = yinit;
            simSaved{abs(n)} = zeros(1,totIdx);
        else
            yinf = ySS(end,:);
        end
    catch
        yinf = yinit;
        fprintf('error in SS computation \n');
    end
    ssQOI = [yinf(cSRBulkIdx), yinf(SLVoltIdx), yinf(sum(juncLocLogic(1:6))),...
             yinf(sum(juncLocLogic(1:7))), yinf(sum(juncLocLogic(1:8))),...
             yinf(sum(juncLocLogic(1:13))), yinf(forceIdx)];

    %% Calculate Dynamics
    t = 0:0.0001:T_max(abs(n));
    % try
        [~,y] = SkelMuscleCa_dydt(t, freq(abs(n)), yinf, pVecCur, tic, n, phosphateAccum);
        if size(y,1) < length(t)
            ySim = y;
            y = zeros(length(t), size(ySim,2));
            y(1:size(ySim,1),:) = ySim;
            for i = 1:size(ySim,2)
                y(size(ySim,1)+1:end,i) = ySim(end,i);
            end
        end

        simSaved{abs(n)} = [t(:),y];
        % Baseline model prediction
        % SS computation
        if Createplot
            pVecStruct = importdata('InputParam1.xlsx');
            p0 = pVecStruct.data;
            pVec0 = p0;
            pVec0(95) = 0; % set SOCE flux to zero
            [~,ySS_baseinit] = SkelMuscleCa_dydt(tSS, 0, yinit, pVec0,tic, n, phosphateAccum);
            cSR0 = ySS_baseinit(end,cSRBulkIdx);
            p0(12) = p0(12) * cSR0;
            [~,ySS_base] = SkelMuscleCa_dydt(tSS, 0, yinit, p0, tic, n, phosphateAccum);
            [~,Y_base] = SkelMuscleCa_dydt(t,freq(abs(n)), ySS_base(end,:), p0, tic, n, phosphateAccum);
        end
        
        CaSol = y(:,ciJuncIdx)*JFrac + y(:,ciBulkIdx)*BFrac;
        VSol = y(:,SLVoltIdx);
        if abs(n) <= 5 % Calcium Calculations
            CompInterp{abs(n)} = CaSol;
            if Createplot
                InterpComp_base{abs(n)} = Y_base(:,ciJuncIdx)*JFrac + Y_base(:,ciBulkIdx)*BFrac;
            end
        elseif abs(n) > 5
            CompInterp{abs(n)} = VSol; % Voltage Calculations
            if Createplot
                InterpComp_base{abs(n)} = Y_base(:,SLVoltIdx);
            end
        end
        
        if any(~isreal(y(:)))
            warning("ode15s returned complex numbers!\n")
            y = real(y);
        end
        MaxCaF = max(CaSol); % Maximum [Ca2+] conc
        MaxVF = max(VSol); % Maximum Voltage
        MaxPost = max(y(:,forceIdx)); % Maximum Force
        AvgF = trapz(t,CaSol) / (t(end)-t(1)); % Area under curve for Calcium
        AvgPost = trapz(t,y(:,forceIdx)) / (t(end)-t(1)); % Area under curve for Force
        AvgVolt = trapz(t,VSol) / (t(end)-t(1)); % Area under curve for Voltage
 
        defaultVal = 0;
        if all(VSol > -60)
            VoltWidth = defaultVal;
        elseif VSol(1) > -60
            VoltWidth = defaultVal;
        else
            low = find(VSol >= -60,1,'first');
            if all(VSol(low:end)>-60)
                VoltWidth = defaultVal;
            else
                high = find(VSol(low:end)<=-60,1,'first');
                high = high + low - 1;
                VoltWidth = interp1(VSol(high-1:high),t(high-1:high),-60) - ...
                            interp1(VSol(low-1:low),t(low-1:low),-60);
            end
        end

        qoiList((n_index-1)*14+1:(n_index*14)) = [ssQOI, MaxCaF, MaxVF, MaxPost, ...
                                                  AvgF, AvgPost, AvgVolt, VoltWidth] ;
    % catch
    %     qoiList((n_index-1)*14+1:(n_index*14)) = [ssQOI, zeros(1,7)];
    %     fprintf('error in dynamics comp \n');
    % 
    % end
end

%% Objective Value Calc
delta = cell(1,9);
sum_delta = zeros(1,9);
Error = cell(1,9);

for j_index = 1 :length(expt_n) %:9
    j = expt_n(j_index);
    weight = length(InterpExpt{j}) ;
    sigma_C = 0.5;
    sigma_V = 5 ;
    delta{j} = CompInterp{j}(:) - InterpExpt{j}(:);
    if j < 6
        Error{j} = ((delta{j} ./ sigma_C ) .^ 2 )./ weight;
    elseif j > 5
        Error{j} = ((delta{j} ./ sigma_V ) .^ 2 )./ weight;
    end
    sum_delta(j) = sum(Error{j});

end

objVal = sum(sum_delta);
if saveProgress
    try
        load(fullfile(progressPath,'objTest.mat'),'objTest')
        if objVal < objTest
            objTest = objVal;
            save(fullfile(progressPath,'objTest.mat'),'objTest')
            save(fullfile(progressPath,'pBest.mat'),'pVec')
        end
    catch
        fprintf('file was busy I guess\n')
    end
end

%% Plots
if Createplot

    figure
    for index = 1:2 % Update according to the number of Calcium expts used
        i = expt_n(index);
        subplot(2,2,index)
        x = 0:0.0001:T_max(i);
        Time_Comp = [x, fliplr(x)];
        Ca_Comp = [InterpExpt{i} + 2*sigma_C, fliplr(InterpExpt{i} - 2*sigma_C)];
        plot(0:0.0001:T_max(i),CompInterp{i},'LineWidth',3, 'color',[0.49,0.18,0.56]) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',3, 'Color',[0.10,0.85,0.83])
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
        plot(0:0.0001:T_max(i),CompInterp{i},'LineWidth',3, 'color',[0.49,0.18,0.56] ) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',3, 'Color',[0.10,0.85,0.83])
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