function [objVal, qoiList, simSaved] = SkelMuscleObj(pVec, varargin)
%% Function for calculating the objective value used in estimation or the set of QOIs in sensitivity analysis
% Input:
%   - pVec: parameter vector
% Optional inputs (stored in varargin):
%   - Createplot (varargin{1}): plot the fit for tests considered in
%   objective fcn, deafult is false
%   - progressPath (varargin{2}): if this provides a path to a folder, then
%   set saveProgress to true, and pVec will save to a file 'pBest.mat' in
%   that folder if and only if the current objective value is better than
%   that stored in the file 'objTest.mat' in the specified folder. This
%   allows for incremental saving during a parameter estimation run
% Output:
%   - objVal: value of objective function (see paper for expression)
%   - qoiList: QOIs over each experiment tested. Consists of appended iterations
%     of the following set of QOIs for each experiment:
%       *ssQOI:
%       *MaxCaF:
%       *MaxVF:
%       *MaxPost:
%       *AvgF:
%       *AvgPost:
%       *AvgVolt:
%       *VoltWidth:
%   - simSaved: cell vector containing t and y for each experiment tested
%
% Note that experiments used in estimation are assigned numbers 1-11, but
% only cases 2, 3, and 8 are used in the final estimation. These correspond
% to:
% 2: Rincon et al 2021, Fig 1B calcium data (5 peaks 100 Hz, IIb muscle)
% 3: Baylor and Hollingworth 2003, Fig 2A (fast twitch curve)
% 8: Miranda et al 2020, Fig 2B (WT) membrane voltage data 
% 
% Each of these experiments has associated data stored in
% Data/Exptdata.mat, extracted from original papers using PlotDigitizer

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

VOnly = false;
if length(pVec) < 106 || max(pVec) < 1000
    load Data/p0Struct.mat p0Struct
    p0 = p0Struct.data;
    if length(pVec) == 28 % then fitting to voltage only
        highSensIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
        VOnly = true;
    else % then fitting to both calcium and voltage
        highSensIdx = 1:106;
        VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
        VOnlyStruct = load('Data/pVec_VOnly.mat', 'pVec');
        pVecVOnly = VOnlyStruct.pVec;
        pVecVOnly = pVecVOnly(:); % be sure it is a column vector
        [VOnlyOnly,onlyIdx] = setdiff(VOnlyIdx, highSensIdx); % non overlapping indices
        p0(VOnlyOnly) = pVecVOnly(onlyIdx).*p0(VOnlyOnly); % set p0 according to previous estimation
    end
    pRef = ones(106,1);
    pRef(highSensIdx) = pVec;
    pVec = pRef(:) .* p0(:);
end
pVec(91) = 0;

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
    pVec(75);       % yinit(26) is the initial condition for 'Pi_Myo'
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
load Data/Exptdata.mat Expt
freq = [100, 100, 67, 67,67,60, 60, 60, 60, 67, 15];
T_max = [0.03 0.1 0.05 0.045 0.08 0.025 0.006 0.012 0.002, 0.08, 0.32];
if VOnly
    expt_n = [8];
else
    expt_n = [3,2,8];
end

%Interpolating experimental values
for m_index = 1 :length(expt_n) %:9
    m = expt_n(m_index);
    % T_max(m) = max(Expt{m}(:,1))/1000;
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
vol_SR = 0.05*vol_Fiber;
SRJ_occupancy = 0.5;
SA_SRJ = SA_TT * SRJ_occupancy;
vol_SRJ = SA_SRJ*diffusion_length;
JSRFrac = vol_SRJ / vol_SR;
BSRFrac = 1 - JSRFrac;
ssPenalty = zeros(length(expt_n),1);

for n_index = 1 :length(expt_n)
    n = -expt_n(n_index); % negative n to indicate this is estimation 
    % (determines EC conc, temperature, and stimulus)
    try
        pVec0 = pVec;
        pVec0(95) = 0; % set SOCE flux to zero
        [~,ySS] = SkelMuscleCa_dydt(tSS, 0, yinit, pVec0,tic,n,false);%phosphateAccum);
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
        [~,ySS,~,~,~,ySSFinal] = SkelMuscleCa_dydt(tSS, 0, yinit, pVecCur, tic, n, false);
        if size(ySS,1) < length(tSS) || any(isnan(ySS(:)))
            yinf = yinit;
            simSaved{abs(n)} = zeros(1,totIdx);
        else
            yinf = ySSFinal;
        end
    catch
        yinf = yinit';
        fprintf('error in SS computation \n');
    end
    ssQOI = [yinf(cSRBulkIdx), yinf(SLVoltIdx), yinf(sum(juncLocLogic(1:6))),...
             yinf(sum(juncLocLogic(1:7))), yinf(sum(juncLocLogic(1:8))),...
             yinf(sum(juncLocLogic(1:13))), yinf(forceIdx)];

    %% Calculate Dynamics
    t = 0:0.0001:T_max(abs(n));
    try
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
            load Data/p0Struct.mat p0Struct
            p0 = p0Struct.data;
            pVec0 = p0;
            pVec0(95) = 0; % set SOCE flux to zero
            [~,ySS_baseinit] = SkelMuscleCa_dydt(tSS, 0, yinit, pVec0,tic, n, phosphateAccum);
            cSR0 = ySS_baseinit(end,cSRBulkIdx);
            p0(12) = p0(12) * cSR0;
            [~,ySS_base] = SkelMuscleCa_dydt(tSS, 0, yinit, p0, tic, n, phosphateAccum);
            [~,Y_base] = SkelMuscleCa_dydt(t,freq(abs(n)), ySS_base(end,:), p0, tic, n, phosphateAccum);
        end
        
        CaSol = y(:,ciJuncIdx)*JFrac + y(:,ciBulkIdx)*BFrac;
        CaSRSol = y(:,cSRJuncIdx)*JSRFrac + y(:,cSRBulkIdx)*BSRFrac;
        if ~VOnly
            ssPenalty(n_index) = ((CaSol(1)-0.1)/0.1)^2 + ((CaSRSol(1)-500)/500)^2;
        end
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
    catch
        qoiList((n_index-1)*14+1:(n_index*14)) = [ssQOI, zeros(1,7)];
        fprintf('error in dynamics comp \n');
        if abs(n) <= 5 % Calcium Calculations
            CompInterp{abs(n)} = ssQOI(5)*ones(size(t));
        elseif abs(n) > 5
            CompInterp{abs(n)} = ssQOI(2)*ones(size(t));
        end
    end
end

%% Objective Value Calc
delta = cell(1,9);
sum_delta = zeros(1,9);
Error = cell(1,9);

for j_index = 1 :length(expt_n) %:9
    j = expt_n(j_index);
    weight = length(InterpExpt{j}) ;
    sigma_C = 0.5;
    sigma_V = 5;
    delta{j} = CompInterp{j}(:) - InterpExpt{j}(:);
    if j < 6
        Error{j} = ((delta{j} ./ sigma_C ) .^ 2 )./ weight;
    elseif j > 5
        Error{j} = ((delta{j} ./ sigma_V ) .^ 2 )./ weight;
    end
    sum_delta(j) = sum(Error{j});

end

objVal = sum(sum_delta);
objVal = objVal + sum(ssPenalty);
if saveProgress
    try
        load(fullfile(progressPath,'objTest.mat'),'objTest')
        if objVal < objTest
            objTest = objVal;
            save(fullfile(progressPath,'objTest.mat'),'objTest')
            save(fullfile(progressPath,'pBest.mat'),'pVec')
        end
    catch
        fprintf('file was busy most likely\n')
    end
end

%% Plots
if Createplot
    figure
    if VOnly
        for index = 1
            i = expt_n(index); %% Update according to the number of Voltage expts used
            subplot(1,2,index)
            x = 0:0.0001:T_max(i);
            Time_Comp = [x, fliplr(x)];
            V_Comp = [InterpExpt{i} + 2*sigma_V, fliplr(InterpExpt{i} - 2*sigma_V)];
            plot(0:0.0001:T_max(i),CompInterp{i},'LineWidth',2, 'color',[0.49,0.18,0.56] ) %'b',
            hold on
            plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',2, 'Color',[0.10,0.85,0.83])
            plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineStyle', '--' ,'LineWidth',2)
            fill(Time_Comp,V_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)       %'color',[0.00,0.00,1.00]
            xlabel('Time (s)');
            ylabel('Membrane Potential (mV)');
            % title('V_{SL} for expt - '+ expt_title(i));
            prettyGraph
        end
        legend('Calibrated','Pre-Calibrated','Experiment', '95% Confidence Interval')
    else
        for index = 1:2 % Update according to the number of Calcium expts used
            i = expt_n(index);
            subplot(2,2,index)
            x = 0:0.0001:T_max(i);
            Time_Comp = [x, fliplr(x)];
            Ca_Comp = [InterpExpt{i} + 2*sigma_C, fliplr(InterpExpt{i} - 2*sigma_C)];
            plot(0:0.0001:T_max(i),CompInterp{i},'LineWidth',2, 'color',[0.49,0.18,0.56]) %'b',
            hold on
            plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',2, 'Color',[0.10,0.85,0.83])
            plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2,'Linestyle','--')
            fill(Time_Comp,Ca_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)
            xlabel('Time (s)');
            % title('[Ca^{2+}] for '+ expt_title(i));
            ylabel('Concentration (μM)');
            prettyGraph
        end
        for index = 3
            i = expt_n(index); %% Update according to the number of Voltage expts used
            subplot(2,2,index)
            x = 0:0.0001:T_max(i);
            Time_Comp = [x, fliplr(x)];
            V_Comp = [InterpExpt{i} + 2*sigma_V, fliplr(InterpExpt{i} - 2*sigma_V)];
            plot(0:0.0001:T_max(i),CompInterp{i},'LineWidth',2, 'color',[0.49,0.18,0.56] ) %'b',
            hold on
            plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',2, 'Color',[0.10,0.85,0.83])
            plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineStyle', '--' ,'LineWidth',2)
            fill(Time_Comp,V_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)       %'color',[0.00,0.00,1.00]
            xlabel('Time (s)');
            ylabel('Membrane Potential (mV)');
            % title('V_{SL} for expt - '+ expt_title(i));
            prettyGraph
        end
        legend('Calibrated','Pre-Calibrated','Experiment', '95% Confidence Interval')
    end
end
end