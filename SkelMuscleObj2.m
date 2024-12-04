function [objVal, simSaved] = SkelMuscleObj2(pVec, varargin)
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
    highSensIdx = 1:105;%[12,31,34,39,41,42,43,59:75,89,91,95];%[2,4,5,6,8,10,14,15,16,24,25,26,28,29,32,33,35,37,40,42,43,60,68,74,77,78,80,81,83,90,91,92];
    pRef = ones(105,1);
    pRef(highSensIdx) = pVec;
    paramStruct = importdata('InputParam1.xlsx');
    p0 = paramStruct.data;
    pVec = pRef(:) .* p0(:);
end
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
    param(75);       % yinit(26) is the initial condition for 'Pi_Myo'
    1300;        % c_o (µM)
    147000.0;   % Na_o (µM)
    4000.0;      % K_o (µM)
    128000.0;   % Cl_o (µM)
    15000; % free CSQ
    ];


tSS = 0:1000;
penaltyVal = 100000;
freq = [1;25;50;75;100;125;150;175;200;250];
simSaved = cell(length(freq),1);
maxForce_withSOCE = zeros(size(freq));
forceRatio_withSOCE = zeros(size(freq));
maxForce_noSOCE = zeros(size(freq));
forceRatio_noSOCE = zeros(size(freq));
phosphateAccum = true;

juncLocLogic = true(1,31);
juncLocLogic(17:21) = false; % cross bridges
bulkLocLogic = true(1,31);
bulkLocLogic([1,4,27:30]) = false; % SOCE, wRyR, extracell ions
yinit = [yinit(juncLocLogic); yinit(bulkLocLogic)];

% save indices for later
totIdx = sum(juncLocLogic) + sum(bulkLocLogic);
cSRBulkIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:2));
forceIdx = sum(juncLocLogic) + sum(bulkLocLogic(1:21));
try
    param0 = param;
    param0(95) = 0; % set SOCE flux to zero
    [~,ySS] = SkelMuscleCa_dydt(tSS, 0, yinit, param0,tic,2,false);
    ySS_noSOCE = ySS(end,:);
    if size(ySS,1) < length(tSS) || any(isnan(ySS(:)))
        cSR0 = yinit(cSRBulkIdx);
    else
        cSR0 = ySS(end,cSRBulkIdx);
    end
catch
    cSR0 = yinit(2);
    ySS_noSOCE = yinit;
end
try
    param(12) = param(12) * cSR0;
    [~,ySS] = SkelMuscleCa_dydt(tSS, 0, yinit, param,tic,1,false);
    if size(ySS,1) < length(tSS) || any(isnan(ySS(:)))
        objVal = penaltyVal;
        return
    end
    ySS_withSOCE = ySS(end,:);
catch
    ySS_withSOCE = yinit;
    fprintf('error in SS computation \n');
end

for n = 1:length(freq)
    %% Calculate Dynamics
    t = 0:0.0001:0.5;
    try
        [~,Y_withSOCE] = SkelMuscleCa_dydt(t, freq(n), ySS_withSOCE, param, tic, 1, phosphateAccum); % compute time-dependent solution
        [~,Y_noSOCE] = SkelMuscleCa_dydt(t, freq(n), ySS_noSOCE, param, tic, 2, phosphateAccum); % compute time-dependent solution
        if size(Y_withSOCE,1) < length(t) || size(Y_noSOCE,1) < length(t)
            objVal = penaltyVal;
            simSaved{n} = zeros(1,totIdx);
            return
        end
        simSaved{n} = [t(:),Y_withSOCE];
        A2_withSOCE = Y_withSOCE(:,forceIdx);
        maxForce_withSOCE(n) = max(A2_withSOCE);
        forceRatio_withSOCE(n) = A2_withSOCE(end) / max(A2_withSOCE); 
        A2_noSOCE = Y_noSOCE(:,forceIdx);
        maxForce_noSOCE(n) = max(A2_noSOCE);
        forceRatio_noSOCE(n) = A2_noSOCE(end) / max(A2_noSOCE); 
    catch
        y = ones(length(t), length(ySS_withSOCE));
        y = y .* ySS_withSOCE';
        simSaved{n} = [t(:),y];
        maxForce_withSOCE(n) = ySS_withSOCE(21);
        forceRatio_withSOCE(n) = 1; 
        maxForce_noSOCE(n) = ySS_noSOCE(21);
        forceRatio_noSOCE(n) = 1; 
    end
end
% normalize to max overall force
maxMaxForce = max([maxForce_withSOCE; maxForce_noSOCE]);
maxForce_withSOCE = maxForce_withSOCE / maxMaxForce;
maxForce_noSOCE = maxForce_noSOCE / maxMaxForce;

%% Objective Value Calc
maxForceExpMean_withSOCE = [0.239024373; 0.321951196; 0.595121908; 0.824390184; 0.951219443;...
                            0.999999927; 0.999999927; 0.975609685; 0.946341394; 0.907317007];
maxForceExpSEM_withSOCE = [0.009756097; 0.014634145; 0.019512194; 0.024390242; 0.024390242;...
                           0.024390242; 0.024390242; 0.024390242; 0.029268291; 0.029268291];
forceRatioExpMean_withSOCE = [0.972839506; 0.992592593; 1.002469136; 1.002469136;...
                               0.997530864; 0.982716049; 0.95308642; 0.92345679];
forceRatioExpSEM_withSOCE = [0.009876543; 0.009876543; 0.004938272; 0.004938272;...
                             0.004938272; 0.004938272; 0.009876543; 0.019753086];
maxForceExpMean_noSOCE =  [0.229268276; 0.297560954; 0.546341423; 0.751219457; 0.824390184;...
                           0.848780426; 0.843902377; 0.829268232; 0.80487799; 0.760975554];
maxForceExpSEM_noSOCE = [0.009756097; 0.014634145; 0.024390242; 0.024390242; 0.019512194;...
                         0.019512194; 0.024390242; 0.019512194; 0.019512194; 0.019512194];
forceRatioExpMean_noSOCE = [0.972839506; 0.967901235; 0.958024691; 0.962962963;...
                            0.958024691; 0.928395062; 0.854320988; 0.725925926];
forceRatioExpSEM_noSOCE = [0.004938272; 0.004938272; 0.009876543; 0.009876543;...
                           0.009876543; 0.009876543; 0.014814815; 0.014814815];

objVal = sum(((maxForce_withSOCE - maxForceExpMean_withSOCE) ./ maxForceExpSEM_withSOCE).^2) +...
    sum(((maxForce_noSOCE - maxForceExpMean_noSOCE) ./ maxForceExpSEM_noSOCE).^2) +...
    sum(((forceRatio_withSOCE(3:end) - forceRatioExpMean_withSOCE) ./ forceRatioExpSEM_withSOCE).^2) +...
    sum(((forceRatio_noSOCE(3:end) - forceRatioExpMean_noSOCE) ./ forceRatioExpSEM_noSOCE).^2);

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

if Createplot
    figure
    subplot(2,1,1)
    plot(freq, maxForce_withSOCE, 'b')
    hold on
    plot(freq, maxForce_noSOCE, 'r')
    errorbar(freq, maxForceExpMean_withSOCE, maxForceExpSEM_withSOCE, maxForceExpSEM_withSOCE, 'b')
    errorbar(freq, maxForceExpMean_noSOCE, maxForceExpSEM_noSOCE, maxForceExpSEM_noSOCE, 'r')
    ylabel('Max force')
    xlabel('Freq')
    prettyGraph
    subplot(2,1,2)
    plot(freq(3:end), forceRatio_withSOCE(3:end), 'b')
    hold on
    plot(freq(3:end), forceRatio_noSOCE(3:end), 'r')
    errorbar(freq(3:end), forceRatioExpMean_withSOCE, forceRatioExpSEM_withSOCE, forceRatioExpSEM_withSOCE, 'b')
    errorbar(freq(3:end), forceRatioExpMean_noSOCE, forceRatioExpSEM_noSOCE, forceRatioExpSEM_noSOCE, 'r')
    ylabel('force ratio')
    xlabel('Freq')
    prettyGraph
end

end