%% particle swarm
load PSOMito_02-Jul-2026.mat pSol
mitoParams = pSol;
% mitoParams(40) = 1;
logLikelihoodMito(mitoParams)
try
    psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
        'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 20, 'SwarmSize', 30,...
        'InitialPoints', mitoParams);
catch
    psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
        'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 20, 'SwarmSize', 30,...
        'InitialSwarmMatrix', mitoParams);
end
objVal = inf;
mitoParams = ones(40,1);
save('bestMitoParam.mat','objVal', 'mitoParams')

delete(gcp('nocreate'))
if psOptions.UseParallel
    parpool(30) %% **CHANGE SWARMSIZE!** and save file location
end
username = getenv('USER');
progressPath1 = sprintf('/tscc/lustre/ddn/scratch/%s/SkelMuscleCaData',username);
progressPath2 = sprintf('/expanse/lustre/scratch/%s/temp_project/SkelMuscleCaData',username);
if isfolder(progressPath1)
    progressPath = progressPath1;
elseif isfolder(progressPath2)
    progressPath = progressPath2;
else
    progressPath = pwd;
end
% lb = 0.2*ones(7,1); % PP, V_Pi, g0_anion, V_anion, KPi_i, Ap, Bp
% ub = 5.0*ones(7,1);
% lb(4) = 0.5;
% ub(4) = 2;
% lb(3) = 0;
% lb(1) = 0.04;
% ub(1) = 25;0
lb = 0.5*ones(40,1);
lb(38) = max([lb(38),1.01*(10000/15000)]); % do not violate bounds for ATP/ADP positivity
lb(39) = max([lb(39),1.01*(50/250)]); % do not violate bounds for NAD/NADH positivity
ub = 2.0*ones(40,1);
lb(17) = 167/177;%177 
ub(17) = 187/177;
lb(23) = 180/190;%190
ub(23) = 200/190;

[pSol,fval,exitflag] = particleswarm(@logLikelihoodMito,length(lb),lb,ub,psOptions);
filename = "PSOMito_" + string(datetime("today")) +".mat";
save(fullfile(progressPath,filename));

%% Load data
file_names = {'Pi0.001.csv','Pi0.01.csv','Pi0.1.csv','Pi1.csv'};
folder_names = {'csv_ca_out', 'csv_mito_ca', 'csv_bound_free', 'csv_voltage'};
ylabels = {'Cyto Ca2+ (µM)', 'Mito Ca2+ (µM)', 'Bound/free Ca2+', 'Mito voltage (mV)'};
errVals = [0.25, 5, 100, 10];
CaOutData =cell(1,4);
figure
yMins = [0,0,0,150];
yMaxes = [6,60,850,230];
defaultColors = colororder;
for i = 1:length(folder_names)
    for j = 1:length(file_names)
        % subplot(4,4,4*(i-1)+j)
        subplot(1,4,i)
        xlim([50 350])
        ylim([yMins(i),yMaxes(i)])
        hold on
        curFile = sprintf('Data/%s/%s', folder_names{i}, file_names{j});
        curDataCell = readcell(curFile);
        curData = cell2mat(curDataCell(2:end,:));
        plot(curData(:,1),curData(:,2))
        tLoop = [curData(:,1); curData(end:-1:1,1)];
        yLoop = [curData(:,2) + errVals(i); curData(end:-1:1,2) - errVals(i)];
        fill(tLoop, yLoop, defaultColors(j,:), 'LineStyle', 'none', 'FaceAlpha', 0.2)
        if i == 1
            CaOutData{j} = curData;
        end
        prettyGraph
        ylabel(ylabels{i})
        xlabel('Time (s)')
    end
end
legend('Pi = 1 µM', 'Pi = 10 µM', 'Pi = 100 µM', 'Pi = 1000 µM')

%% plot fit
load PSOMito_11-Jul-2026.mat pSol
pMito = pSol;
% load pSol_fullWithMito/pBestMito_Jul11.mat mitoParams
% pMito = mitoParams;
% pMito(40) = 2.0;
% first test driven ca case
tTest = 0:400.0;
freqs = [0,0.1,1,10,20];%,10,20];%logspace(-1,2,20);
phosVals = [0,1,10,100,1000];
figure
% yinit = [0.1; 10000; 1000; 50; 150]; % Pi might be closer to 5000 in mito
yinit = [0.1; 10000; 0; 50; 150; 0.0]; % Pi might be closer to 5000 in mito
yinit(7) = pSol(38)*15000 - yinit(2);
yinit(8) = pSol(39)*250 - yinit(4);
% initOpts = optimoptions("fsolve", 'MaxFunctionEvaluations', 10000);
% ySSfun = @(y0) (mito_dyn(0.0, y0, [tTest', 0.1*ones(size(tTest))', 5000*ones(size(tTest))', 100*ones(size(tTest))', 1000.0*0.1*ones(size(tTest))']));
% yinit = fsolve(ySSfun, yinit);
options = odeset('RelTol',1e-3,'NonNegative',[1:4,6:8],'MaxStep',0.01);
for j = 1:length(phosVals)
    % freq = freqs(j);
    PiCur = phosVals(j);
    if PiCur == 0
        tStarts = [];
        c_i = zeros(size(tTest))' + CaOutData{j}(1,2);
        tVals = tTest';
    else
        % tStarts = tTest(1):(1/freq):tTest(end);
        % c_i = 0.01+7*(heaviside(tTest-120));%-heaviside(tTest-15));
        cData = CaOutData{j-1};
        % c_i = cData(:,2);
        % tVals = cData(:,1);
        c_i = interp1(cData(:,1),cData(:,2),tTest');
        tVals = tTest';
        % c_i = smooth(c_i, 1001);
    end
    
    % c_i = zeros(size(tTest)) + 0.1;
    % cMag = 5.0;
    % for i = 1:length(tStarts)
    %     [~,startIdx] = min((tTest-tStarts(i)).^2);
    %     c_i(startIdx:end) = c_i(startIdx:end) + cMag*exp(-(tTest(startIdx:end)-tStarts(i))/0.02);
    % end
    % c_i(tTest >= 20) = 0.1;
    
    ATPvals = 0.01*ones(size(c_i)); %5000
    ADPvals = 0.01*ones(size(c_i)); %100
    Pivals = PiCur*ones(size(c_i));
    Pivals(tVals < 240) = 0;%
    cytoData = [tVals, c_i, ATPvals, ADPvals, Pivals];
    pCur = ones(40,1);
    pCur(:) = pMito;
    % pCur(32) = 1.0;
    % varyIdx = [1,32,36,37,38,39,40]; % PP, V_Pi, g0_anion, V_anion, KPi_i, Ap, Bp
    % pCur(36)= 10;
    % pCur(varyIdx) = pSol;
    % pCur(42) = 0.01;
    [t,y] = ode15s(@mito_dyn, tTest, yinit, options, cytoData, pCur, tic);

    if PiCur == 0 & j==1
        yinit = y(end,:);
        continue
    else
        % subplot(4,4,j-1)
        subplot(1,4,1)
        plot(tVals, c_i)
        hold on
        % xlim([tTest(1),tTest(end)])
        xlim([50 350])
        ylim([yMins(1), yMaxes(1)])
        ylabel(ylabels{1})
        xlabel('Time (s)')
    end
    fluxes = zeros(length(t),10);
    for i = 1:length(t)
        [~,fluxes(i,:)] = mito_dyn(t(i),y(i,:),cytoData,pCur,tic);
    end
    prettyGraph
    % subplot(4,4,4+j-1)
    subplot(1,4,2)
    plot(t, y(:,1))
    hold on
    prettyGraph
    xlim([50 350])
    ylim([yMins(2), yMaxes(2)])
    ylabel(ylabels{2})
    xlabel('Time (s)')
    % subplot(4,4,8+j-1)
    subplot(1,4,3)
    % V_F1FO = 35000; % uM/s
    % p6 = 200; %mV
    % p7 = 8.5; %mV
    % K_iATP = 10000; %uM
    % J_F1FO = V_F1FO*(1./(1+exp((p6-y(:,5))/p7))).*(K_iATP./(K_iATP+y(:,2)));
    % plot(t, y(:,3))
    bM = pCur(2)*0.01;
    cratio = (y(:,1)*(1-bM) + y(:,6)*bM)./(y(:,1)*bM);
    plot(t, cratio)%1/(pCur(2)*0.01) + y(:,6)./y(:,1))%cumtrapz(t,fluxes(:,1)-fluxes(:,2)+fluxes(:,3)))
    hold on
    prettyGraph
    xlim([50 350])
    ylim([yMins(3), yMaxes(3)])
    ylabel(ylabels{3})
    xlabel('Time (s)')
    % subplot(4,4,12+j-1)
    subplot(1,4,4)
    plot(t, y(:,5))
    hold on
    drawnow
    prettyGraph
    xlim([50 350])
    ylim([yMins(4), yMaxes(4)])
    ylabel(ylabels{4})
    xlabel('Time (s)')
    % figure
    % a1 = 20;
    % a2 = 3.43;
    % a3 = 0.5;
    % plot(t,a1*fluxes(:,5) - a2*fluxes(:,7) + a3*fluxes(:,10) - fluxes(:,9) - fluxes(:,2) - 2*fluxes(:,1) - 2*fluxes(:,3) - fluxes(:,8) - fluxes(:,6))
end

%%
logLikelihoodMito(pMito)

%% plot ANT
ATP_totvec = 1000:2000:15000;%15000];%1:10000;
V_Mvec = 0:10:200;
ATP_M = 1:15000;%[1000,2000,4000,8000,15000];%5000;
figure
hold on
for i = 1:length(ATP_totvec)
    ATP_tot = ATP_totvec(i);
    ADP = 15000 - ATP_tot;
    % ATP_M = ATP_Mvec(i);
    ADP_M = 15000 - ATP_M;
    alpha_c = 0.111;
    alpha_m = 0.139;
    F = 96485.3321;
    V_M = 150;
    R = 8314.46261815;
    T = 273.15 + 22;
    V_ANT = 5000;
    f_ANT = 0.5;
    R_c = alpha_c*ATP_tot./ADP;
    R_m = ADP_M./(alpha_m*ATP_M);
    ANT_num = (1-R_c*R_m*exp(-F*V_M/(R*T)));
    ANT_denom = (1+(alpha_c*ATP_tot/ADP)*exp(-f_ANT*F*V_M/(R*T)))*(1+ADP_M./(ATP_M*alpha_m));
    J_ANT = V_ANT.*ANT_num./ANT_denom;
    plot(ATP_M, J_ANT)
end

function [dydt,fluxes] = mito_dyn(t, y, cytoData, p, timer)
    % [~,idx] = min((cytoData(:,1)-t).^2);
    if toc(timer) > 20
        error('Took too long!!')
    end
    idx = find(cytoData(:,1) > t, 1, 'first');
    if isempty(idx) || idx <= 1
        if t==cytoData(end,1)
            idx = length(cytoData(:,1)) - 1;
            tInc = 1;
        else
            error('Driver data did not extend far enough in time')
        end
    else
        tInc = (t - cytoData(idx-1,1))/diff(cytoData(idx-1:idx,1));
    end
    c_i = cytoData(idx-1,2) + tInc*diff(cytoData(idx-1:idx,2));
    % c_i = interp1(cytoData(:,1),cytoData(:,2),t);
    ATP_i = cytoData(idx-1,3) + tInc*diff(cytoData(idx-1:idx,3));
    ADP_i = cytoData(idx-1,4) + tInc*diff(cytoData(idx-1:idx,4));
    Pi_i = cytoData(idx-1,5) + tInc*diff(cytoData(idx-1:idx,5));
    ANPtot = 15000*p(38); %uM
    NADtot = 250*p(39); %uM
    c_M = y(1); % mito calcium
    ATP_M = y(2);
    ADP_M = y(7); %ANPtot - ATP_M;
    Pi_M = y(3);
    NADH = y(4);
    NAD = y(8); %NADtot - NADH; % NAD+ in mito is implicit due to mass cons
    if abs(NAD + NADH - NADtot) > 0.01*NADtot
        % warning("NAD and NADH are not staying consistent")
        NAD = NADtot - NADH;
    elseif abs(ATP_M + ADP_M - ANPtot) > 0.01*ANPtot
        % warning("ATP and ADP are not staying consistent in mito")
        ADP_M = ANPtot - ATP_M;
    end
    if ADP_M < -1e-3
        error("ADP cannot be negative!!!!")
    elseif NAD < -1e-3
        error("NAD cannot be negative!!!!")
    end
    V_M = y(5); % IM potential mito
    PiCa_M = y(6); % PiCa precipitate in mito
    % protons = y(7); 
    % b_M = 0.0003; % mito ca buffering capacity
    Kdb = 500.0*p(1);%10.0;
    % PP = 1e4*p(1);
    b_M = 0.01*p(2);%(1 + Kdb*Pi_M/(Kdb+c_M)^2 + 50)^(-1);
    vol_tot = 1000; % arbitrary for fixed cyl geometry
    V_SA_mito = 3.002e-3;
    volFracMito = 0.05;
    vol_mito = vol_tot * volFracMito;
    F = 96485.3321; % Faraday's constant in pC/pmol
    R = 8314.46261815; % gas constant in fJ/(pmol K)
    T = 273.15 + 22; % temperature in K
    N_A = 6.022e11; % molecules per pmol 
    J_to_I = 602.2*vol_mito*F/N_A; % convert flux in uM/s to total current 
    V_MCU = 0.0006*p(3); % uM/s 
    K_trans = 6*p(4); %uM CORRECTED
    K_act = 0.38*p(5); %uM CORRECTED
    n_au = 2.3; %uniporter cooperativity factor
    L = 50*p(6);
    p1 = 0.1*p(7); %1/mV
    V_NCX_M = 0.35*p(8); % uM/s ** DOUBLE CHECK THIS
    p2 = .016*p(9); %1/mV ** DOUBLE CHECK THIS
    k_mPTP = 0.008*p(10); % 1/s
    p3 = 0.05*p(11); % 1/mV
    k_GLY = 450*p(12); % uM/s
    K_MCa = 0.1*p(13); % uM
    K_mNAD = 1*p(14);
    V_O = 600*p(15); % uM/s
    K_O = 100*p(16); % uM IS THIS 100 OR 1000??
    p4 = 177*p(17); % mV
    p5 = 5*p(18); % mV
    V_AGC = 25*p(19); % uM/s 
    K_AGC = 0.14*p(20); % uM
    K_iCa = 0.1*p(21); % uM OR MAYBE 0.3
    V_F1FO = 35000*p(22); % uM/s
    p6 = 190*p(23); %mV 
    p7 = 8.5*p(24); %mV
    K_iATP = 10000*p(25); %uM
    V_ANT = 5000*p(26); %uM/s
    alpha_c = 0.111*p(27);
    alpha_m = 0.139*p(28);
    R_c = alpha_c*ATP_i/ADP_i;
    R_m = ADP_M/(alpha_m*ATP_M);
    f_ANT = 0.5*p(29);
    p8 = 2*p(30); %uM/s per mV
    p9 = -30*p(31); %uM/s
    % V_Pi = 500.0*p(32); % UNKNOWN!!
    V_Pi = 5000*p(32); % estimated from measurements in literature
    % C_p = 1e-5; % pC/mV per sq micron
    % SA_mito = vol_mito/V_SA_mito;
    C_pEff = 1.8*p(33); % uM/mV effective capacitance CHECK!!!
    a1 = 20*p(34);
    a2 = 3.43*p(35);
    % g0_anion = 1000*p(36);
    % V_anion = 230*p(37);
    % KPi_i = 1*p(36);
    KPi_i = 1000*p(36);
    % KPi_M = 1000*p(39);
    g_anion = 0;%g0_anion*Pi_M/1000;%y(7);
    % ginf_anion = Pi_M/1000;%(Pi_M + KPi_M);
    % tau_anion = 100.0;%*p(42);
    koffB = 1.0*p(37);
    konB = koffB / Kdb;
    KPi_M = 10000 * p(40);
    % Ap = 1e-6 * p(39);
    % Bp = 1e-7 * p(40);

    MCUnum = (c_i/K_trans)*(1+c_i/K_trans)^3;
    MCUdenom = (1+c_i/K_trans)^4 + L/(1+c_i/K_act)^n_au;
    J_MCU = V_MCU*(MCUnum/MCUdenom)*exp(p1*V_M);
    J_NCX_M = V_NCX_M*(c_M/c_i)*exp(p2*V_M);
    J_mPTP = k_mPTP*(c_i-c_M)*exp(p3*V_M);
    J_PDH = k_GLY*(c_M/(K_MCa + c_M))*(1/(K_mNAD + NADH/NAD));
    J_O = V_O*(NADH/(K_O+NADH))/(1+exp((V_M-p4)/p5));
    J_AGC = V_AGC*(c_i/(K_AGC+c_i))*(K_iCa/(K_iCa+c_M))*exp(V_M*0.01); % MAYBE REMOVE ELECTROSENS
    J_F1FO = V_F1FO*(1/(1+exp((p6-V_M)/p7)))*(K_iATP/(K_iATP+ATP_M))*(Pi_M/1000);
    ANT_num = (1-R_c*R_m*exp(-F*V_M/(R*T)));
    ANT_denom = (1+(alpha_c*ATP_i/ADP_i)*exp(-f_ANT*F*V_M/(R*T)))*(1+ADP_M/(ATP_M*alpha_m));
    J_ANT = V_ANT*ANT_num/ANT_denom;
    J_Hleak = p8*V_M + p9;
    a3 = 0;%.5; % negative charge transport due to PiP - electrically neutral!
    % J_Pi = V_Pi*Pi_i/(Pi_i + KPi_i);%*(Pi_i - Pi_M);%*exp(V_M*0.01);%(V_M-190)/100;
    J_Pi = V_Pi*(Pi_i/KPi_i - 0.01*Pi_M/KPi_i)/(1+ Pi_i/KPi_i + 0.01*Pi_M/KPi_i);

    % define other charge transport variables
    J_anion = -0;%*g_anion*(V_M - V_anion);%(Pi_M/(Pi_M+KPi_M));%H_Pi; % effective variable for anion transport
    J_MCU_charge = J_MCU;
    J_F1FO_charge = J_F1FO;
    
    % calcium phosphate buffering in mito
    % PC_tot = Pi_M * c_M;
    % if PC_tot >= PP
    %     dPrecip = Ap* (PC_tot)*(PC_tot - PP);
    % else
    %     dPrecip = -Bp*PiCa_M*(PP - PC_tot) ;
    % end
    dPrecip = Pi_M * c_M * konB - koffB * PiCa_M;
    % labelConc = 1.0;
    % KdLabel = 10;
    % koffLabel = 10;

    dydt = zeros(size(y));
    dydt(1) = b_M*(J_MCU - J_NCX_M + J_mPTP - dPrecip);%konPrecip*Pi_M*c_M + konPrecip*Kdb*PiCa_M); % mito ca
    dydt(2) = J_F1FO - J_ANT; % mito ATP
    dydt(3) = -J_F1FO + J_Pi - dPrecip;%konPrecip*Pi_M*c_M + konPrecip*Kdb*PiCa_M; % mito Pi
    dydt(4) = J_PDH - J_O + J_AGC; % NADH - AGC is indirect?
    % dydt(5) = J_to_I*(a1*J_O - a2*J_F1FO - J_Hleak - J_NCX_M - 2*J_MCU - 2*J_mPTP - J_ANT - J_AGC)/(C_p*SA_mito);
    dydt(5) = (a1*J_O - a2*J_F1FO_charge + a3*J_Pi - J_Hleak - J_NCX_M - 2*J_MCU_charge - 2*J_mPTP - J_ANT - J_AGC + J_anion)/(C_pEff);
    dydt(6) = dPrecip;%konPrecip*Pi_M*c_M - konPrecip*Kdb*PiCa_M;
    dydt(7) = -dydt(2); % for ADP
    dydt(8) = -dydt(4); % for NAD
    % dydt(7) = (koffLabel/KdLabel)*c_M*(labelConc-y(7)) - koffLabel*y(7);
    % dydt(7) = (ginf_anion - y(7))/tau_anion;
    fluxes = [J_MCU,J_NCX_M,J_mPTP,J_PDH,J_O,J_AGC,J_F1FO,J_ANT,J_Hleak,J_Pi];
end

function [totError,structOut] = logLikelihoodMito(mitoParams)
    file_names = {'Pi0.001.csv','Pi0.01.csv','Pi0.1.csv','Pi1.csv'};
    folder_names = {'csv_ca_out', 'csv_mito_ca', 'csv_bound_free', 'csv_voltage'};
    dataCells = cell(4,4);
    for i = 1:length(folder_names)
        for j = 1:length(file_names)
            curFile = sprintf('Data/%s/%s', folder_names{i}, file_names{j});
            curDataCell = readcell(curFile);
            dataCells{i,j} = cell2mat(curDataCell(2:end,:));
        end
    end
    CaOutData = dataCells(1,:);
    tTest = 0:.01:400.0;
    phosVals = [0,1,10,100,1000];
    yinit = [0.1; 10000; 0; 50; 150; 0.0; 0.0; 0.0]; % Pi might be closer to 5000 in mito
    yinit(7) = mitoParams(38)*15000 - yinit(2);
    yinit(8) = mitoParams(39)*250 - yinit(4);
    if any(yinit < 0)
        error("Cannot have negative initial conditions")
    end
    options = odeset('RelTol',1e-3,'NonNegative',[1:4,6:8],'MaxStep',0.01);
    mse_Ca = 0;
    mse_buffer = 0;
    mse_voltage = 0;
    fullParams = ones(40,1);
    % varyIdx = [1,32,36,37,38,39,40]; % PP, V_Pi, g0_anion, V_anion, KPi_i, Ap, Bp
    fullParams(:) = mitoParams;%varyIdx) = mitoParams;
    structOut = struct();
    for j = 1:length(phosVals)
        PiCur = phosVals(j);
        if PiCur == 0
            c_i = zeros(size(tTest))' + CaOutData{j}(1,2);
            tVals = tTest';
        else
            cData = CaOutData{j-1};
            c_i = cData(:,2);
            tVals = cData(:,1);
        end
        
        ATPvals = 0.01*ones(size(c_i)); %5000
        ADPvals = 0.01*ones(size(c_i)); %100
        Pivals = PiCur*ones(size(c_i));
        Pivals(tVals < 240) = 0;%
        cytoData = [tVals, c_i, ATPvals, ADPvals, Pivals];
        try [t,y] = ode15s(@mito_dyn, tTest, yinit, options, cytoData, fullParams, tic);
        catch
            totError = 1e6;
            return
        end
        if t(end) < tTest(end)
            totError = 1e6;
            return
        end
        if PiCur == 0 & j==1
            yinit = y(end,:);
            continue
        end
        logicCa = dataCells{2,j-1}(:,1) >= tTest(1) & dataCells{2,j-1}(:,1) <= tTest(end);
        Cainterp = interp1(t,y(:,1),dataCells{2,j-1}(logicCa,1));
        relErrorCa = (Cainterp - dataCells{2,j-1}(logicCa,2))./5;%dataCells{2,j-1}(logicCa,2);
        mse_Ca = mse_Ca + mean(relErrorCa.^2);
        logicRatio = dataCells{3,j-1}(:,1) >= tTest(1) & dataCells{3,j-1}(:,1) <= tTest(end);
        bM = 0.01*fullParams(2);
        Cainterp2 = interp1(t,y(:,1),dataCells{3,j-1}(logicRatio,1));
        CaPhosInterp = interp1(t,y(:,6),dataCells{3,j-1}(logicRatio,1));
        CaRatiointerp = (Cainterp2*(1-bM) + CaPhosInterp*bM)./(Cainterp2*bM);
        relRatio = dataCells{3,j-1}(logicRatio,2);% - dataCells{3,j-1}(1,2);
        relErrorRatio = (CaRatiointerp - relRatio)./100;%relRatio;
        mse_buffer = mse_buffer + mean(relErrorRatio.^2);
        logicVoltage = dataCells{4,j-1}(:,1) >= tTest(1) & dataCells{4,j-1}(:,1) <= tTest(end);
        Voltageinterp = interp1(t,y(:,5),dataCells{4,j-1}(logicVoltage,1));% - y(1,5);
        relVoltage = dataCells{4,j-1}(logicVoltage,2);% - dataCells{4,j-1}(1,2);
        relErrorVoltage = (Voltageinterp - relVoltage)./10;%dataCells{4,j-1}(logicVoltage,2);
        mse_voltage = mse_voltage + mean(relErrorVoltage.^2);
        structOut(j).Cainterp = [dataCells{2,j-1}(logicCa,1), Cainterp];
        structOut(j).CaRatiointerp = [dataCells{3,j-1}(logicRatio,1), CaRatiointerp];
        structOut(j).Voltageinterp = [dataCells{4,j-1}(logicVoltage,1), Voltageinterp];
    end
    totError = mse_Ca + mse_buffer + mse_voltage;
    if isnan(totError)
        totError = 1e6;
    end
    try
        load bestMitoParam.mat objVal
        if totError < objVal
            objVal = totError;
            save('bestMitoParam.mat','objVal', 'mitoParams')
        end
    catch
        fprintf("Problem saving file\r\n")
    end
    % if totError < 10
    %     fprintf("This is a good one\r\n")
    % end
end