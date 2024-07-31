%% Sweep parameters for quick testing of objective function
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
    3000        % yinit(26) is the initial condition for 'p_i_SR'
    0           % yinit(27) is the initial condition for 'PiCa'
    3000        % yinit(28) is the initial condition for 'Pi_Myo'
    ];

% Importing parameters 
param = importdata('InputParam1.xlsx');

%% sweep across ranges for sensitive parameters
% you can test different values of default parameters here as well!
p0 =  param.data;
numParam = length(p0);
highSensIdx = [1,3,5,6,8,9,10,11,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,43,44,45,46,51,52,53,69,70]; % a vector listing the indices of all parameters we are still including (higher sensitivity values)
% p0([20,42,43]) = 0.2*p0([20,42,43]); % NCX, SERCA, PMCA
% p0(44) = 0.1*p0(44); % leak SR1
pVec = ones(1,numParam);
samples = 100; % number of random samples to generate
sigmaTest = 0.2; % geometric standard deviation controlled extent of random changes in parameters
randPop = exp(sigmaTest*randn([samples, length(highSensIdx)]));
objVals = zeros(samples, 1);
simSaved = cell(samples, 1);
figure
hold on
for i = 1:samples
    [objVals(i), simSavedCur] = pToObj(randPop(i,:), p0, yinit, false);
    if isempty(simSavedCur{2})
        continue
    end
    % look at experiment 2 results used in fitting
    simSaved{i} = simSavedCur{2};
    % note that the index of 9 actually indicates variable number 8
    % (calcium), the first column of simSaved{i} is time
    plot(simSaved{i}(:,1), simSaved{i}(:,9))
    drawnow
    fprintf("%d\n", i)
end

% plot the best solution over an extended time
[~,bestIdx] = min(objVals);
pBest = p0(:);
pBest = randPop(bestIdx,:).*pBest(highSensIdx);
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pBest, tic, 2);
tSol = 0:.0001:10;
[Time,Y] = SkelMuscleCa_dydt(tSol, 100, 0, ySS(end,:), pBest, tic, 2);
figure
plot(Time, Y(:,8))

%% plot the solution from PSO over time compared to default parameters
p0 =  param.data;
load PSO_30-Jul-2024.mat pSol
% pSol(12) = pSol(12)*0.2;
highSensIdx = [1,3,5,6,8,9,10,11,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,43,44,45,46,51,52,53,69,70]; % a vector listing the indices of all parameters we are still including (higher sensitivity values)
pPSO = p0(:);
pPSO(highSensIdx) = pSol(:).*pPSO(highSensIdx);
fprintf("Objective value from PSO is %.3f\n", pToObj(pSol, p0, yinit, true))

%% plot the solution from PSO over an extended time
[TimeSS,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, p0, tic, 2);
tSol = 0:.0001:10;
[Time,Y] = SkelMuscleCa_dydt(tSol, 20, 0, ySS(end,:), p0, tic, 1);

figure
plot(Time, Y(:,8))
figure
plot(Time, Y(:,5))

%% Function for calculating the objective value for estimation

function [objVal, simSaved] = pToObj(pVec, p, yinit, Createplot)
%% Function for calculating the objective value for estimation
% Input:
%       pVec - Vector of particles
% Output:
%       objVal - Minimum error between the experimental and model
%       output

%Initialize values
T_max = zeros(1,9);
InterpExpt = cell(1,9);
Expt_t = cell(1,9);
InterpComp = cell(1,9);
InterpComp_base = cell(1,9);
CompV = cell(1,5);
CompC = cell(1,5);
StartTimer = tic;
param = p(:); %.* pVec(:); % initialize all parameter values to defaults
highSensIdx = [1,3,5,6,8,9,10,11,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,40,43,44,45,46,51,52,53,69,70]; % a vector listing the indices of all parameters we are still including (higher sensitivity values)
param(highSensIdx) = param(highSensIdx) .* pVec(:);

tSS = 0:1000;
penaltyVal = 100000;

%Experimental Data
%Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C], [K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_T MJ_V]};
load Exptdata.mat Expt
freq = [100, 100, 67, 67,67,60, 60, 60,60];
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

simSaved = cell(8,1);

for n_index = 1 :length(expt_n)
    n = expt_n(n_index);

    % Compute SS with ode15s
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
    % elseif ySS(end,2) < 800 || ySS(end,2) > 2000 || param(12) > 0.5*ySS(end,2) %changed from 0.4
    %     objVal = penaltyVal;
    %     simSaved{n} = zeros(1,29);
    %     return
    end
    yinf = ySS(end,:);

    % SS plots

    % if Createplot
    %     figure
    %     SS_index = [6,7,13,5,8];
    %     for a = 1:length(SS_index)
    %         var = SS_index(a);
    %         if a < 4
    %             plot(TimeSS,ySS(:,var),'LineWidth',2)
    %             hold on
    %             ylabel('Concentration (\mu M)');
    %             legend('[Na^+]', '[Cl^-]', '[K^+]');
    %         elseif a > 3
    %             plot(TimeSS,ySS(:,var),'LineWidth',2)
    %             hold on
    %             legend('V_{SL} (mV)', '[Ca^{2+}] (\muM');
    %         end
    %         title('Variables during steady state', 'FontSize',16);
    %         xlabel('Time (s)', 'FontSize',15);
    %     end
    %     hold off
    % end


    %% Calculate Dynamics
    t = 0:0.0001:T_max(n);
    [~,y] = SkelMuscleCa_dydtEst(t,freq(n), 0, yinf, param,StartTimer,n);
    if size(y,1) < length(t)
        objVal = penaltyVal;
        simSaved{n} = zeros(1,29);
        return
    end

    simSaved{n} = [t(:),y];

    % Baseline model prediction
    % SS computation
    if Createplot
        [TimeSS_base,ySS_base] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, p(:),StartTimer,n);
        [~,Y_base] = SkelMuscleCa_dydtEst(t,freq(n), 0, ySS_base(end,:), p(:), StartTimer,n);
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
        plot(0:0.0001:T_max(i),InterpComp{i},'LineWidth',2, 'color',[0.00,0.00,1.00]) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',2, 'Color',[0.10,0.85,0.83])
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2,'Linestyle','--')
        fill(Time_Comp,Ca_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)
        xlabel('Time (s)', 'FontSize',18);
        title('[Ca^{2+}] for '+ expt_title(i), 'FontSize',18);
        ylabel('Concentration (uM)', 'FontSize',18);
        prettyGraph

    end

    for index = 3:4
        i = expt_n(index); %% Update according to the number of Voltage expts used
        subplot(2,2,index)
        x = 0:0.0001:T_max(i);
        Time_Comp = [x, fliplr(x)];
        V_Comp = [InterpExpt{i} + 2*sigma_V, fliplr(InterpExpt{i} - 2*sigma_V)];
        plot(0:0.0001:T_max(i),InterpComp{i},'LineWidth',2, 'color',[0.49,0.18,0.56] ) %'b',
        hold on
        plot(0:0.0001:T_max(i),InterpComp_base{i},'LineWidth',2, 'Color',[0.10,0.85,0.83])
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineStyle', '--' ,'LineWidth',2)
        fill(Time_Comp,V_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)       %'color',[0.00,0.00,1.00]
        xlabel('Time (s)', 'FontSize',18);
        ylabel('Membrane Potential (mV)', 'FontSize',18);
        title('V_{SL} for expt - '+ expt_title(i), 'FontSize',18);
        prettyGraph

    end
    legend('Calibrated','Pre-Calibrated','Experiment', '95% Confidence Interval', 'FontSize',12)
end
end