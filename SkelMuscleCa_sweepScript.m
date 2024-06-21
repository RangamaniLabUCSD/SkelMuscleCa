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
    ];

% Parameter values
param = importdata('InputParam1.xlsx');

%%
p0 =  param.data;
numParam = length(p0);
p0([20,42,43]) = 0.2*p0([20,42,43]); % NCX, SERCA, PMCA
p0(44) = 0.1*p0(44); % leak SR1
pVec = ones(1,numParam);
samples = 100;
randPop = exp(0.2*randn([samples, length(p0)]));
objVals = zeros(samples, 1);
simSaved = cell(samples, 1);
figure
hold on
for i = 1:samples
    [objVals(i), simSaved, fluxesSaved] = pToObj(randPop(i,:), p0, yinit, false);
    simSaved{i} = simSaved{7};
    fluxesSaved{i} = fluxesSaved{7};
    plot(simSaved{i}(:,1), simSaved{i}(:,2))
    drawnow
    fprintf("%d\n", i)
end

% plot the best solution over an extended time
[~,bestIdx] = min(objVals);
pBest = randPop(bestIdx,:).*p0';
[~,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pBest, tic, 2);
tSol = 0:.0001:10;
[Time,Y] = SkelMuscleCa_dydt(tSol, 100, 0, ySS(end,:), pBest, tic, 2);
figure
plot(Time, Y(:,8))

%% plot the solution from PSO over time
param = importdata('InputParam1.xlsx');
p0 =  param.data;
load PSO_25-Apr-2024.mat pSol
% pSol(12) = pSol(12)*0.2;
pPSO = pSol.*p0';
fprintf("Objective value from PSO is %.3f\n", pToObj(pSol, p0, yinit, false))
[TimeSS,ySS] = SkelMuscleCa_dydt([0 1000],0, 0, yinit, pPSO, tic, 2);
tSol = 0:.0001:10;
[Time,Y] = SkelMuscleCa_dydt(tSol, 10, 0, ySS(end,:), pPSO, tic, 1);

figure
plot(Time, Y(:,8))

%% Function for calculating the objective value for estimation

function [objVal, simSaved, fluxesSaved] = pToObj(pVec, p, yinit, Createplot)
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
CompV = cell(1,5);
CompC = cell(1,5);
StartTimer = tic;
param = p(:) .* pVec(:);
tSS = 0:1000;
penaltyVal = 100000;

%Experimental Data
%Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C], [K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_T MJ_V]};
load Exptdata.mat Expt
freq = [100, 100, 67, 67,67, 60, 60, 60, 60];
expt_title = ["Rincon","Rincon", "Baylor & Hollingworth", "Hollingworth", "Baylor & Hollingworth", "Yonemura","Bibollet", "Miranda","Wallinga"];
expt_n = [1 2 7 8]; % Index of the experimental used for estimation

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

fluxesSaved = cell(size(expt_title));
simSaved = cell(size(expt_title));
for n_index = 1 :length(expt_n)
    n = expt_n(n_index);

    % Compute SS with ode15s
    [TimeSS,ySS] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, param,StartTimer,n);
    if size(ySS,1) < length(tSS)
        objVal = penaltyVal;
        return
    end
    yinf = ySS(end,:);

    % SS plots

    if Createplot
        figure
        SS_index = [6,7,13,5,8];
        for a = 1:length(SS_index)
            var = SS_index(a);
            if a < 4
                plot(TimeSS,ySS(:,var),'LineWidth',2)
                hold on
                ylabel('Concentration (\mu M)');
                legend('[Na^+]', '[Cl^-]', '[K^+]');
            elseif a > 3
                plot(TimeSS,ySS(:,var),'LineWidth',2)
                hold on
                legend('V_{SL} (mV)', '[Ca^{2+}] (\muM');
            end
            title('Variables during steady state', 'FontSize',16);
            xlabel('Time (s)', 'FontSize',15);
        end
        hold off
    end


    % Calculate Dynamics
    t = 0:0.0001:T_max(n);
    [~,y,~,fluxes] = SkelMuscleCa_dydtEst(t,freq(n), 0, yinf, param,StartTimer,n);
    fluxesSaved{n} = fluxes;
    if n <= 5
        InterpComp{n} = y(:,8) ; % Calcium Calculations
        CompV{n} = y(:,5) ;
    elseif n > 5
        InterpComp{n} = y(:,5); % Voltage Calculations
        CompC{n} = y(:,8);
    end
    simSaved{n} = [t', InterpComp{n}];
end

% Objective Value Calc
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
        Error{j} = ((delta{j} ./ sigma_C ) .^ 2 )./ weight; %     (((InterpComp{j}' - InterpExpt{j})./sigma_C).^2)./weight ;
    elseif j > 5
        Error{j} = ((delta{j} ./ sigma_V ) .^ 2 )./ weight; % (((InterpComp{j}' - InterpExpt{j}) ./ sigma_V).^2) ./ weight;
    end
    sum_delta(j) = sum(Error{j});

end

objVal = sum(sum_delta);

% Plots
if Createplot

    figure
    for index = 1:2 % Update according to the number of Calcium expts used
        i = expt_n(index);
        subplot(2,2,index)
        plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        hold on
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        xlabel('Time (s)', 'FontSize',15);
        legend('Computational','Experimental', 'FontSize',15)
        title('Computational vs expt \Delta[Ca^{2+}] concentration for expt - '+ expt_title(i), 'FontSize',16);
        ylabel('\Delta[Ca^{2+}] Concentration (uM)', 'FontSize',15);

    end

    for index = 3:4
        i = expt_n(index); %% Update according to the number of Voltage expts used
        subplot(2,2,index)
        plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        hold on
        plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        xlabel('Time (s)', 'FontSize',15);
        legend('Computational','Experimental', 'FontSize',15)
        ylabel('V_{PM} (mV)', 'FontSize',15);
        title('Computational vs expt V_{PM} for expt - '+ expt_title(i), 'FontSize',16);

    end
end
end