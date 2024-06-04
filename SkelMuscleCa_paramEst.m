%% Function for estimating parameters using Particle swarm estimation 
% Input:
%       yinit - Default state variable values
%       expPeaks - Experimental values for ion concentrations
%       lb - Lower bound for estimation
%       ub - upper bound for estimation
%       p_est - Default input parameters
%       Createplot - Logic command to output plots. 1 - Creates plot. 0 -
%       Omits plots
%
% Output:
%       pSol - Particle swarm solution
%       fval - Objective function value
%       exitflag - Reason for ending estimation (Values between -5 and 1.
%                  Explanation of corresponding exit condition avialabel at
%                  https://www.mathworks.com/help/gads/particleswarm.html)

function [pSol,fval,exitflag] = SkelMuscleCa_paramEst(~,lb,ub,yinit,p,Createplot)

psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations', 50, 'SwarmSize', 24);

% pSol Results
% load PSO_25-Apr-2024.mat pSol

% To check default parameter behavior
numParam = length(lb);
pVec = ones(1,numParam) ; 
pToObj(pVec);

delete(gcp('nocreate'))
% parpool(50)
[pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pToObj(pSol)
filename = "PSO_" + date +".mat";
save(filename);

    function objVal = pToObj(pVec)
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
        param = p(:) .* pVec(:);
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


        for n_index = 1 :length(expt_n)
            n = expt_n(n_index);

            % Compute SS with ode15s
            [TimeSS,ySS] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, param,StartTimer,n);
            if size(ySS,1) < length(tSS)
                objVal = penaltyVal;
                return
            end
            if any(isnan(ySS))
                objVal = penaltyVal;
                return
            elseif ySS(end,2) < 800 || ySS(end,2) > 2000
                objVal = penaltyVal;
                return
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
                return
            end

            % Baseline model prediction
            % SS computation
            [TimeSS_base,ySS_base] = SkelMuscleCa_dydtEst(tSS,0, 0, yinit, p(:),StartTimer,n);
            [~,Y_base] = SkelMuscleCa_dydtEst(t,freq(n), 0, ySS_base(end,:), p(:), StartTimer,n);

            if n <= 5
                InterpComp{n} = y(:,8) ; % Calcium Calculations
                CompV{n} = y(:,5) ;
                InterpComp_base{n} = Y_base(:,8);
            elseif n > 5
                InterpComp{n} = y(:,5); % Voltage Calculations
                CompC{n} = y(:,8);
                InterpComp_base{n} = Y_base(:,5);
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

            % figure
            % subplot(1,2,1)
            % plot(0:0.0001:T_max(1),CompV{1},'b','LineWidth',2)
            % subplot(1,2,2)
            % plot(0:0.0001:T_max(5),CompV{5},'b','LineWidth',2)
            % %plot (Time,InterpComp{1},'r')
            % hold off
            % xlabel('Time (s)', 'FontSize',15);
            % title('V_{SL} for expt 1 vs 5', 'FontSize',16);



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
                hold on
                plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2,'Linestyle','--')
                fill(Time_Comp,Ca_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)
                xlabel('Time (s)', 'FontSize',18);                
                title('[Ca^{2+}] for '+ expt_title(i), 'FontSize',18);
                ylabel('Concentration (uM)', 'FontSize',18);

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
                hold on
                plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineStyle', '--' ,'LineWidth',2)
                fill(Time_Comp,V_Comp, 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2)       %'color',[0.00,0.00,1.00]
                xlabel('Time (s)', 'FontSize',18);
                %legend('Computational','Experimental', '95% Confidence Interval', 'FontSize',18)
                ylabel('Membrane Potential (mV)', 'FontSize',18);
                title('V_{SL} for expt - '+ expt_title(i), 'FontSize',18);

            end
            legend('Calibrated','Pre-Calibarated','Experiment', '95% Confidence Interval', 'FontSize',18)
        end
    end
end

