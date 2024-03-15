%% Function definition
% Input:
%       tSpan - Time span
%       freqVals - Frequency
%       expPeaks - Experimental values for ion concentrations
%       lb - Lower bound for estimation
%       ub - upper bound for estimation
%       p_est - baseline input parameters
%
% Output:
%       pSol - Particle swarm solution

function [pSol,fval,exitflag] = SkelMuscleCa_paramEst_LM(~,lb,ub,yinit,p_est)

psOptions = optimoptions('particleswarm','UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations',100);%, 'SwarmSize',20);

numParam = length(lb);
pVec = ones(1,numParam); %41
% load PSO_11-Mar-2024_2.mat pSol
% pVec = pSol;
% error_i = pToObj(pVec);
%%
% pVec = 0.8 + (1.25 - 0.8)* rand(100,numParam);
% for a = 1:100
%     error(a) = pToObj(pVec(a,:));
% end
%%
% pVec = [1.13221318715608	0.910784829890607	1.08766990881102	0.960520363703168	0.881220949480694	1.14980749044848	1.01784366722258	0.832650430164826	1.13204434984506	1.21134064952321	0.853989752515383	1.02556941419551	1.04505055423832	1.24069003552003	1.20312237401022	1.05626370998520	1.21757747171939	0.906544295531512	0.811685340553057	0.948058482316822	1.18105151245671	1.00575732420527	0.872659170917175	1.24424664541645	1.20403751216216	1.15221082903596	1.11655793288450	0.878720899125972	1.06486411347693	1.14367568901613	0.890078146487194	1.24850588834879	1.24985485041859	0.811471290427414	1.03705854781936	1.02675718342494	0.917570025970154	1.14125840581366	1.10821470759528	0.890444308368329	1.17040333644423];
% pToObj(pVec)
delete(gcp('nocreate'))
parpool(50)
[pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pToObj(pSol)
filename = "PSO_" + date +"_2.mat";
save(filename);


    function objVal = pToObj(pVec) % Obj function should be scalar

        %Initialize values
        %fprintf(' ------ Code Started ------ ');
        T_max = zeros(1,9);        
        InterpExpt = cell(1,9);
        Expt_t = cell(1,9);
        InterpComp = cell(1,9);
        CompV = cell(1,5);
        CompC = cell(1,5);
        %yinf_ratio = zeros(9,17);

        load Exptdata.mat Expt
        %Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C], [K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_T MJ_V]};
        StartTimer = tic;
        freq = [100, 100, 67, 67,67, 60, 60, 60,60];
        expt_title = ["Rincon","Rincon", "Baylor & Hollingworth", "Hollingworth", "Baylor & Hollingworth", "Yonemura","Bibollet", "Miranda","Wallinga"];
        param = p_est(:) .* pVec(:);
        tSS = 0:1000;

        expt_n = [1 5 7 8]; % 1:9; % [1 8];%

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

        penaltyVal = 100000;
        count = 0;
        for n_index = 1 :length(expt_n) %:9 
            n = expt_n(n_index);
            %% Colmpute SS yinit with LM algorithm.
            %initialGuess = yinit;
            %options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIter', 10000); %,'UseParallel',true);
            %lb1 = 0.5 * initialGuess;
            %ub1 = 2 * initialGuess;
            %lb1(5) = 2 * initialGuess(5);
            %ub1(5) = 0.5 * initialGuess(5);

            %[y_LM] = lsqnonlin(@(initialGuess) objectivefn(initialGuess,n,param),initialGuess, lb1, ub1, options);
            %if any(isnan(y_LM))
            %    objVal = penaltyVal;
            %    return
            %end
             
            %% Compute SS with ode15s         
            [TimeSS,ySS] = SkelMuscleCa1_SS(tSS,0, 0, yinit, param,StartTimer,n);
            if size(ySS,1) < length(tSS)
                objVal = penaltyVal;
                count = count+1;
                return
            end
            if any(isnan(ySS))
                objVal = penaltyVal;
                count = count+1;    
                return
            end
            yinf = ySS(end,:);   
             
            % figure
            % plot(TimeSS,ySS(:,6),'b','LineWidth',2)
            % hold on
            % plot(TimeSS,ySS(:,7),'r','LineWidth',2)
            % hold on
            % plot(TimeSS,ySS(:,13),'g','LineWidth',2)
            % xlabel('Time (s)')
            % legend('[Na^+]', '[Cl^-]', '[K^+]');
            % hold off
            % figure
            % plot(TimeSS,ySS(:,5),'b','LineWidth',2)
            % hold on
            % plot(TimeSS,ySS(:,8),'r','LineWidth',2)
            % hold off
            % xlabel('Time (s)')
            % legend('V_{SL}', '[Ca^{2+}]');

        
       
          
            %% Calculate Dynamics
            t = 0:0.0001:T_max(n);
            [Time,y] = SkelMuscleCa1_SS(t,freq(n), 0, yinf, param,StartTimer,n);
            %[Time,y] = SkelMuscleCa1_SS([0 t(n)],freq(n), 0, yinit, p,tic,n);
            if size(y,1) < length(t) 
                objVal = penaltyVal;
                count = count+1;
                return
            end
            if n <= 5
                InterpComp{n} = y(:,8) ; % Calcium Calculations
                CompV{n} = y(:,5) ;
            elseif n > 5
                InterpComp{n} = y(:,5); % Voltage Calculations
                CompC{n} = y(:,8);
            end
            %InterpComp{n} = interp1(Time,Comp,0:0.0001:T_max(n));            
            %yinit = yinf;
        end

        %% Plots   
        % figure
        % subplot(1,2,1)
        % plot(0:0.0001:T_max(1),CompV{1},'b','LineWidth',2)
        % subplot(1,2,2)
        % plot(0:0.0001:T_max(5),CompV{5},'b','LineWidth',2)
        % %plot (Time,InterpComp{1},'r')
        % hold off
        % xlabel('Time (s)');
        % title('V_{SL} for expt 1 vs 5'); %', num2str(i));

        % figure
        % for index = 1:2 %
        %     i = expt_n(index);
        %     subplot(2,2,index)
        %     plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        %     hold on
        %     plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        %     xlabel('Time (s)');
        %     legend('Computational','Experimental')
        %     title('Computational vs expt \Delta[Ca^{2+}] concentration for expt - '+ expt_title(i)); % num2str(i));
        %     ylabel('\Delta[Ca^{2+}] Concentration (uM)');
        %     fontsize(15,"points")
        % end
        % 
        % 
        % % figure
        % for index = 3:4 
        %     i = expt_n(index); %:9
        %     subplot(2,2,index)
        %     plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        %     hold on
        %     plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        %     xlabel('Time (s)');
        %     legend('Computational','Experimental')
        %     ylabel('V_{PM} (mV)');
        %     title('Computational vs expt V_{PM} for expt - '+ expt_title(i)); %');
        %     fontsize(15,"points")
        % end

        %% Objective Value Calc       
        delta = cell(1,9);
        sum_delta = zeros(1,9);
        Error = cell(1,9);
       
        for j_index = 1 :length(expt_n) %:9
            j = expt_n(j_index);
            weight = length(InterpExpt{j}) ;
            sigma_C = (0.05 * InterpExpt{j}); %
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
    end
end

