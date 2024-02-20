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

psOptions = optimoptions('particleswarm','SwarmSize',100,'UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations',50); 

numParam = length(lb);
pVec = ones(1,numParam); %41
% load LM_PSO_2_16_1.mat pSol_LM
% pVec = pSol_LM;
pToObj(pVec)
delete(gcp('nocreate'))
parpool(50)
[pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pToObj(pSol)
save('LM_PSO_2_16_1.mat')

    function objVal = pToObj(pVec) % Obj function should be scalar

        %Initialize values
        %fprintf(' ------ Code Started ------ ');
        T_max = zeros(1,9);        
        InterpExpt = cell(1,9);
        Expt_t = cell(1,9);
        InterpComp = cell(1,9);
        %yinf_ratio = zeros(9,17);

        load Exptdata.mat Expt
        %Expt = {[R_t R_C],[R_MP_t R_MP_C] [HB_t HB_C], [H_t H_C],[HB_MP_t HB_MP_C], [K_t K_V], [B_t B_V] , [M_t M_V], [W_t W_V], [MJ_T MJ_V]};
        StartTimer = tic;
        freq = [100, 100, 67, 67,67, 60, 60, 60,60];
        param = p_est(:) .* pVec(:);
        tSS = [0 1000];

        %Interpolating experimental values
        for m = 1 %:9
            T_max(m) = max(Expt{m}(:,1))/1000;
            Expt_t{m} = Expt{m}(:,1)/1000;
            InterpExpt{m} = interp1(Expt_t{m},Expt{m}(:,2),0:0.0001:T_max(m));
        end

        penaltyVal = 100000;
        count = 0;
        for n = 1 %:9         
            %% Colmpute SS yinit with LM algorithm.
            initialGuess = yinit;
            options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIter', 10000); %,'UseParallel',true);
            lb1 = 0.5 * initialGuess;
            ub1 = 2 * initialGuess;
            lb1(5) = 2 * initialGuess(5);
            ub1(5) = 0.5 * initialGuess(5);

            [y_LM] = lsqnonlin(@(initialGuess) objectivefn(initialGuess,n,param),initialGuess, lb1, ub1, options);
            if any(isnan(y_LM))
                objVal = penaltyVal;
                return
            end
             
            %% Compute SS with ode15s         
            [~,ySS] = SkelMuscleCa1_SS(tSS,0, 0, yinit, param,StartTimer,n);

            if any(isnan(ySS))
                objVal = penaltyVal;
                count = count+1;    
                return
            end
            yinf = ySS(end,:);            
       
          
            %% Calculate Dynamics
            t = [0 T_max(n)];
            [Time,y] = SkelMuscleCa1_SS(t,freq(n), 0, yinf, param,StartTimer,n);
            %[Time,y] = SkelMuscleCa1_SS([0 t(n)],freq(n), 0, yinit, p,tic,n);
            if size(y,1) < 2
                objVal = penaltyVal;
                count = count+1;
                return
            end
            if n <= 5
                Comp = y(:,8); % Calcium Calculations
            elseif n > 5
                Comp = y(:,5); % Voltage Calculations
            end
            InterpComp{n} = interp1(Time,Comp,0:0.0001:T_max(n));            
            %yinit = yinf;
        end

        %%
        % Plots       
        % figure
        % for i = 1
        %     plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        %     hold on
        %     plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        %     xlabel('Time (s)');
        %     legend('Computational','Experimental')
        %     title('Computational vs expt \Delta[Ca^{2+}] concentration for expt - Rincon',i); % num2str(i));
        %     ylabel('\Delta[Ca^{2+}] Concentration (uM)');
        %     fontsize(15,"points")
        % end

        % figure
        % for i = 6:9        
        %     subplot(2,2,i-5)            
        %     plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        %     hold on
        %     plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        %     xlabel('Time (s)');
        %     legend('Computational','Experimental')
        %     ylabel('V_{PM} (mV)');
        %     title('Computational vs expt V_{PM} for expt %d',i); %', num2str(i));
        %     fontsize(15,"points")
        % end 
        
        %% Objective Value Calc       
        delta = cell(1,9);
        sum_delta = zeros(1,9);

        for j = 1 %:9
             weight = length(InterpExpt{j});
             % sigma_C = 1.1; %(0.1 + (0.05 * InterpExpt{j})).^2;
             if j < 6
                 delta{j} = ((InterpComp{j} - InterpExpt{j}).^2); % ./ (weight * sigma_C) ;
            else
                if j > 5
                    delta{j} = ((InterpComp{j} - InterpExpt{j}).^2)/(weight* 25) ;
                end
             end
                sum_delta(j) = sum(delta{j});
           
         end 
        
        objVal = sum(sum_delta);
    end
end
