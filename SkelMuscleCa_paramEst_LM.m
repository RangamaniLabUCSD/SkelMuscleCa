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

psOptions = optimoptions('particleswarm','SwarmSize',500,'UseParallel',true,'HybridFcn',@fmincon,...
    'PlotFcn','pswplotbestf','Display','iter','MaxStallIterations',100); %%Increased from 15

numParam = length(lb);
pVec = ones(1,numParam); %41
pToObj(pVec)
delete(gcp('nocreate'))
%parpool(12)
[pSol,fval,exitflag] = particleswarm(@pToObj,numParam,lb,ub,psOptions);
pToObj(pSol)


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
        param = p_est .* pVec';
        tSS = [0 1000];

        %Interpolating experimental values
        for m = 3 %:9
            T_max(m) = max(Expt{m}(:,1))/1000;
            Expt_t{m} = Expt{m}(:,1)/1000;
            InterpExpt{m} = interp1(Expt_t{m},Expt{m}(:,2),0:0.0001:T_max(m));
        end

        penaltyVal = 100000;
        count = 0;
        for n = 3 %:9         
            %% Colmpute SS yinit with LM algorithm.
            initialGuess = yinit;
            options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIter', 10000); %,'UseParallel',true);
            lb1 = 0.5 * initialGuess;
            ub1 = 2 * initialGuess;
            lb1(5) = -140; %2 * initialGuess(5);
            ub1(5) = 0.5 * initialGuess(5);

            [y_LM] = lsqnonlin(@(initialGuess) objectivefn(initialGuess,n,param),initialGuess, lb1, ub1, options);
            if any(isnan(y_LM))
                objVal = penaltyVal;
                return
            end
           
            %% Compute SS with ode15s            
            % [~,ySS] = SkelMuscleCa1_SS(tSS,0, 0, y_LM, param,StartTimer,n);
            % 
            % if any(isnan(ySS))
            %     objVal = penaltyVal;
            %     count = count+1;    
            %     return
            % end
            yinf = y_LM; %ySS(end,:);
             
            %     yinf_ratio(n,:) = yinf ./ y_LM';
            % else
            %     yinf_ratio(n,:) = yinf ./ y_LM;
            % end
            
            %fprintf('----- SS Complete ------ ');           
          
            %% Calculate Dynamics
            t = [0 T_max(n)];
            [Time,y] = SkelMuscleCa1_SS(t,freq(n), 0, yinf, param,StartTimer,n);
            %[Time,y] = SkelMuscleCa1_SS([0 t(n)],freq(n), 0, yinit, p,tic,n);
            if size(y,1)<2
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
            yinit = yinf;
        end

        %fprintf('----- Dynamics Complete ------ ');
        %fprintf(' ----- %d ----- \n',count); 
        %toc(StartTimer)

        %%
        % Plots       
        % figure
        % idx = [6,2,7,8,9,14,15];
        % idx_Title = ["Voltage_{PM}", "Ca_{SR}" ,"[Na^+]_i" ,"[Cl^-]_i" ,"[Ca^{2+}]", "[K^+]_i", "[CaParv]"];
        % for j = 1:length(idx)
        % 
        %     subplot(4,2,j)
        %     index = idx(j);           
        %     for i = 1:9
        %         if j == 1
        %             plot(t_SS{i},ySS{i}(:,index),'LineWidth',2);
        %         else
        %         semilogy(t_SS{i},ySS{i}(:,index),'LineWidth',2);
        %         end
        %         hold on
        %     end
        %     %hold on
        %         title(idx_Title(j));
        %         xlabel('Time (s)');                
        %         if j == 1
        %             ylabel('V_{PM} mV')
        %         else
        %             ylabel('Conc (uM)')
        %         end            
        % end

        % figure
        % %for i = 1:5
        % % subplot(3,2,i)
        % i = 3;
        % plot(0:0.0001:T_max(i),InterpComp{i},'b','LineWidth',2)
        % hold on
        % plot(0:0.0001:T_max(i),InterpExpt{i},'r','LineWidth',2)
        % xlabel('Time (s)');
        % legend('Computational','Experimental')
        % title('Computational vs expt \Delta[Ca^{2+}] concentration for expt - Rincon',i); % num2str(i));
        % ylabel('\Delta[Ca^{2+}] Concentration (uM)');
        % fontsize(15,"points")
        %end
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
        % delta = cell(1,9);
        sum_delta = zeros(1,9);

        % for j = 3 %:9
        %     weight = length(InterpExpt{j});
        %     sigma_C = 1; %(0.1 + (0.05 * InterpExpt{j})).^2;
        %     if j < 6
        %         delta{j} = ((InterpComp{j} - InterpExpt{j}).^2) ./ (weight * sigma_C) ;
        %     else
        %         if j > 5
        %             delta{j} = ((InterpComp{j} - InterpExpt{j}).^2)/(weight* 25) ;
        %         end
        %         sum_delta(j) = sum(delta{j});
        %     end
        % end 
        delta = ((InterpComp{3} - InterpExpt{3}).^2) ./ 1.1; %./length(InterpExpt{1});% sigma_C = 1
        objVal = sum(delta);
    end
end

