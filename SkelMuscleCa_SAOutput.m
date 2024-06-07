% Output Function for Senstivity Analysis
function QOI = SkelMuscleCa_SAOutput(param)

%param is a matrix of available parameters
% Read in parameters
input = importdata('InputParam.xlsx');
parameters = input.data; % Parameter values
load PSO_25-Apr-2024.mat yinit 
% parameters = p;

QOI = zeros(size(param,1),9);
parpool(50);

   parfor (i = 1 : length(param(:,1))) %, np) %i is the number of samples
    for i = 1 : length(param(:,1))
        p = param(i,:)'.*parameters;
        t = [0 1000];
        tSpan = [0 1];
        freqVec = 100;
        StartTimer = tic;
    
         try
            [tSS,ySS,currtimeSS] = SkelMuscleCa_dydt(t,0, false, yinit, p,StartTimer,2); % SS calculation
           
            yInf = ySS(end,:);   
            %Start ODE soln with SS as initial condition
    
            %% Dynamics - Peak Ca2+, Peak Voltage, Ca2+ Area under curve
            [t_F,y_F,~] = SkelMuscleCa_dydt(tSpan, freqVec, false, yInf, p,StartTimer,2);      % Control.
            for j = 1 : size(y_F,1)
                for k = 1:size(y_F,2)
                    if ~isreal(y_F(j,k))
                        y_F(j,k) = 0;
                    end
                end
            end
            MaxCaF = max(y_F(:,8));                                            % Maximum [Ca2+] conc for control
            MaxVF = max(y_F(:,5));                                             % Maximum Voltage for control
            AreaF = trapz(t_F,y_F(:,8));                                       % Area under curve for control
            QOI(i,:) = [yInf(2), yInf(5), yInf(6), yInf(7),yInf(8), yInf(13),MaxCaF,MaxVF,AreaF]; % Quantitites of interest         
            fprintf('Session %d of %d , time = %.2f s \n',i,size(param,1),currtimeSS);
    
        catch
            yInf = yinit;
            [t_F,y_F,~] = SkelMuscleCa_dydt(tSpan, freqVec, false, yInf, p,tic,2);         % Control.
            for j = 1 : size(y_F,1)
                for k = 1:size(y_F,2)
                    if ~isreal(y_F(j,k))
                        y_F(j,k) = 0;
                    end
                end
            end
            MaxCaF = max(y_F(:,8));                                            % Maximum [Ca2+] conc for control
            MaxVF = max(y_F(:,5));                                             % Maximum Voltage for control
            AreaF = trapz(t_F,y_F(:,8));                                       % Area under curve for control
            QOI(i,:) = [0,0,0,0,0,0,MaxCaF,MaxVF,AreaF];                       % Quantitites of interest
            fprintf('Session %d of %d, with yinit \n',i,size(param,1));

        end
    end
end
