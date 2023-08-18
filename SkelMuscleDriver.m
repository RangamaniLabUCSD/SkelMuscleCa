%% start parameter optimization
lb = 0.5*ones(3,1); 
ub = 2*ones(3,1);
tSpan = [0 1];
freqVals = 0;

% Add arrays for Ion concentration values 

Ca_i_exp = [0.1,0.05,0.1];          %uM
Ca_SR_exp = [1500,1300,1300];       %uM
Na_i_exp = [14700,12700, 12700];    %uM
K_i_exp = [154500, 180300, 150900]; %uM
Voltage_PM_exp = [-88, -60, -79];   %mV

expPeaks = [Ca_i_exp; Ca_SR_exp; Na_i_exp; K_i_exp; Voltage_PM_exp]; %Value
%pSol = SkelMuscleCa_paramEst(tSpan, 0, expPeaks, lb, ub);

param = importdata('InPutParam.xlsx');
p = param.data; % Parameter values

selParam = ones(length(p),1);

yinit = SkelMuscleCa_SS(tSpan,false,selParam); %Calculating SS conc with psol estimates.
%[~,~,SSParam] = SkelMuscleCa(tSpan, 0, false, yinit, p); %Solving ODE for state variables.

%% Dynamics 
tSpan = [0 5];
freqVec = logspace(-1,2,100);    

% [Ca2+] conc over time for different ATP level conditions.

figure
hold on
yStored = cell(size(freqVec)); % LowATP condition = ~50% of ATP
for i = 1:length(freqVec)
    [t_T,y_T] = SkelMuscleCa(tSpan, freqVec(i), 1, yinit, selParam); %Solving ODE for state variables.
    plot(t_T,y_T(:,9))
    title(sprintf('%s',"Intracellular Ca^{2+} concentration"+ max(tSpan)+" s, Low ATP")) 
    xlabel('Time')
    ylabel('Intracellular Ca^{2+} concentration [\muM]')
    xlim(tSpan)
    ylim([0 10])
    yStored{i} = [t_T, y_T];
    
end 
saveas(gcf,'Cavstime_True.png','png')

figure
hold on 
yStoredF = cell(size(freqVec)); % Control ATP condition
for i = 1:length(freqVec)
    [t_F,y_F] = SkelMuscleCa(tSpan, freqVec(i), 0, yinit, selParam);
    plot(t_F,y_F(:,9))
    title(sprintf('%s',"Intracellular Ca^{2+} concentration"+ max(tSpan)+" s, Control")) 
    xlabel('Time')
    ylabel('Intracellular Ca^{2+} concentration [\muM]')
    xlim(tSpan)
    ylim([0 10])
    yStoredF{i} = [t_F,y_F];
end 
saveas(gcf,'Cavstime_False.png','png')


%% Maximum [Ca2+] conc & Semilog plot
figure 
MaxCa = zeros(size(freqVec)); 
MaxCaF = zeros(size(freqVec));
for i = 1:length(freqVec) 
    MaxCa(i) = max(yStored{i}(:,10)); 
    MaxCaF(i) = max(yStoredF{i}(:,10));
end

semilogx(freqVec,MaxCa,'-s');
hold on 
semilogx(freqVec,MaxCaF,'-s');
title(sprintf('%s', "Peak cytosolic Ca^{2+} concentration for " + max(tSpan)+" s")) % + ATP_level))
ylabel('Peak cytosolic Ca^{2+} concentration [\muM]') 
xlabel('Stimulus frequency [Hz]')
legend('Low ATP','Control')
saveas(gcf,'SemilogMaxCa.png','png')
 
%% Area under curve for diff Low ATP level inputs
Area = zeros(size(freqVec));
AreaF = zeros(size(freqVec));
for i = 1:length(freqVec)
    Area(i) = trapz(yStored{i}(:,1),yStored{i}(:,10));
   AreaF(i) = trapz(yStoredF{i}(:,1),yStoredF{i}(:,10));
end
 
figure
semilogx(freqVec,AreaF,'b')
hold on
semilogx(freqVec,Area) 
title(sprintf('%s',"Area vs Frequency for "+ max(tSpan)+" s")) % + ATP_level))
xlabel('Stimulus frequency [Hz]')
ylabel('Area [\muM s]')
legend('Control','Low ATP')
saveas(gcf,'Area.png','png') %Saves the area under the curve plot. 

%% Time when max calcium conc is achieved 
freq_idx = find(MaxCa == max(MaxCa)); 
Ca_idx = find(yStored{freq_idx}(:,10) == max(MaxCa));                                        
time = yStored{freq_idx}(Ca_idx,1);
