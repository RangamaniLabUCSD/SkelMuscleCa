clear 
load('MorrisResults10-24.mat')

% QOI(i,:) = [yInf(2), yInf(5), yInf(6), yInf(7),yInf(8), yInf(13), yInf(23), MaxCaF,MaxVF, MaxPost, AvgF, AvgPost, AvgVolt, VoltWidth];

% graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion","SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force", "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force","Avg Voltage", "AP Width"};
graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion","SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force", "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force","Avg Voltage", "AP Width"};
expt_names = repmat(graph_names,4)';

n = 1;
muVec = zeros(size(MorrisAnalysis.Results.MuStar));

ActionPotential = [1:6, 8:10, 13,14,16:30,33,36:40,45,76:82, 85:88];
CaInflux_SR_Release = [7,11,15,32,35,83,90:92];
CaBuffering = [46:54,58,68:72,74,75,96,97];
CaEfflux_SOCE = [12,31,34,41:44,89,95];
Crossbridge_Cycle = [55:57,59:67,73,84,93,94];


for k= 57 %1:length(graph_names)
    figure
    scatter(MorrisAnalysis.Results.MuStar(ActionPotential,k),MorrisAnalysis.Results.Std(ActionPotential,k),'filled','MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(CaInflux_SR_Release,k),MorrisAnalysis.Results.Std(CaInflux_SR_Release,k),'filled','MarkerFaceColor',[0.4660 0.6740 0.1880])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(CaBuffering,k),MorrisAnalysis.Results.Std(CaBuffering,k),'filled','MarkerFaceColor',[0 0.4470 0.7410])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(CaEfflux_SOCE,k),MorrisAnalysis.Results.Std(CaEfflux_SOCE,k),'filled','MarkerFaceColor',	[1 0 0])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(Crossbridge_Cycle,k),MorrisAnalysis.Results.Std(Crossbridge_Cycle,k),'filled','MarkerFaceColor',[0.4940 0.1840 0.5560])
    xlabel('mu*')
    ylabel('simga')
    yticklabels({})
    % title(graph_names{k})
    hold on
    ten_percent = max(MorrisAnalysis.Results.MuStar(:,k))*0.1;
    xline(ten_percent , '--b' )
    legend('Action Potential', 'Calcium Influx and SR Release','Calcium Buffering','Calcium Efflux and SOCE','Cross-Bridge Cycle','10% Threshold')
    
    % labels_above_ten = cell(size(length(graph_names)));
    
    for i = 1:length(MorrisAnalysis.Results.MuStar(:,k))
        mu_star = MorrisAnalysis.Results.MuStar(i,k);

        if mu_star >= ten_percent
            muVec(i,k) = true;
            % mu_star = labels_above_ten{k};
        else
            muVec(i,k)=false;
        end
    end

 
    label_length = length(MorrisAnalysis.Results.VariableNames(1,:));
    label(1:label_length,n) = MorrisAnalysis.Results.VariableNames(1,:).';
    new_label(1:label_length,n) = MorrisAnalysis.Results.VariableNames(1,:).';
   
    % labels_above_ten = MorrisAnalysis.Results.VariableNames(1,:).' ;
    for j= 1:length(muVec(:,k))
        if muVec(j,k) ==1
            text(MorrisAnalysis.Results.MuStar(j,k), MorrisAnalysis.Results.Std(j,k), label(j,n));
        end
      
        if muVec(j,k) ==0
            new_label(j,n) = strrep(label(j,n),label(j,n),'0');
        end
    
    end

    prettyGraph
    hold off 
    % n= n+1;
end
