clear 
load('MorrisResults9-16.mat')

% QOI(i,:) = [yInf(2), yInf(5), yInf(6), yInf(7),yInf(8), yInf(13), yInf(23), MaxCaF,MaxVF, MaxPost, AvgF, AvgPost, AvgVolt, VoltWidth];

% graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion","SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force", "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force","Avg Voltage", "AP Width"};
graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion","SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force", "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force","Avg Voltage", "AP Width"};
expt_names = repmat(graph_names,4)';

n = 1;
muVec = zeros(size(MorrisAnalysis.Results.MuStar));

for k= 57 %1:length(graph_names)
    figure
    scatter(MorrisAnalysis.Results.MuStar(:,k),MorrisAnalysis.Results.Std(:,k),'filled','MarkerFaceColor',[0.64,0.08,0.18])
    xlabel('mu*')
    ylabel('simga')
    % title(graph_names{k})
    hold on
    ten_percent = max(MorrisAnalysis.Results.MuStar(:,k))*0.1;
    xline(ten_percent , '--b' )


    
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
