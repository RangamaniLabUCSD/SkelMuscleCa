clear 
load('MorrisResults_1e6NEW.mat')
graph_names = {"Steady State SR Calcium", "SS Voltage_{PM}", "SS Sodium Ion",...
    "SS Chlorine Ion","SS Myoplasmic Calcium","SS Potassium Ion","SS Force",...
    "Max Calcium Myo","Max Voltage","Max Force","Avg Myo Calcium","Avg Force",...
    "Avg Voltage", "AP Width"};
expt_names = repmat(graph_names,1,4)';
expt_names = [expt_names; "Objective Value"];

n = 1;
muVec = zeros(size(MorrisAnalysis.Results.MuStar));

ActionPotential = [1:6, 8:10, 13,14,16:30,33,36:40,45,76:82, 85:88];
CaInflux_SR_Release = [7,11,15,32,35,83,90:92];
CaBuffering = [46:54,58,68:72,74,75,96,97,98];
CaEfflux_SOCE = [12,31,34,41:44,89,95];
Crossbridge_Cycle = [55:57,59:67,73,84,93,94];
Diffusion = 98:105; % new parameters for diffusion between junctional and bulk compartments

figure
for k=1:length(expt_names)
    scatter(MorrisAnalysis.Results.MuStar(ActionPotential,k),...
        MorrisAnalysis.Results.Std(ActionPotential,k),'filled',...
        'MarkerFaceColor',[0.9290 0.6940 0.1250])
    hold on
    scatter(MorrisAnalysis.Results.MuStar(CaInflux_SR_Release,k),...
        MorrisAnalysis.Results.Std(CaInflux_SR_Release,k),'filled',...
        'MarkerFaceColor',[0.4660 0.6740 0.1880])
    scatter(MorrisAnalysis.Results.MuStar(CaBuffering,k),...
        MorrisAnalysis.Results.Std(CaBuffering,k),'filled',...
        'MarkerFaceColor',[0 0.4470 0.7410])
    scatter(MorrisAnalysis.Results.MuStar(CaEfflux_SOCE,k),...
        MorrisAnalysis.Results.Std(CaEfflux_SOCE,k),'filled',...
        'MarkerFaceColor',	[1 0 0])
    scatter(MorrisAnalysis.Results.MuStar(Crossbridge_Cycle,k),...
        MorrisAnalysis.Results.Std(Crossbridge_Cycle,k),'filled',...
        'MarkerFaceColor',[0.4940 0.1840 0.5560])
    scatter(MorrisAnalysis.Results.MuStar(Diffusion,k),...
        MorrisAnalysis.Results.Std(Diffusion,k),'filled',...
        'MarkerFaceColor','g')
    xlabel('mu*')
    ylabel('simga')
    yticklabels({})
    title(expt_names{k})
    hold on
    ten_percent = max(MorrisAnalysis.Results.MuStar(:,k))*0.1;
    xline(ten_percent , '--b' )
    legend('Action Potential', 'Calcium Influx and SR Release',...
        'Calcium Buffering','Calcium Efflux and SOCE','Cross-Bridge Cycle',...
        'Diffusion','10% Threshold')
    
    for i = 1:length(MorrisAnalysis.Results.MuStar(:,k))
        mu_star = MorrisAnalysis.Results.MuStar(i,k);
        if mu_star >= ten_percent
            muVec(i,k) = true;
        else
            muVec(i,k)=false;
        end
    end
    prettyGraph
    drawnow
    hold off 
end

% store names of significant parameters
anySig = any(muVec,2); % if any of the considered QOIs show high sensitivity
labels_above_ten = MorrisAnalysis.Results.VariableNames(1,:).';
for i = find(~anySig)'
    labels_above_ten{i} = '0';
end