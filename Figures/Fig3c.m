function Fig3c(X, Y1, Y2)
% Function to plot the opening probability of Orai1 and SR Ca^2+
%  X:  Time (s)
%  Y1:  Opening probability of Orai1 channel
%  Y2:  SR Ca^2+ concentration (uM)

figure; %('OuterPosition',[1289 392 1086 872]);

axes1 = axes('Position',[0.109439252336449 0.111283697047497 0.775 0.815]);
hold(axes1,'on');

yyaxis left %(axes1,'left');
plot(X,Y1,'Color',[0 0.447058823529412 0.741176470588235]);
ylabel('p_{Orai1}^{open}','FontSize',16,'Color',[0 0.447058823529412 0.741176470588235]);


yyaxis right %(axes1,'right');
plot(X,Y2)%,'Color',[0.494117647058824 0.184313725490196 0.556862745098039]);
ylabel('[Ca^{2+}]_{SR}','FontSize',16) %,'Color',[0.494117647058824 0.184313725490196 0.556862745098039]);


xlabel({'Time (s)'},'FontWeight','bold');
box(axes1,'on');
set(axes1,'FontSize',14);
hold(axes1,'off');

