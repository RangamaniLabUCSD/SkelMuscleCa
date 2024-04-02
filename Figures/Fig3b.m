function Fig3b(X, Y1, Y2)
% Zoomed in plot of V_SL and Myoplasmic Ca^2+
%  X:  Time (s)
%  Y1: V_SL (mV)
%  Y2: Myoplasmic Ca^2+ (uM)

figure1 = figure;
axes1 = axes('Parent',figure1,'FontSize',14,'Box','off','FontSmoothing','on'); 

yyaxis(axes1,'left');
set(axes1,'YColor',[0 0 0]);
plot(X,Y1,'LineWidth',2,'Color',[0.188235294117647 0.819607843137255 0.788235294117647]);
ylabel('V_{SL} (mV)','FontWeight','bold','Color',[0.188235294117647 0.819607843137255 0.788235294117647],'FontSmoothing','on');

yyaxis(axes1,'right');
set(axes1,'YColor',[0 0 0]);
plot(X,Y2,'LineWidth',2,'Color',[0.635294117647059 0.0784313725490196 0.184313725490196]); 
ylabel('[Ca^{2+}]_{myo} (\muM)','FontWeight','bold','FontSize',16,'Color',[0.635294117647059 0.0784313725490196 0.184313725490196],'FontSmoothing','on');

xlabel('Time (s)','FontWeight','bold','FontSize',16,'FontSmoothing','on');
xlim(axes1,[0.05 0.15]);

