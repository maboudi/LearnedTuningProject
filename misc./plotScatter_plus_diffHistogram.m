    `1function plotScatter_plus_diffHistogram(xdata, ydata)

plotheight = 250; %% in points
plotwidth  = 250;

nor = 1; %% number of events per row
noc = 1; %% number of events per column

mainPanelDimension = 80;
gapr = 10;     gapc = 10;

leftSpaceH = plotwidth - (noc*mainPanelDimension + (noc-1)*gapc);
leftSpaceV = plotheight - (nor*mainPanelDimension + (nor-1)*gapr);

leftedge = leftSpaceH/2;    rightedge = leftSpaceH/2;    topedge = leftSpaceV/2;    bottomedge = leftSpaceV/2;



fontsize = 6;
mainPanelsPos  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); 


figure;

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



%% scatter plot 

ax1 = axes('position', mainPanelsPos{1,1},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

scatter(xdata, ydata, 2, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5)

xlim([-1 1])
ylim([-1 1])

set(gca, 'xtick', [-1 0 1], 'ytick', [-1 0 1], 'box' , 'off')

annotation('line', [mainPanelsPos{1,1}(1) sum(mainPanelsPos{1,1}([1,3]))+0.04], [mainPanelsPos{1,1}(2) sum(mainPanelsPos{1,1}([2,4]))+0.04], 'color', 'r')



%% diametrical histogram  


diamMainPanel = norm(mainPanelsPos{1,1}([3 4]))/sin(pi/4);

diamHistPos(1) = sum(mainPanelsPos{1,1}([1 3]))-diamMainPanel/2 + mainPanelsPos{1,1}(3)/2 + 0.005;
diamHistPos(2) = sum(mainPanelsPos{1,1}([2 4]))-diamMainPanel/2 + mainPanelsPos{1,1}(4)/2 + 0.005; 
diamHistPos(3) = diamMainPanel;
diamHistPos(4) = diamMainPanel;


ax2 = axes('position',diamHistPos,'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


bins = -2:0.2:2;

xyDiff = histc(xdata - ydata, bins);
xyDiff(end) = [];
xyDiff = xyDiff/length(xdata);

bar(xyDiff, 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)


ylim([0 1])
xlim([0.5 20.5])
camroll(-45)
set(gca, 'visible', 'off')

