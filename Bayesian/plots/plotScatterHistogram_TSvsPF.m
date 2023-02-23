function plotScatterHistogram_TSvsPF(BDscores, acceptedEvts, fileinfo, directory)


plotheight = 25; %% in cm
plotwidth  = 25;

nor = 1; %% number of events per row
noc = 3; %% number of events per column

mainPanelDimension = 2;
gapr = 3;     gapc = 4;

leftSpaceH = plotwidth - (noc*mainPanelDimension + (noc-1)*gapc);
leftSpaceV = plotheight - (nor*mainPanelDimension + (nor-1)*gapr);

leftedge = leftSpaceH/2;    rightedge = leftSpaceH/2;    topedge = leftSpaceV/2;    bottomedge = leftSpaceV/2;



fontsize = 5;
mainPanelsPos  = subplot_pos(plotwidth, plotheight, leftedge, rightedge, bottomedge, topedge, noc, nor, gapr, gapc); 


xAxis  = 'within-PBE time swap';
yAxis  = 'place field unit ID shuffle';



periods = fieldnames(BDscores);
BDReplayMetrics = {'weightedCorr';'replayScore'};

for jj = 1: length(BDReplayMetrics)
    
    
    figure;

    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [plotwidth plotheight]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

    
    for ii = 1: length(periods)
        

        period = periods{ii};
        BDReplayMetric = BDReplayMetrics{jj};

        xdata = BDscores.(period).data.wPBEtimeswap.(BDReplayMetric).prctilescore(acceptedEvts.(period), 1);
        ydata = BDscores.(period).data.unitIDshuffle.(BDReplayMetric).prctilescore(acceptedEvts.(period), 1);


        xHist = hist(xdata, 10);
        xHist = xHist/length(xdata);
        xHistChance = 0.1*ones(size(xHist));

        yHist = hist(ydata, 10);
        yHist = yHist/length(xdata);
        yHistChance = 0.1*ones(size(yHist));




        % edges = {0:10:100 0:10:100};
        % edges{1}(end) = 101;
        % edges{2}(end) = 101;
        % 
        % counts       = hist3([xdata ydata], 'Edges', edges);
        % counts(:,11) = []; 
        % counts(11,:) = [];
        % 
        % ratios = counts./sum(counts(:)); 
        % 
        % ratios = ratios'; 
        % 
        % diagHist = zeros(19, 1);
        % diagHist(10) = sum(diag(ratios));
        % for k = 1:9
        %     
        %     diagHist(10-k) = sum(diag(ratios, -k));
        %     diagHist(10+k) = sum(diag(ratios, k));
        %     
        % end




        %% scatter plot 


        if length(xdata) > 2000
            subSampleIdx = randperm(length(xdata), 2000);
            xdata2 = xdata(subSampleIdx);
            ydata2 = ydata(subSampleIdx);
        else
            xdata2 = xdata;
            ydata2 = ydata;
        end




        ax1 = axes('position',mainPanelsPos{1,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

        plot(xdata2, ydata2, '.', 'markersize', 3, 'color', 'k')
        xlim([0 100])
        ylim([0 100])

        set(gca, 'xtick', [0 100], 'ytick', [0 100], 'xticklabel', [], 'yticklabel', [], 'box' , 'off')

        annotation('line', [mainPanelsPos{1,ii}(1) sum(mainPanelsPos{1,ii}([1,3]))+0.04], [mainPanelsPos{1,ii}(2) sum(mainPanelsPos{1,ii}([2,4]))+0.04], 'color', 'r')


        title({period, sprintf('(n=%d)', length(xdata)), '', '', ''}, 'fontsize', 10)



        %% diametrical histogram  


        diamMainPanel = norm(mainPanelsPos{1,ii}([3 4]))/sin(pi/4);

        diamHistPos(1) = sum(mainPanelsPos{1,ii}([1 3]))-diamMainPanel/2 + mainPanelsPos{1,ii}(3)/2 + 0.005;
        diamHistPos(2) = sum(mainPanelsPos{1,ii}([2 4]))-diamMainPanel/2 + mainPanelsPos{1,ii}(4)/2 + 0.005; 
        diamHistPos(3) = diamMainPanel;
        diamHistPos(4) = diamMainPanel;


        ax2 = axes('position',diamHistPos,'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


        bins = -100:20:100;
        bins(end) = 100.5;

        xyDiff = histc(xdata - ydata, bins);
        xyDiff(end) = [];
        xyDiff = xyDiff/length(xdata);


        xyDiffChance = histc(xdata(randperm(length(xdata), length(xdata))) - ydata(randperm(length(ydata), length(ydata))), bins);
        xyDiffChance(end) = [];
        xyDiffChance = xyDiffChance/length(xdata);

        bar(xyDiff, 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)
        hold on
        bar(xyDiffChance, 'Facecolor', 'none', 'Edgecolor', [0.7 0.7 0.7])
        % bar(diagHist, 'Facecolor', 'b',% 'Edgecolor', 'none', 'facealpha', 0.5)

        ylim([0 1])
        xlim([0.5 10.5])
        camroll(-45)
        set(gca, 'visible', 'off')



        %% y histogram

        leftHistwidth  = 1/plotwidth;
        leftHistheight = mainPanelsPos{1,ii}(4);

        leftHistPos(1) = mainPanelsPos{1,ii}(1) - leftHistwidth - 0.005;
        leftHistPos(2) = mainPanelsPos{1,ii}(2);
        leftHistPos(3) = leftHistwidth;
        leftHistPos(4) = leftHistheight;

        ax3 = axes('position', leftHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');


        bar(fliplr(yHist), 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)
        hold on
        bar(fliplr(yHistChance), 'Facecolor', 'none', 'Edgecolor', [0.7 0.7 0.7])
        % bar(flipud(sum(ratios, 2)), 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)
        camroll(-90)

        xlim([0.5 10.5])
        ylim([0 max([xHist(:);yHist(:)])])
        set(gca, 'box', 'off', 'ytick', [], 'ycolor', 'none', 'xtick', [0.5 10.5], 'xticklabel', {'100', '0'})
        xlabel(yAxis)




        %% x histogram

        bottomHistwidth  = mainPanelsPos{1,ii}(3);
        bottomHistheight = 1/plotwidth;

        bottomHistPos(1) = mainPanelsPos{1,ii}(1);
        bottomHistPos(2) = mainPanelsPos{1,ii}(2) - bottomHistheight - 0.005;
        bottomHistPos(3) = bottomHistwidth;
        bottomHistPos(4) = bottomHistheight;


        ax4 = axes('position', bottomHistPos, 'XGrid', 'off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');

        bar(xHist, 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)
        hold on
        bar(fliplr(xHistChance), 'Facecolor', 'none', 'Edgecolor', [0.7 0.7 0.7])
        % bar(sum(ratios, 1), 'Facecolor', 'b', 'Edgecolor', 'none', 'facealpha', 0.5)
        xlim([0.5 10.5])
        ylim([0 max([xHist(:);yHist(:)])])
        set(gca, 'box', 'off', 'ytick', [], 'yColor', 'none', 'xtick', [0.5 10.5], 'xticklabel', {'0', '100'})
        xlabel(xAxis)

    end

    dim = [0.15 0.6 0.5 .1];
    annotation('textbox',dim,'String',{fileinfo.name; BDReplayMetric}, 'FitBoxToText','on', 'Interpreter', 'none');

    
    filename = fullfile(directory, ['tsVSPF_' BDReplayMetric]);
    print(gcf, filename, '-dpdf', '-r0')
    
  
end





