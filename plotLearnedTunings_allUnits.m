
% plot the learning tunings across time for example units in an example session

iSess = 1;
selectUnits = 1:30;

epochNames = {'pre'; 'run'; 'post'};

% KL divergence
% assemblyTuningPFKLdiv{iSess};
% 
% assemblyTuningKLdiv_time{iSess};
% 
% epochsKLdiv{iSess};
% epochsKLdiv_zscore{iSess};



%% plots

close all

cnctbinCenters{iSess} = [binCenters{iSess}.pre binCenters{iSess}.run binCenters{iSess}.post];


plotwidth  = 500;
plotheight = 850;
fontsize   = 4;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 30;
noc = 1;

leftmargin = 15;  rightmargin = 15;    topmargin = 60;    bottommargin = 60;

gapc = 15;
gapr = 5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);



% for each unit

plotwidth_unit  = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_unit = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_unit = 1;
noc_unit = 20;


leftmargin_unit = 0;  rightmargin_unit = 0;     bottommargin_unit = 0;    topmargin_unit = 0;
gapr_unit = 0;    gapc_unit = 10;

% ax_3 = zeros(nor*noc, 1);
% ax_4 = zeros(nor*noc, 1);


KLmin = inf;
KLmax = -inf;

for iUnit = 1:numel(selectUnits)
    currUnit = selectUnits(iUnit); 
    
    
    sub_pos_unit = subplot_pos(plotwidth_unit, plotheight_unit, leftmargin_unit, rightmargin_unit, bottommargin_unit, topmargin_unit, noc_unit, nor_unit, gapr_unit, gapc_unit); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1:noc_unit
        sub_pos_unit{ii}(1) = (sub_pos_unit{ii}(1) * plotwidth_unit + sub_pos{iUnit}(1) * plotwidth) / plotwidth;
        sub_pos_unit{ii}(2) = (sub_pos_unit{ii}(2) * plotheight_unit + sub_pos{iUnit}(2) * plotheight) / plotheight;
        sub_pos_unit{ii}(3) = sub_pos_unit{ii}(3) * plotwidth_unit / plotwidth;
        sub_pos_unit{ii}(4) = sub_pos_unit{ii}(4) * plotheight_unit / plotheight;
    end
    
    
    % plot the learned tunings in 15 minutes time windows
    
    position = [sub_pos_unit{1}(1:2) sum(sub_pos_unit{4}([1 3]))-sub_pos_unit{1}(1) sum(sub_pos_unit{4}([2 4]))-sub_pos_unit{1}(2)];
    axx = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    
    hold on

    currSpatialTuning = spikes_pooled{iSess}(currUnit).spatialTuning_smoothed.uni;
    nPosBins = numel(currSpatialTuning);
    
    learnedTuning = cell(3,1);
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time{iSess}.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters{iSess}.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); 

    end
    colormap(axx, 'jet')
    caxis([0 1])
    
    plot(0.1*currSpatialTuning + startT{iSess}.run , 1:nPosBins, 'linewidth', 1, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(gca);
    
    if iUnit == 1
        line(0.2*[0 2] + (endT{iSess}.post+startT{iSess}.post)/2, [yl(2)-10 yl(2)-10], 'linewidth', 1, 'color', [0 0 0 0.6])
        text((endT{iSess}.post+startT{iSess}.post)/2, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+8, ['unit ' num2str(currUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters{iSess}.post(end)])
    
    
    if ismember(iUnit, [30 60])
       for iEpoch = 1:3
           currEpoch      = epochNames{iEpoch};
           currbinCenters{iSess} = binCenters{iSess}.(currEpoch);
           text(median(currbinCenters{iSess}), yl(1)-50, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if ismember(iUnit, 1)
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, 30)
        xlabel('time (hour)')
    end
    
    
    if ~ismember(iUnit, 30)
       set(gca, 'XTickLabel', {}) 
    end
    
    if iUnit == 31
       h = colorbar;
       h.Label.String = 'normalized learned tuning';
       h.Location = 'northoutside'; 
       h.FontSize = fontsize;
       
       h.Position(2) = h.Position(2) + 0.02;
       h.Position(3) = h.Position(3)/3;
       h.Position(4) = h.Position(4)/2; 
    end
    
    
    hold off

    set(gca, 'position', position, 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', fontsize)

    
    
    
    % plot the cross-correlation matrix
    
    sub_pos_unit{5}(1) = sub_pos_unit{5}(1) - 0.2*sub_pos_unit{5}(3);
    sub_pos_unit{5}(3) = sub_pos_unit{5}(3)*1.4;
    sub_pos_unit{5}(2) = sub_pos_unit{5}(2) - 0.2*sub_pos_unit{5}(4);
    sub_pos_unit{5}(4) = sub_pos_unit{5}(4)*1.4;
    
    axx = axes('position', sub_pos_unit{5},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters{iSess}.(epochNames{iEpoch}), binCenters{iSess}.(epochNames{jEpoch}), corrMat)
 
        end
    end
    caxis([-1 1])
    colormap(axx, 'jet')
    

    ylim([0 binCenters{iSess}.post(end)])
    xlim([0 binCenters{iSess}.post(end)])
    
  
    set(gca, 'YTickLabel', {})
    set(gca, 'XTickLabel', [])
    
    if iUnit == 30
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};           
           
           h = text(median(binCenters{iSess}.(currEpoch)), binCenters{iSess}.pre(1)-50 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           set(h, 'rotation', 45)

       end
    end
    
    tickLabels = 0:2:cnctbinCenters{iSess}(end);
    xticks(tickLabels)
    yticks(tickLabels)
    
    axis square
    
    if iUnit == 1
        h = colorbar;
        h.Label.String = 'correlation';
        h.Location = 'northoutside';
       h.FontSize = fontsize;
       
       oldPos = h.Position;

    end
    set(gca, 'box', 'off', 'position', sub_pos_unit{5}, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    
    
    
    % plot the distribution of cross-correlations and its comparison with
    % unit-ID surrogates
    
    sub_pos_unit{6}(1) = sub_pos_unit{6}(1) + sub_pos_unit{6}(3)/2;
    sub_pos_unit{6}(3) = sub_pos_unit{6}(3)/2;
    
    
    position = [sub_pos_unit{6}(1:2) sum(sub_pos_unit{7}([1 3]))-sub_pos_unit{6}(1) sum(sub_pos_unit{7}([2 4]))-sub_pos_unit{6}(2)];
    
    ax_3(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    tt = 4;
    for iEpoch = 1:3
        
        for jEpoch = iEpoch:3
            
            if jEpoch == iEpoch
                idx = iEpoch;
            else
                idx = tt;
                tt = tt+1;
            end
            
            % data
            data = epochsXCorr{iSess}.data(iEpoch, jEpoch, currUnit, 1);
            
            scatter(idx-0.15, data, 1, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
            
            % ui surrogate
            ui = squeeze(epochsXCorr{iSess}.ui(iEpoch, jEpoch, :));
            [xData, yData] = myViolin(ui, [-1 1]);
            patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
            try
                pval = signrank(ui, data, 'tail', 'left');
                signn = significanceSign(pval);
                text(idx, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

            catch
            end

        end
    end
    
    xlim([0.5 6.5])

    tmp = gca;
    tmp.YGrid = 'on';
    
    
    if ismember(iUnit, [30 60])
        
        set(ax_3(iUnit), 'xticklabel', {'pre'; 'run'; 'post'; 'pre-run'; 'pre-post'; 'run-post'})
        
    else
        set(ax_3(iUnit), 'xticklabel', {})
    end
    
    set(ax_3(iUnit), 'fontsize', fontsize, 'position', position, 'box', 'off', 'xtick', 1:6, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax_3(iUnit), 45)
    
    
    if iUnit == 1
        ttl = title({'Pearson correlation'; ''; 'learned tuning stability'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = -0.1; 
        ttl.HorizontalAlignment = 'left';
    end
    
    
    

    
    % plot the PF Pearson correlation of 15minutes-time window learned tunings
    
    data = assemblyTuningPFcorr_time{iSess}.pre.data(currUnit, :);
    if ~isempty(data(~isnan(data)))
    
        ax_4(iUnit) = axes('position', sub_pos_unit{8},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        hold on

        for iEpoch = 1:3
            idx = iEpoch;
            currEpoch = epochNames{iEpoch};
            
            % data
            data = assemblyTuningPFcorr_time{iSess}.(currEpoch).data(currUnit, :);

            if iEpoch == 2
                scatter((idx-0.15)*ones(numel(data), 1), data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
            else
                [xData, yData] = myViolin(data, [-1 1]);

                patch(idx-0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none') 
            end

            % shuffle

            ui = assemblyTuningPFcorr_time{iSess}.(currEpoch).ui;
            ui = ui(:);
            [xData, yData] = myViolin(ui, [-1 1]);

            patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')


            try
               
                
%                 pval1 = numel(find(data > ui))/numel(ui);
%                 pval2 = numel(find(data < ui))/numel(ui);
%                 
                
                pval = ranksum(data, ui);
                signn = significanceSign(pval);
                text(idx, 1.1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

                zeroMedianPval = signrank(data);
                signn = significanceSign2(zeroMedianPval);
                text(idx, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

            catch

            end

        end
    
    
        xlim([0 4])

        set(ax_4(iUnit), 'fontsize', fontsize, 'position', sub_pos_unit{8}, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])


        if ismember(iUnit, [30 60])

            set(ax_4(iUnit), 'xticklabel', {'pre'; 'run'; 'post'})

        else
            set(ax_4(iUnit), 'xticklabel', {})
        end

        xtickangle(ax_4(iUnit), 45)

        tmp = gca;
        tmp.YGrid = 'on';

        if iUnit == 1
            ttl = title({'PF matching score'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = -0.1; 
            ttl.HorizontalAlignment = 'left';
        end

    end
    
    
    
    
    
    % plot the distribution of place field correlation for the learned tunings calculated during the whole period
    
    position = sub_pos_unit{9};
    position(1) = position(1)-0.2*position(3);
    
    
    if ~isnan(assemblyTuningPFcorr{iSess}.pre.data(currUnit))
        ax_6(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        hold on

        for iEpoch = 1:3
            idx = iEpoch;
            currEpoch = epochNames{iEpoch};

            % data

            data = assemblyTuningPFcorr{iSess}.(currEpoch).data(currUnit);

            h1 = scatter(idx-0.15, data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);


            % shuffle

            ui = assemblyTuningPFcorr{iSess}.(currEpoch).ui;
            ui = ui(:);
            [xData, yData] = myViolin(ui, [-1 1]);

            h2 = patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')

            try
                pval = signrank(ui, data);
                signn = significanceSign(pval);
                text(idx, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

            catch

            end
        end
    
        xlim([0 4])

    %     if iUnit == 1
    %        legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northeastout', 'box', 'off')
    %     end

        set(ax_6(iUnit), 'fontsize', fontsize, 'position', position, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])


        if ismember(iUnit, [30 60])

            set(ax_6(iUnit), 'xticklabel', {'pre'; 'run'; 'post'})

        else
            set(ax_6(iUnit), 'xticklabel', {})
        end

        xtickangle(ax_6(iUnit), 45)

        tmp = gca;
        tmp.YGrid = 'on';

        if iUnit == 1
            ttl = title('entire', 'fontsize', fontsize, 'fontweight', 'normal'); 
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = 0; 
            ttl.HorizontalAlignment = 'left';
        end
    
    end
    

    
    
    % plot KL divergence stuff
    
    % matrix
    
    sub_pos_unit{10}(1) = sub_pos_unit{10}(1) - 0.2*sub_pos_unit{10}(3);
    sub_pos_unit{10}(3) = sub_pos_unit{10}(3)*1.4;
    sub_pos_unit{10}(2) = sub_pos_unit{10}(2) - 0.2*sub_pos_unit{10}(4);
    sub_pos_unit{10}(4) = sub_pos_unit{10}(4)*1.4;
    
    ax_7(iUnit) = axes('position', sub_pos_unit{10},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            KLdiv = calKLDivergence(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            KLdiv(isnan(KLdiv)) = inf;
            imagesc(binCenters{iSess}.(epochNames{iEpoch}), binCenters{iSess}.(epochNames{jEpoch}), KLdiv)
 
        end
    end
    
    KLmin = min(KLmin, min(KLdiv(:)));
    KLmax = max(KLmax, max(KLdiv(:)));
    
    caxis([0 1])
    colormap(ax_7(iUnit), flipud(jet))
    

    ylim([0 binCenters{iSess}.post(end)])
    xlim([0 binCenters{iSess}.post(end)])
    
  
    set(gca, 'YTickLabel', {})
    set(gca, 'XTickLabel', [])
    
    if iUnit == 30
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};           
           
           h = text(median(binCenters{iSess}.(currEpoch)), binCenters{iSess}.pre(1)-50 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           set(h, 'rotation', 45)

       end
    end
    
    tickLabels = 0:2:cnctbinCenters{iSess}(end);
    xticks(tickLabels)
    yticks(tickLabels)
    
    axis square
    
    if iUnit == 1
        h = colorbar;
        h.Label.String = 'KL divergence';
        h.Location = 'northoutside';
        h.FontSize = fontsize;
        h.Ticks = [0 0.5 1]; 
        h.TickLabels = {'0'; '0.5'; '>1'};
       
       oldPos = h.Position;

    end
    set(gca, 'box', 'off', 'position', sub_pos_unit{10}, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    
    
    
    
    % plot the distribution of KL divergence and its comparison with
    % unit-ID surrogates
    
    sub_pos_unit{11}(1) = sub_pos_unit{11}(1) + sub_pos_unit{11}(3)/2;
    sub_pos_unit{11}(3) = sub_pos_unit{11}(3)/2;
    
    
    position = [sub_pos_unit{11}(1:2) sum(sub_pos_unit{12}([1 3]))-sub_pos_unit{11}(1) sum(sub_pos_unit{12}([2 4]))-sub_pos_unit{11}(2)];
    
    ax_8(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    tt = 4;
    for iEpoch = 1:3
        
        for jEpoch = iEpoch:3
            
            if jEpoch == iEpoch
                idx = iEpoch;
            else
                idx = tt;
                tt = tt+1;
            end
            
            
            % data
            data = epochsKLdiv{iSess}.data(iEpoch, jEpoch, currUnit);
            scatter(idx-0.15, data, 1, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
            
            % ui surrogate
            ui = squeeze(epochsKLdiv{iSess}.ui(iEpoch, jEpoch, :));
            [xData, yData] = myViolin(ui, [0 3]);
            patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
            try
                pval = ranksum(ui, data, 'tail', 'right');
                signn = significanceSign(pval);
                text(idx, 3, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')
                
            catch
            end

        end
    end
    
    xlim([0.5 6.5])
    ylim([0 3])
    
    tmp = gca;
    tmp.YGrid = 'on';
    
    
    if ismember(iUnit, [30 60])
        
        set(ax_8(iUnit), 'xticklabel', {'pre'; 'run'; 'post'; 'pre-run'; 'pre-post'; 'run-post'})
        
    else
        set(ax_8(iUnit), 'xticklabel', {})
    end
    
    set(ax_8(iUnit), 'fontsize', fontsize, 'position', position, 'box', 'off', 'xtick', 1:6, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax_8(iUnit), 45)
    
    
    if iUnit == 1
        ttl = title({'KL divergence'; '';'learned tuning stability'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = -0.1; 
        ttl.HorizontalAlignment = 'left';
    end
    
    
    
    
    % plot the PF KL divergence of 15minutes-time window learned tunings
    
    data = assemblyTuningPFKLdiv_time{iSess}.pre.data(currUnit, :);
    if ~isempty(data(~isnan(data)))
    
        ax_9(iUnit) = axes('position', sub_pos_unit{13},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        hold on

        
        
        for iEpoch = 1:3
            idx = iEpoch;
            currEpoch = epochNames{iEpoch};
            
            
            
            
            % data
            data = assemblyTuningPFKLdiv_time{iSess}.(currEpoch).data(currUnit, :);
            
            if iEpoch == 2
                scatter((idx-0.15)*ones(numel(data), 1), data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
            else
                [xData, yData] = myViolin(data, [0 3]);

                patch(idx-0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none') 
            
            end

            
            
            
            % shuffle
            ui = assemblyTuningPFKLdiv_time{iSess}.(currEpoch).ui;
            ui = ui(:);
            
            [xData, yData] = myViolin(ui, [0 3]);

            patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')


            try
                pval = ranksum(data, ui, 'tail', 'left');
                signn = significanceSign(pval);
                text(idx, 3, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

            catch

            end
            
            

        end
    
        ylim([0 3])

        xlim([0 4])

        set(ax_9(iUnit), 'fontsize', fontsize, 'position', sub_pos_unit{13}, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])


        if ismember(iUnit, [30 60])

            set(ax_9(iUnit), 'xticklabel', {'pre'; 'run'; 'post'})

        else
            set(ax_9(iUnit), 'xticklabel', {})
        end

        xtickangle(ax_9(iUnit), 45)

        tmp = gca;
        tmp.YGrid = 'on';

        if iUnit == 1
            ttl = title({'PF matching score'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = -0.1; 
            ttl.HorizontalAlignment = 'left';
        end

    end
    
    
    
    
    % plot the distribution of KL divergence for learned tuning calculated during the entire period
    
    position = sub_pos_unit{14};
    position(1) = position(1)-0.2*position(3);
    
    
    if ~isnan(assemblyTuningPFKLdiv{iSess}.pre.data(currUnit))
        ax_10(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        hold on

        for iEpoch = 1:3
            idx = iEpoch;
            currEpoch = epochNames{iEpoch};

            % data

            data = assemblyTuningPFKLdiv{iSess}.(currEpoch).data(currUnit);

            h1 = scatter(idx-0.15, data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);


            % shuffle

            ui = assemblyTuningPFKLdiv{iSess}.(currEpoch).ui;
            ui = ui(:);
            [xData, yData] = myViolin(ui, [0 3]);

            h2 = patch(idx+0.15+xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')

            try
                pval = signrank(ui, data, 'tail', 'right');
                signn = significanceSign(pval);
                text(idx, 3, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')

            catch

            end
        end
        ylim([0 3])

        xlim([0 4])

    %     if iUnit == 1
    %        legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northeastout', 'box', 'off')
    %     end

        set(ax_10(iUnit), 'fontsize', fontsize, 'position', position, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])


        if ismember(iUnit, [30 60])

            set(ax_10(iUnit), 'xticklabel', {'pre'; 'run'; 'post'})

        else
            set(ax_10(iUnit), 'xticklabel', {})
        end

        xtickangle(ax_10(iUnit), 45)

        tmp = gca;
        tmp.YGrid = 'on';

        if iUnit == 1
            ttl = title('entire', 'fontsize', fontsize, 'fontweight', 'normal'); 
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = 0; 
            ttl.HorizontalAlignment = 'left';
        end
    
    end
    


end

linkaxes(ax_3, 'y')
linkaxes(ax_4, 'y')
% linkaxes(ax_5, 'y')
linkaxes(ax_6, 'y')

% 
% for ii = 1:numel(ax_7)
%     caxis(ax_7(ii), [0 5])
% end

    
%% sub-functions



function KLdiv = calKLDivergence(learnedTunings, spatialTunings)


% the individual learned tunnigs and spatial tunings should be in columns


% this function works well when when the learnedTunings matrix contains multipels instance of learned tunings (for examples if they belong to different time windows or different instances of unit identity shuffle)

nPosBins = size(spatialTunings, 1);

spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [nPosBins 1]);
spatialTunings = spatialTunings + eps;

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [nPosBins 1]);
learnedTunings = learnedTunings + eps;


nSTs = size(spatialTunings, 2);
nLTs = size(learnedTunings, 2); % number of learned tunings


KLdiv = nan(nLTs, nSTs);

for ii = 1:nSTs
    
    currSpatialTuning = spatialTunings(:, ii);

    spatialTuningTerm = repmat(currSpatialTuning, [1 nLTs]);
    KLdiv(:, ii) = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);

end

end



function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end


function [xData, yData] = myViolin(data, theRange)

% binEdges = linspace(min(data)-0.1*range(data), max(data)+0.1*range(data), 50);

binEdges = linspace(theRange(1), theRange(2), 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


if size(count, 1) > size(count, 2)
    count = count';
end

gw = gausswindow(2,6);
count = conv(count, gw, 'same');
count = count/sum(count);

xData = [count -fliplr(count)];
yData = [binCenters fliplr(binCenters)];

end

function signn = significanceSign(pval)

signn = '';

if pval < 0.001
    signn = '***';
elseif pval < 0.01
    signn = '**';
elseif pval < 0.05
    signn = '*';
end

end

function signn = significanceSign2(pval)

signn = '';

if pval < 0.001
    signn = '###';
elseif pval < 0.01
    signn = '##';
elseif pval < 0.05
    signn = '#';
end

end