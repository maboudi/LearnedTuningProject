clear
clc


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);

epochNames = {'pre'; 'run'; 'post'};





sessionNumber = 3;
sessionName   = rr(sessionNumber+2).name;

basePath      = fullfile(parentDir, sessionName);




% spikes and behavior data

load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
behavior = fileInfo.behavior;
behavior.time = behavior.time/3600;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

nUnits = numel(spikes);


% bidirectional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end



% load learned tunings
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', 'assemblyTunings_time', 'activeUnits');
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTunings', 'assemblyTuningPFcorr');



for iEpoch = 1:numel(epochNames)
    
    currEpoch = epochNames{iEpoch};
    
    binCenters.(currEpoch) = binCenters.(currEpoch)/3600;

end


 
% % For the unit-ID shuffle do independent randomization of unit-ID for each
% % time window
% 
% 
% for iEpoch = 1:3
%     currEPoch = epochNames{iEpoch};
%     
%     allLTs = assemblyTunings_time.(currEpoch).ui;
%     
%     nTimeBins = size(allLTs, 3);
%     nInstances = size(allLTs, 4);
%     allLTs2 = nan(size(allLTs));
%     
%     for iTimeBin = 1:nTimeBins
%         allLTs2(:, :, iTimeBin, :) = allLTs(:, :, iTimeBin, randperm(nInstances));
%     end
%     
%     assemblyTunings_time.(currEpoch).ui = allLTs2;
%     
% end



% alternative unit-ID shuffle by shuffling the unit identities of the learned
% tunings rather than the PBE spikings

okUnits = intersect(intersect(activeUnits.pre, activeUnits.post), activeUnits.run);

for iEpoch = 1:3
    epochName = epochNames{iEpoch};
    
    
    currData = assemblyTunings_time.(epochName).data;
    nTimeBins = size(currData, 3);
    
    assemblyTunings_time.(epochName).ui = nan([size(currData), 100]);
    for inst = 1:100
        
        for iTimeBin = 1:nTimeBins

            for iUnit = 1:numel(okUnits)

                currUnit = okUnits(iUnit);
                otherUnits = setdiff(okUnits, currUnit);
                assemblyTunings_time.(epochName).ui(currUnit, :, iTimeBin, inst) = currData(otherUnits(randperm(numel(otherUnits), 1)), :, randperm(nTimeBins, 1));
            end
        end
    end
end




% remove the time windows with overlap
    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap.(epochName) = binCenters.(epochName)(1:3:end);
    
    assemblyTunings_time_nonOverlap.(epochName).data = assemblyTunings_time.(epochName).data(:, :, 1:3:end);
    assemblyTunings_time_nonOverlap.(epochName).ui  = assemblyTunings_time.(epochName).ui(:, :, 1:3:end, :); 

end



% truncate the long recording sessions to the first 4 hours of post sleep
idx = find(binCenters.post < endT.post);
binCenters.post = binCenters.post(idx);
assemblyTunings_time.post.data = assemblyTunings_time.post.data(:, :, idx);
assemblyTunings_time.post.ui   = assemblyTunings_time.post.ui(:, :, idx, :);


idx = find(binCenters_nonOverlap.post < endT.post);
binCenters_nonOverlap.post = binCenters_nonOverlap.post(idx);
assemblyTunings_time_nonOverlap.post.data = assemblyTunings_time_nonOverlap.post.data(:, :, idx);
assemblyTunings_time_nonOverlap.post.ui   = assemblyTunings_time_nonOverlap.post.ui(:, :, idx, :);






% calculate the PF correlations for data and ui



% dataTypes = {'data'; 'ui'};
% for iDataType = 1:2
%     currDataType = dataTypes{iDataType};
%     
%     for iEpoch = 1:3
%         currEpoch = epochNames{iEpoch};
%         
%         allLTs = assemblyTunings_time_nonOverlap.(currEpoch).(currDataType);
%         
%         nTimeBins = size(allLTs, 3);
%         nInstances = size(allLTs, 4); 
%         
%         assemblyTuningPFcorr_time.(currEpoch).(currDataType) = nan(nUnits, nTimeBins, nInstances);
%         
%         for ins = 1:nInstances
%         
%             for iTimeBin = 1:nTimeBins 
%                 currLT = allLTs(:, :, iTimeBin, ins);
%                 corrMat = corr(currLT', spatialTunings_merge', 'type', 'pearson');
%                 assemblyTuningPFcorr_time.(currEpoch).(currDataType)(:, iTimeBin, ins) = diag(corrMat);
%             end
%         end
% 
%     end
%     
% end

% % % second thought

%% data

dataTypes = {'data'; 'ui'};
for iDataType = 1:2
    currDataType = dataTypes{iDataType};
    
    for iEpoch = 1:3
        currEpoch = epochNames{iEpoch};
        
        allLTs = assemblyTunings_time_nonOverlap.(currEpoch).(currDataType);
        
        nTimeBins = size(allLTs, 3);
        nInstances = size(allLTs, 4); 
        
        assemblyTuningPFcorr_time.(currEpoch).(currDataType) = nan(nUnits, nTimeBins, nInstances);
        
        for ins = 1:nInstances
        
            for iTimeBin = 1:nTimeBins 
                currLT = allLTs(okUnits, :, iTimeBin, ins);
                currPF = spatialTunings_merge(okUnits, :);
                
                if iDataType == 1
                    corrMat = corr(currLT', currPF', 'type', 'pearson');
                elseif iDataType == 2
                    corrMat = corr(currLT', currPF(randperm(numel(okUnits)), :)', 'type', 'pearson');
                end
                assemblyTuningPFcorr_time.(currEpoch).(currDataType)(okUnits, iTimeBin, ins) = diag(corrMat);
            end
        end

    end
    
end




%% for assemblyTunings_time calculate the correlation matrix of learned tunings

dataTypes = {'data'; 'ui'};
epochNames = {'pre'; 'run'; 'post'};


for idataType = 1:2
    currDataType = dataTypes{idataType}; 
    
    nInstances = size(assemblyTunings_time_nonOverlap.pre.(currDataType), 4);       
    epochsXCorr.(currDataType) = cell(3,3, nUnits, nInstances);  
    epochsXCorr_entire.(currDataType) = nan(3,3, nUnits, nInstances);  


    for iUnit = 1:nUnits

        for inst = 1:nInstances
            for iEpoch = 1:3
                tuning1 = squeeze(assemblyTunings_time_nonOverlap.(epochNames{iEpoch}).(currDataType)(iUnit, :, :, inst));

                for jEpoch = iEpoch:3
                    tuning2 = squeeze(assemblyTunings_time_nonOverlap.(epochNames{jEpoch}).(currDataType)(iUnit, :, :, inst));

                    corrMat = corr(tuning1, tuning2);

                    if jEpoch == iEpoch
                        
                        allOffDiag = [];
                        for ii = 1:size(corrMat, 1)-1
                            allOffDiag = [allOffDiag; diag(corrMat, ii)];
                        end
                        corrMat = allOffDiag;
                        
%                         corrMat = corrMat - eye(size(corrMat, 1)).*corrMat;
                    end
                    epochsXCorr.(currDataType){iEpoch, jEpoch, iUnit, inst} = corrMat(:);
                    
                end
            end
        end
        
        
        
        %the same but for the assembly tunings based on the entire period

        for inst = 1:nInstances
            for iEpoch = 1:3
                tuning1 = squeeze(assemblyTunings.(epochNames{iEpoch}).(currDataType)(iUnit, :, inst));
                
                for jEpoch = iEpoch+1:3
                    tuning2 = squeeze(assemblyTunings.(epochNames{jEpoch}).(currDataType)(iUnit, :, inst));
                    
                    epochsXCorr_entire.(currDataType)(iEpoch, jEpoch, iUnit, inst) = corr(tuning1', tuning2');
                end
            end
        end
        

        
    end    
end

% epochsXCorr_zscore = (epochsXCorr.data) - mean(abs(epochsXCorr.ui), 4)./std(abs(epochsXCorr.ui), [], 4);





%% plot the learning tunings across time for example units in an example session


close all


okUnits = intersect(intersect(activeUnits.pre, activeUnits.post), activeUnits.run);
selectUnits = okUnits(1:end);
% selectUnits = [60 19 40 45];


cnctBinCenters = [binCenters.pre binCenters.run binCenters.post];


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
noc = 2;

leftmargin = 15;  rightmargin = 15;    topmargin = 60;    bottommargin = 60;

gapc = 15;
gapr = 5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);




% for each unit

plotwidth_unit = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_unit = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_unit = 1;
noc_unit = 10;

leftmargin_unit = 0;  rightmargin_unit = 0;     bottommargin_unit = 0;    topmargin_unit = 0;
gapr_unit = 0;    gapc_unit = 10;

ax_3 = zeros(nor*noc, 1);
ax_4 = zeros(nor*noc, 1);

for iUnit = 1:numel(selectUnits)
    currUnit = selectUnits(iUnit); 
    
    
    sub_pos_unit = subplot_pos(plotwidth_unit, plotheight_unit, leftmargin_unit, rightmargin_unit, bottommargin_unit, topmargin_unit, noc_unit, nor_unit, gapr_unit, gapc_unit); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1:noc_unit
        sub_pos_unit{ii}(1) = (sub_pos_unit{ii}(1) * plotwidth_unit + sub_pos{iUnit}(1) * plotwidth) / plotwidth;
        sub_pos_unit{ii}(2) = (sub_pos_unit{ii}(2) * plotheight_unit + sub_pos{iUnit}(2) * plotheight) / plotheight;
        sub_pos_unit{ii}(3) = sub_pos_unit{ii}(3) * plotwidth_unit / plotwidth;
        sub_pos_unit{ii}(4) = sub_pos_unit{ii}(4) * plotheight_unit / plotheight;
    end
    
    
    % plot the learned tunings across time
    
    position = [sub_pos_unit{1}(1:2) sum(sub_pos_unit{4}([1 3]))-sub_pos_unit{1}(1) sum(sub_pos_unit{4}([2 4]))-sub_pos_unit{1}(2)];

    axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    
    hold on

    currSpatialTuning = spikes(currUnit).spatialTuning_smoothed.uni;
    
    learnedTuning = cell(3,1);
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

    end
    caxis([0 1])
    
    plot(0.1*currSpatialTuning + startT.run , 1:nPosBins, 'linewidth', 1, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(gca);
    
    if iUnit == 1
        line(0.2*[0 2] + (endT.post+startT.post)/2, [yl(2)-10 yl(2)-10], 'linewidth', 1, 'color', [0 0 0 0.6])
        text((endT.post+startT.post)/2, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+8, ['unit ' num2str(currUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters.post(end)])
    
    
    if ismember(iUnit, [30 60])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters.(currEpoch);
           text(median(currBinCenters), yl(1)-50, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if ismember(iUnit, [1])
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, [30 60])
        xlabel('time (hour)')
    end
    
    
    if ~ismember(iUnit, [30 60])
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

    
    
    
    % plot the cross-correlation matrices
    
    sub_pos_unit{5}(1) = sub_pos_unit{5}(1) - 0.2*sub_pos_unit{5}(3);
    sub_pos_unit{5}(3) = sub_pos_unit{5}(3)*1.4;
    sub_pos_unit{5}(2) = sub_pos_unit{5}(2) - 0.2*sub_pos_unit{5}(4);
    sub_pos_unit{5}(4) = sub_pos_unit{5}(4)*1.4;
    
    axes('position', sub_pos_unit{5},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters.(epochNames{iEpoch}), binCenters.(epochNames{jEpoch}), corrMat)
 
        end
    end
    caxis([-1 1])
    colormap('jet')
    

    ylim([0 binCenters.post(end)])
    xlim([0 binCenters.post(end)])
    
  
    set(gca, 'YTickLabel', {})
    set(gca, 'XTickLabel', [])
    
    if iUnit == 1
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};           
           
           h = text(median(binCenters.(currEpoch)), binCenters.post(1)+5 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           set(h, 'rotation', 45)

       end
    end
    
    tickLabels = 0:2:cnctBinCenters(end);
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
            
            data = squeeze(epochsXCorr.data{iEpoch, jEpoch, currUnit, 1});
            
            if numel(data) > 1
            
                [xData, yData] = myViolin(data);

                h1 = patch(idx-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
            elseif numel(data) == 1
                scatter(idx-0.15, data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
            end
            
            
            % ui surrogate
            
            ui = cell2mat(squeeze(epochsXCorr.ui(iEpoch, jEpoch, currUnit, :)));
            [xData, yData] = myViolin(ui);
            h2 = patch(idx+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            
            scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
            try
                pval = ranksum(data, ui);
                signn = significanceSign(pval);
                text(idx, 1.1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')
                
                zeroMedianPval = signrank(data);
                signn = significanceSign2(zeroMedianPval);
                text(idx, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')
                
                
            catch
            end

        end
    end
    
    xlim([0.5 6.5])

    
%     if iUnit == 31
%        legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northeastout', 'box', 'off')
%     end
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
        ttl = title({'learned tuning consistency'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = -0.1; 
        ttl.HorizontalAlignment = 'left';
    end
    
    
    % plot the correlation with PF in each epoch for each unit
    
    ax_4(iUnit) = axes('position', sub_pos_unit{9},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    for iEpoch = 1:3
        idx = iEpoch;
        currEpoch = epochNames{iEpoch};
        % data
  
        data = assemblyTuningPFcorr_time.(currEpoch).data(currUnit, :);
        
        
        if iEpoch == 2
            scatter((idx-0.15)*ones(numel(data), 1), data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
        else
            [xData, yData] = myViolin(data);

            patch(idx-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none') 
        end

        
        % shuffle
        
        ui = assemblyTuningPFcorr_time.(currEpoch).ui(currUnit, :, :);
        ui = ui(:);
        [xData, yData] = myViolin(ui);

        patch(idx+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
        
        try
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

    set(ax_4(iUnit), 'fontsize', fontsize, 'position', sub_pos_unit{9}, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])

    
    if ismember(iUnit, [30 60])
        
        set(ax_4(iUnit), 'xticklabel', {'pre'; 'run'; 'post'})
        
    else
        set(ax_4(iUnit), 'xticklabel', {})
    end
    
    xtickangle(ax_4(iUnit), 45)
    
    tmp = gca;
    tmp.YGrid = 'on';
    
    if iUnit == 1
        ttl = title({'Place Field correlation'; ''; '15min-wins'}, 'fontsize', fontsize, 'fontweight', 'normal');
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = -0.1; 
        ttl.HorizontalAlignment = 'left';
    end

    
    
    
    % plot the distribution of consistency between epochs for the assembly tuning calculated during the whole period
    
    position = sub_pos_unit{8};
    position(1) = position(1) - 0.2*position(3);
    
    ax_5(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    tt = 1;
    for iEpoch = 1:3
        
        for jEpoch = iEpoch+1:3
            
            
            % data
            data = squeeze(epochsXCorr_entire.data(iEpoch, jEpoch, currUnit, 1));            
            scatter(tt-0.15, data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

            
            % ui surrogate
            
            ui = squeeze(epochsXCorr_entire.ui(iEpoch, jEpoch, currUnit, :));
            [xData, yData] = myViolin(ui);
            patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            
            scatter(tt+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
            try
                pval = signrank(ui, data);
                signn = significanceSign(pval);
                text(tt, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')
                
            catch
            end
            
            tt = tt +1;

        end
    end
    
    xlim([0.5 3.5])
    
    tmp = gca;
    tmp.YGrid = 'on';
    
    
    if ismember(iUnit, [30 60])
        
        set(ax_5(iUnit), 'xticklabel', {'pre-run'; 'pre-post'; 'run-post'})
        
    else
        set(ax_5(iUnit), 'xticklabel', {})
    end
    
    set(ax_5(iUnit), 'fontsize', fontsize, 'position', position, 'box', 'off', 'xtick', 1:3, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax_5(iUnit), 45)
    
    
    if iUnit == 1
        ttl = title('entire', 'fontsize', fontsize, 'fontweight', 'normal');
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = 0; 
        ttl.HorizontalAlignment = 'left';
    end
    

    
    % plot the distribution of place field correlation for the assembly tuning calculated during the whole period
    
    position = sub_pos_unit{10};
    position(1) = position(1)-0.2*position(3);
    
    ax_6(iUnit) = axes('position', position,'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    for iEpoch = 1:3
        idx = iEpoch;
        currEpoch = epochNames{iEpoch};
        
        % data
  
        data = assemblyTuningPFcorr.(currEpoch).data(currUnit);
        
        
        h1 = scatter(idx-0.15, data, 0.5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
        

        % shuffle
        
        ui = assemblyTuningPFcorr.(currEpoch).ui(currUnit, :);
        ui = ui(:);
        [xData, yData] = myViolin(ui);

        h2 = patch(idx+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
        
        try
            pval = signrank(ui, data);
            signn = significanceSign(pval);
            text(idx, 1, signn, 'fontsize', 2, 'color', 'k', 'HorizontalAlignment', 'center', 'FontName', 'Arial')
            
        catch
            
        end
    end
    
    xlim([0 4])
    
    if iUnit == 1
       legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northeastout', 'box', 'off')
    end
    
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

linkaxes(ax_3, 'y')
linkaxes(ax_4, 'y')
linkaxes(ax_5, 'y')
linkaxes(ax_6, 'y')

    
%% sub-functions


function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end


function [xData, yData] = myViolin(data)

binEdges = linspace(min(data)-0.1*range(data), max(data)+0.1*range(data), 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
if size(count, 1) > size(count, 2)
    count = count';
end

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

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
