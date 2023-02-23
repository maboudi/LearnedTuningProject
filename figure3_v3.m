clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


% currSessions = [1:5 7:9 10 13:15 6];
currSessions = 13;
nSessions    = numel(currSessions);


startT = cell(nSessions, 1);
endT   = cell(nSessions, 1);


spikes_pooled     = cell(nSessions, 1);
binCenters_pooled = cell(nSessions, 1);


assemblyTunings_time_pooled        = cell(nSessions, 1);
assemblyTuningPFcorr_zscore_pooled = cell(nSessions, 1);
assemblyTuningPFcorr_pooled        = cell(nSessions, 1);

epochsXCorr         = cell(nSessions, 1);
epochsXCorr_zscore  = cell(nSessions, 1);
epochsXCorr_pval    = cell(nSessions, 1);
epochsXCorr_ui_pval = cell(nSessions, 1);


activeUnits = cell(nSessions, 1);

epochNames = {'pre'; 'run'; 'post'};


for iSess = 1:numel(currSessions)

sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name;

basePath      = fullfile(parentDir, sessionName);



% spikes and behavior data
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
behavior = fileInfo.behavior;
behavior.time = behavior.time/3600;

startT{iSess}.pre  = behavior.time(1,1); endT{iSess}.pre  = behavior.time(2,1); 
startT{iSess}.run  = behavior.time(2,1); endT{iSess}.run  = behavior.time(2,2); 
startT{iSess}.post = behavior.time(2,2); endT{iSess}.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

nUnits = numel(spikes);


% bidirectional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end

spikes_pooled{iSess} = spikes;




% load learned tunings
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', 'assemblyTunings_time', 'activeUnits');

binCenters            = s.binCenters;
assemblyTunings_time  = s.assemblyTunings_time;
activeUnits{iSess}    = s.activeUnits;



for iEpoch = 1:numel(epochNames)
    binCenters.(epochNames{iEpoch}) = binCenters.(epochNames{iEpoch})/3600;
end

binCenters_pooled{iSess}           = binCenters;
assemblyTunings_time_pooled{iSess} = assemblyTunings_time;


% load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore');
% 
% 
% for iEpoch = 1:numel(epochNames)
%     
%     currEpoch = epochNames{iEpoch};
%     
%     assemblyTuningPFcorr_pooled{iSess}.(currEpoch)        = assemblyTuningPFcorr.(currEpoch).data;
%     assemblyTuningPFcorr_zscore_pooled{iSess}.(currEpoch) = assemblyTuningPFcorr_zscore.(currEpoch);
%     
%     binCenters.(currEpoch) = binCenters.(currEpoch)/3600;
% 
% end


% remove the time windows with overlap
    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap.(epochName) = binCenters.(epochName)(1:3:end);
    
    assemblyTunings_time_nonOverlap.(epochName).data = assemblyTunings_time.(epochName).data(:, :, 1:3:end);
    assemblyTunings_time_nonOverlap.(epochName).ui  = assemblyTunings_time.(epochName).ui(:, :, 1:3:end, :); 
end


% truncate the long recording sessions to the first 4 hours of post sleep
idx = find(binCenters.post < endT{iSess}.post);

binCenters.post = binCenters.post(idx);
assemblyTunings_time.post.data = assemblyTunings_time.post.data(:, :, idx);
assemblyTunings_time.post.ui   = assemblyTunings_time.post.ui(:, :, idx, :);


idx = find(binCenters_nonOverlap.post < endT{iSess}.post);
binCenters_nonOverlap.post = binCenters_nonOverlap.post(idx);
assemblyTunings_time_nonOverlap.post.data = assemblyTunings_time_nonOverlap.post.data(:, :, idx);
assemblyTunings_time_nonOverlap.post.ui   = assemblyTunings_time_nonOverlap.post.ui(:, :, idx, :);



% calculate the PF correlations for data and ui

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
                currLT = allLTs(:, :, iTimeBin, ins);
                corrMat = corr(currLT', spatialTunings_merge', 'type', 'pearson');
                assemblyTuningPFcorr_time.(currEpoch).(currDataType)(:, iTimeBin, ins) = diag(corrMat);
            end            
        end

    end
    
end



%% for assemblyTunings_time calculate the correlation matrix of learned tunings

dataTypes = {'data'; 'ui'};
epochNames = {'pre'; 'run'; 'post'};

epochsXCorr{iSess}.data    = cell(3,3, nUnits);
epochsXCorr{iSess}.ui      = cell(3,3, nUnits, 100);
epochsXCorr_ui_pval{iSess} = nan(3,3, nUnits);
epochsXCorr_pval{iSess}    = nan(3,3, nUnits);

nInstances = struct('data', 1, 'ui', 100);

for iUnit = 1:nUnits
    
        
    for iEpoch = 1:3

        for jEpoch = iEpoch:3

            for idataType = 1:2
                currDataType = dataTypes{idataType}; 

                for inst = 1:nInstances.(currDataType)

                    tuning1 = squeeze(assemblyTunings_time_nonOverlap.(epochNames{iEpoch}).(currDataType)(iUnit, :, :, inst));
                    tuning2 = squeeze(assemblyTunings_time_nonOverlap.(epochNames{jEpoch}).(currDataType)(iUnit, :, :, inst));

                    corrMat = corr(tuning1, tuning2);

                    if jEpoch == iEpoch

                        allOffDiag = [];
                        for ii = 1:size(corrMat, 1)-1
                           allOffDiag = [allOffDiag; diag(corrMat, ii)]; 
                        end

                        corrMat = allOffDiag;
                    end

                    epochsXCorr{iSess}.(currDataType){iEpoch, jEpoch, iUnit, inst} = corrMat(:);

                end
            end
                

            % calculate the significance level for the calculated
            % correlations

            data = squeeze(epochsXCorr{iSess}.data{iEpoch, jEpoch, iUnit, 1});
            ui   = cell2mat(squeeze(epochsXCorr{iSess}.ui(iEpoch, jEpoch, iUnit, :)));

            try
                epochsXCorr_pval{iSess}(iEpoch, jEpoch, iUnit) = signrank(data);
            catch
                
            end
            
            
            try
                if nanmedian(data) > 0
                    epochsXCorr_ui_pval{iSess}(iEpoch, jEpoch, iUnit) = ranksum(data, ui, 'tail', 'right');

                else
                    epochsXCorr_ui_pval{iSess}(iEpoch, jEpoch, iUnit) = ranksum(data, ui, 'tail', 'left');

                end
                
            catch
                
            end

        end
    end    
 
end


end





%% plot the learning tunings across time for example units in an example session


iSess = 1;


spikes = spikes_pooled{iSess};

selectedUnits = [72 117 79 63 116 90]; %42


plotwidth  = 360;
plotheight = 235;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 3;
noc = 2;

leftmargin = 30;  rightmargin = 30;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 8;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);

nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

for iUnit = 1:numel(selectedUnits)
    currUnit = selectedUnits(iUnit); 
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on

    currSpatialTuning = spikes(currUnit).spatialTuning_smoothed.uni;

    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time_pooled{iSess}.(currEpoch).ui(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters_pooled{iSess}.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

    end
    
    plot(0.1*currSpatialTuning + startT{iSess}.run , 1:nPosBins, 'linewidth', 1.5, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(ax(iUnit));

    
    if iUnit == 4
        line(0.2*[0 2] + startT{iSess}.run, [yl(2)-10 yl(2)-10], 'linewidth', 1.5, 'color', [0 0 0 0.6])
        text(startT{iSess}.run, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+12, ['unit ' num2str(iUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters_pooled{iSess}.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    if ismember(iUnit, [1 4])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters_pooled{iSess}.(currEpoch);
           text(median(currBinCenters), yl(2)+20, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if iUnit == 2
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, [3 6])
        xlabel('time (hour)')
    end
    
    if iUnit >= 4
        set(ax(iUnit), 'YTickLabel', {})
    end
    
    if ~ismember(iUnit, [3 6])
       set(ax(iUnit), 'XTickLabel', {}) 
    end
    
    hold off
    
end




%% plot cross-correlation matrix of learned tunings for the entire session


plotwidth  = 375;
plotheight = 175;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


cnctBinCenters = [binCenters_pooled{iSess}.pre binCenters_pooled{iSess}.run binCenters_pooled{iSess}.post];
tickLabels = 0:2:cnctBinCenters(end);

for iUnit = 1:numel(selectedUnits)
    
    
    currUnit = selectedUnits(iUnit);
    
    
    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time_pooled{iSess}.(currEpoch).data(currUnit, :, :, 1));

    end
    
    tmpConcat = cell2mat(learnedTuning');

    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters_pooled{iSess}.(epochNames{iEpoch}), binCenters_pooled{iSess}.(epochNames{jEpoch}), corrMat)
 
        end
    end
    caxis([-1 1])
    colormap('jet')
    

    ylim([0 cnctBinCenters(end)])
    xlim([0 cnctBinCenters(end)])
    
    yl = ylim;
    text(0, yl(2)+0.5, sprintf('unit %d', iUnit), 'fontsize', fontsize, 'HorizontalAlignment', 'left')
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    if iUnit == 1
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           h = text(-4, median(binCenters_pooled{iSess}.(currEpoch)), currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           h = text(-2, median(cnctBinCenters), 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           
           text(median(binCenters_pooled{iSess}.(currEpoch)), -4 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           
          text(median(cnctBinCenters), -2 , 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 

       end
    end
    
    
    xticks(tickLabels)
    yticks(tickLabels)
    axis square
    
    if iUnit == 6
        h = colorbar;
        h.Label.String = 'correlation';
    end
    set(gca, 'box', 'off', 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'TickDir', 'out','TickLength',[0.01, 0.01])
   
    
end



%% plot the distribution of mean cross-correlations for the example units


% for the example units (distribution of correlation for real data and ui)

plotwidth  = 330;
plotheight = 150;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


for iUnit = 1:numel(selectedUnits)
    
    currUnit = selectedUnits(iUnit);
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
    
    
    tt = 1;
    for iEpoch = [1 3] % pre and post within-epoch consistency of learned tunings over time
            
        % data
        data = squeeze(epochXCorr{iSess}.data{iEpoch, jEpoch, currUnit, 1});
        
        [xData, yData] = myViolin(data);

        patch(idx-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')

        % ui
        ui = cell2mat(squeeze(epochsXCorr{iSess}.ui(iEpoch, iEpoch, currUnit, :)));
        
        [xData, yData] = myViolin(ui);
        patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
              
        tt = tt + 1;
    end
    
    
    
    % similarity of learned tunings between pre and post
    
    % data
    data = squeeze(epochXCorr{iSess}.data{1, 3, currUnit, 1});
        
    [xData, yData] = myViolin(data);

    patch(idx-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    scatter(idx-0.15, nanmedian(data), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')    
    
    
    
    % ui
    ui = cell2mat(squeeze(epochsXCorr{iSess}.ui{1, 3, currUnit, :}));

    [xData, yData] = myViolin(ui);    
    patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    
    scatter(idx+0.15, nanmedian(ui), 0.5, 'filled', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'none')
    if iUnit == 1
        ylabel({'learned tuning'; 'stability'}, 'fontsize', fontsize)
    end
    
    if iUnit == 6
       legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')
    end
    
    xlim([0.5 3.5])
%     ylim([-1 1])
%     grid on
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    set(ax(iUnit), 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'linewidth', 1, 'box', 'off', 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax(iUnit), 45)
    set(gca, 'ytick', -1:0.5:1)
    
    aa = gca;
    aa.YGrid = 'on';
    
end

linkaxes(ax, 'y')




%% plot the distribution of cross-correlation for an example session (panel D)

iSess = 1;

plotwidth  = 70;
plotheight = 90;
fontsize   = 6;


okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

tt = 1;
for iEpoch = [1 3]
    
    
    currMeds = nan(numel(okUnits), 1);
    for iUnit = 1:numel(okUnits)
        currData = squeeze(epochsXCorr{iSess}.data{iEpoch, iEpoch, okUnits(iUnit), 1});  
        currMeds(iUnit) = nanmedian(currData);
    end
        
    signIdx  = squeeze(epochsXCorr_ui_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4 & squeeze(epochsXCorr_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4;
    
    
    scatter(tt+0.05*randn(numel(find(~signIdx)), 1), currMeds(~signIdx), 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    scatter(tt+0.05*randn(numel(find(signIdx)), 1), currMeds(signIdx), 1, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);

    
    [xData, yData] = myViolin(currMeds);
    patch(tt+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    
    tt = tt + 1;
end

xlim([0.5 2.5])

set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
xtickangle(gca, 45)
ylabel({'median cross-correlation'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



%% correlation between .data and _zscore

iEpoch = 1;
jEpoch = 2;

iSess = 6;
okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


figure; scatter(epochsXCorr{iSess}.data(iEpoch, jEpoch, okUnits, 1), epochsXCorr_zscore{iSess}(iEpoch, jEpoch, okUnits, 1), 5, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')


grid on


%% plot the distribution of cross-corrrelation for all sessions

plotwidth  = 90;
plotheight = 90;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

nSessions = size(epochsXCorr_zscore, 1);

tt = 1;
for iEpoch = [1 3]
        
    currData = cell(nSessions, 1);
    for iSess = 1:nSessions
        curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


        currData{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));

        med = nanmedian(currData{iSess});
        lq  = nanprctile(currData{iSess}, 25);
        hq  = nanprctile(currData{iSess}, 75);
        data_min = min(currData{iSess});
        data_max = max(currData{iSess});


        scatter(curr_x, med, 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

        % whiskers
        line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
        line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

        % box
        patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    end
    
    pooledData = cell2mat(currData);

    med = nanmedian(pooledData);
    lq  = nanprctile(pooledData, 25);
    hq  = nanprctile(pooledData, 75);
    data_min = min(pooledData);
    data_max = max(pooledData);

    curr_x = tt+(iSess+1-ceil(nSessions/2))*0.85/nSessions;

    scatter(curr_x, med, 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    
    tt = tt + 1;
end

xlim([0.5 2.5])

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre';'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
xtickangle(gca, 45)
ylabel({'mean'; 'cross-correlation(z)'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



%% correlation between LT consistency during an epoch and its correlation with PFs


plotwidth  = 90;
plotheight = 130;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nSessions = size(epochsXCorr_zscore, 1);


clear ax med LTcons_pooled PFcorrs_pooled

tt = 1;
for iEpoch = 1
    
    hold on
        
    LTcons  = cell(nSessions, 1); % LT consistencybinCenters_pooled
    PFcorrs = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTcons{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(epochNames{iEpoch})(okUnits);
        
    end
    
    LTcons_pooled.(epochNames{iEpoch}) = cell2mat(LTcons);
    PFcorrs_pooled.(epochNames{iEpoch}) = cell2mat(PFcorrs);
    
    [xData, yData] = myViolin_oneSided(PFcorrs_pooled.(epochNames{iEpoch})(LTcons_pooled.(epochNames{iEpoch}) > 2));
    patch(tt+3*xData, yData, [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    med.(epochNames{iEpoch}) = [nanmedian(PFcorrs_pooled.(epochNames{iEpoch})(LTcons_pooled.(epochNames{iEpoch}) < 2)) nanmedian(PFcorrs_pooled.(epochNames{iEpoch})(LTcons_pooled.(epochNames{iEpoch}) > 2))];

    
    [xData, yData] = myViolin_oneSided(PFcorrs_pooled.(epochNames{iEpoch})(LTcons_pooled.(epochNames{iEpoch}) < 0));
    patch(tt-3*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)

    
    tt = tt+1;
    

end

xlim([0.5 2.5])


xlabel('learned tuning consistency', 'fontsize', fontsize)
ylabel('PF correlation', 'fontsize', fontsize)

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', [1 2], 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 






%% correlation between LT consistency during an epoch and its correlation with PFs


plotwidth  = 200;
plotheight = 250;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



nor = 2;
noc = 1;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 0;
gapr = 5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


nSessions = size(epochsXCorr_zscore, 1);


clear ax 

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    
    ax(tt) = axes('position',sub_pos{tt},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
        
    LTcons          = cell(nSessions, 1); % LT consistency
    PFcorrs         = cell(nSessions, 1);
    PFcorrs_zscored = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTcons{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_pooled{iSess}.(currEpoch)(okUnits);
        PFcorrs_zscored{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(currEpoch)(okUnits);
    end
    
    LTcons_pooled.(currEpoch)          = cell2mat(LTcons);
    PFcorrs_pooled.(currEpoch)         = cell2mat(PFcorrs);
    PFcorrs_zscored_pooled.(currEpoch) = cell2mat(PFcorrs_zscored);
    
    
    consisUnits = find(abs(LTcons_pooled.(currEpoch)) > 2);
    
    
%     allCorr = PFcorrs_pooled.(currEpoch)(consisUnits);
    signIdx = abs(PFcorrs_zscored_pooled.(currEpoch)(consisUnits)) > 1 & PFcorrs_zscored_pooled.(currEpoch)(consisUnits).*PFcorrs_pooled.(currEpoch)(consisUnits) > 0; 
    
    
    binEdges   = -1:0.05:1;
    binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;
    
    allCounts    = histc(PFcorrs_pooled.(currEpoch)(consisUnits), binEdges); allCounts(end) = [];
    nonSigCounts = histc(PFcorrs_pooled.(currEpoch)(consisUnits(~signIdx)), binEdges); nonSigCounts(end) = [];
    

    bar(binCenters, allCounts, 'FaceColor', 'k', 'FaceAlpha', 0.8, 'EdgeColor', 'k')
    bar(binCenters, nonSigCounts, 'FaceColor', 'w', 'EdgeColor', 'k')
    
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
    
    tt = tt + 1;

end

% 
% xlim([0.5 2.5])
% 
% 
% xlabel('learned tuning consistency', 'fontsize', fontsize)
% ylabel('PF correlation', 'fontsize', fontsize)




    
%%
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
count = count';

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

xData = [count -fliplr(count)];
yData = [binCenters fliplr(binCenters)];

end


function [xData, yData] = myViolin_oneSided(data)

binEdges = linspace(min(data)-0.1*range(data), max(data)+0.1*range(data), 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
count = count';

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

xData = [count zeros(1, numel(count))];
yData = [binCenters fliplr(binCenters)];

end
