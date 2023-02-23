clear
clc


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
nSessions    = numel(currSessions);

assemblyTuningPFcorr_zscore_pooled = cell(nSessions, 1);

epochsXCorr        = cell(nSessions, 1);
epochsXCorr_zscore = cell(nSessions, 1);
activeUnits        = cell(nSessions, 1);


for iSess = 1:numel(currSessions)

sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name;

basePath      = fullfile(parentDir, sessionName);



% spikes and behavior data
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
behavior = fileInfo.behavior;
behavior.time = behavior.time/3600;

startT.pre  = behavior.time(1,1); endT.pre  = behavior.time(2,1); 
startT.run  = behavior.time(2,1); endT.run  = behavior.time(2,2); 
startT.post = behavior.time(2,2); endT.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 


% load learned tunings
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', 'assemblyTunings_time', 'activeUnits');

binCenters            = s.binCenters;
assemblyTunings_time  = s.assemblyTunings_time;
activeUnits{iSess}    = s.activeUnits;


load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTuningPFcorr_zscore');
assemblyTuningPFcorr_zscore_pooled{iSess} = assemblyTuningPFcorr_zscore;



epochNames = {'pre'; 'run'; 'post'};
for iEpoch = 1:numel(epochNames)
    binCenters.(epochNames{iEpoch}) = binCenters.(epochNames{iEpoch})/3600;
end



% truncate the long recording sessions to the first 4 hours of post sleep
idx = find(binCenters.post < endT.post);

binCenters.post = binCenters.post(idx);
assemblyTunings_time.post.data = assemblyTunings_time.post.data(:, :, idx);
assemblyTunings_time.post.ui   = assemblyTunings_time.post.ui(:, :, idx, :);



spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

nUnits = numel(spikes);


% bidirectional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end



%% for assemblyTunings_time calculate the correlation matrix of learned tunings

dataTypes = {'data'; 'ui'};
epochNames = {'pre'; 'run'; 'post'};


for idataType = 1:2
    currDataType = dataTypes{idataType}; 
    
    nInstances = size(assemblyTunings_time.pre.(currDataType), 4);       
    epochsXCorr{iSess}.(currDataType) = nan(3,3, nUnits, nInstances);  


    for iUnit = 1:nUnits

        for inst = 1:nInstances
            for iEpoch = 1:3 
                tuning1 = squeeze(assemblyTunings_time.(epochNames{iEpoch}).(currDataType)(iUnit, :, :, inst));

                for jEpoch = iEpoch:3
                    tuning2 = squeeze(assemblyTunings_time.(epochNames{jEpoch}).(currDataType)(iUnit, :, :, inst));

                    corrMat = corr(tuning1, tuning2);

                    if jEpoch == iEpoch
                        corrMat = corrMat - eye(size(corrMat, 1)).*corrMat;
                    end
                    epochsXCorr{iSess}.(currDataType)(iEpoch, jEpoch, iUnit, inst) = nanmean(corrMat(:));
                    
                end
            end
        end
    end    
end

epochsXCorr_zscore{iSess} = (epochsXCorr{iSess}.data - mean(epochsXCorr{iSess}.ui, 4))./std(epochsXCorr{iSess}.ui, [], 4);

end



%% plot the distribution of cross-correlations for example units


selectedUnits = [33 54 51 60 17 2];


% for the example units (data and the distribution of ui)

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
            h1= scatter(idx-0.1, epochsXCorr{iSess}.data(iEpoch, jEpoch, currUnit), 5, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);

            % ui
            currXcorrData = epochsXCorr{iSess}.ui(iEpoch, jEpoch, currUnit, :);

            h2=line([idx+0.1 idx+0.1], [nanprctile(currXcorrData, 25) nanprctile(currXcorrData, 75)], 'linewidth', 0.25, 'color', [0 0 0 0.5]);
%             h2=scatter(idx+0.2, nanmedian(currXcorrData), 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
            
            line([idx idx+0.2], [nanprctile(currXcorrData, 25) nanprctile(currXcorrData, 25)], 'linewidth', 0.5, 'color', [0 0 0 0.5])
            line([idx idx+0.2], [nanprctile(currXcorrData, 75) nanprctile(currXcorrData, 75)], 'linewidth', 0.5, 'color', [0 0 0 0.5])

        end
    end
    
    if iUnit == 1
        ylabel({'mean';'cross-correlation'}, 'fontsize', fontsize)
    end
    
    if iUnit == 6
       legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')

    end
    
    xlim([0.5 6.5])
%     ylim([-1 1])
%     grid on
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    set(ax(iUnit), 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'linewidth', 1, 'box', 'off', 'xtick', 1:6, 'xticklabel', {'pre'; 'run'; 'post'; 'pre-run'; 'pre-post'; 'run-post'})
    xtickangle(ax(iUnit), 45)
    set(gca, 'ytick', -1:0.5:1)
    
    grid on
    
end
 
linkaxes(ax, 'y')





%% plot the dsitribution of cross-correlation for an example session (panel D)


iSess = 6;

plotwidth  = 130;
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
tt = 4;
for iEpoch = 1:3

    for jEpoch = iEpoch:3

        if jEpoch == iEpoch
            idx = iEpoch;
        else
            idx = tt;
            tt = tt+1;
        end

        currData = epochsXCorr_zscore{iSess}(iEpoch, jEpoch, okUnits);

        scatter(idx+0.1*randn(numel(okUnits), 1), currData, 2, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
                
        scatter(idx+0.4, nanmedian(currData), 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k') 
        
        line([idx+0.4 idx+0.4], [nanprctile(currData, 25) nanprctile(currData, 75)], 'linewidth', 1, 'color', [0 0 0 0.8]);
        
        line([idx+0.3 idx+0.5], [nanprctile(currData, 25) nanprctile(currData, 25)], 'linewidth', 1, 'color', [0 0 0 0.8])
        line([idx+0.3 idx+0.5], [nanprctile(currData, 75) nanprctile(currData, 75)], 'linewidth', 1, 'color', [0 0 0 0.8])

    end
end  

set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:6, 'xticklabel', {'pre'; 'run'; 'post'; 'pre-run'; 'pre-post'; 'run-post'})
xtickangle(gca, 45)
ylabel({'mean'; 'cross-correlation(z)'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';





%% plot the distribution of cross-corrrelation for all sessions

plotwidth  = 230;
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

tt = 4;
for iEpoch = 1:3
        
    for jEpoch = iEpoch:3

        if jEpoch == iEpoch
            idx = iEpoch;
        else
            idx = tt;
            tt = tt+1;
        end

        
        currData = cell(nSessions, 1);
        for iSess = 1:nSessions
            curr_x = idx+(iSess-ceil(nSessions/2))*0.85/nSessions;
            okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);

            
            currData{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, jEpoch, okUnits));
            
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
        
        curr_x = idx+(iSess+1-ceil(nSessions/2))*0.85/nSessions;
        
        scatter(curr_x, med, 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)
        
        % whiskers
        line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
        line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

        % box
        patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)


    end
    
    
end
xlim([0.5 6.6])

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:6, 'xticklabel', {'pre'; 'run'; 'post'; 'pre-run'; 'pre-post'; 'run-post'})
xtickangle(gca, 45)
ylabel({'mean'; 'cross-correlation(z)'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



%% correlation between LT consistency during an epoch and its correlation PFs


plotwidth  = 230;
plotheight = 130;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 3;

leftmargin = 30;  rightmargin = 30;    topmargin = 10;    bottommargin = 30;

gapc = 5;
gapr = 8;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);



nSessions = size(epochsXCorr_zscore, 1);


clear ax
for iEpoch = 1
    
    
    ax(iEpoch) = axes('position',sub_pos{iEpoch},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','off','Layer','top');
    hold on
        
    LTcons  = cell(nSessions, 1); % LT consistency
    PFcorrs = cell(nSessions, 1);
    for iSess = 1:nSessions
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTcons{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(epochNames{iEpoch})(okUnits);
    end

    LTcons  = cell2mat(LTcons);
    PFcorrs = cell2mat(PFcorrs);

%     figure;
    scatter(LTcons, PFcorrs, 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
        
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize)
    
    if iEpoch == 1
        ylabel('PF correlation', 'fontsize', fontsize)
    end
    
    if iEpoch >= 2
        set(gca, 'yticklabel', {})
    end
    
    if iEpoch == 2
        xlabel('learned tuning consistency', 'fontsize', fontsize)
    end
    
    title(epochNames{iEpoch}, 'fontsize', fontsize, 'fontweight', 'normal')
    grid on
end
    
linkaxes(ax, 'y') 



%% plots

selectedUnits = [33 54 51 60 17 2];


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



for iUnit = 1:numel(selectedUnits)
    currUnit = selectedUnits(iUnit); 
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on

    currSpatialTuning = spikes(currUnit).spatialTuning_smoothed.uni;

    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

    end
    
    plot(0.1*currSpatialTuning + startT.run , 1:numel(currSpatialTuning), 'linewidth', 1.5, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(ax(iUnit));
    
    if iUnit == 4
        line(0.2*[0 2] + startT.run, [yl(2)-10 yl(2)-10], 'linewidth', 1.5, 'color', [0 0 0 0.6])
        text(startT.run, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+12, ['unit ' num2str(iUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1])
    
    if ismember(iUnit, [1 4])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters.(currEpoch);
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




%% plot cross-correlation matrix of learned tunings in different time
% windows

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


cnctBinCenters = [binCenters.pre binCenters.run binCenters.post];
tickLabels = 0:2:cnctBinCenters(end);

for iUnit = 1:numel(selectedUnits)
    
    
    currUnit = selectedUnits(iUnit);
    
    
    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time.(currEpoch).data(currUnit, :, :, 1));

    end
    
    tmpConcat = cell2mat(learnedTuning');

    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
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
           h = text(-4, median(binCenters.(currEpoch)), currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           h = text(-2, median(cnctBinCenters), 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           
           text(median(binCenters.(currEpoch)), -4 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           
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
    set(gca, 'box', 'off', 'fontsize', fontsize, 'position', sub_pos{iUnit})
    
end
    

function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end
