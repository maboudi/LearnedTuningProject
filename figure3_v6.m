

clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
nSessions    = numel(currSessions);

nShuffles = 1000;

startT = cell(nSessions, 1);sq
endT   = cell(nSessions, 1);


activeUnits = cell(nSessions, 1);
epochNames = {'pre'; 'run'; 'post'};


% initialize the variables
spikes_pooled         = cell(nSessions, 1);
binCenters            = cell(nSessions, 1);
binCenters_nonOverlap = cell(nSessions, 1);

assemblyTunings_time            = cell(nSessions, 1);
assemblyTunings_time_nonOverlap = cell(nSessions, 1);



% Pearson correlation; between learned tunnings of different units and between learned tunings
% and the correponding place fields of individual units

assemblyTuningPFcorr        = cell(nSessions, 1);
assemblyTuningPFcorr_zscore = cell(nSessions, 1);

assemblyTuningPFcorr_time        = cell(nSessions, 1);
% assemblyTuningPFcorr_time_zscore = cell(nSessions, 1);

epochsXCorr         = cell(nSessions, 1);
epochsXCorr_zscore  = cell(nSessions, 1);


% KL  diveregenc (the same manner as the pearson correlation)

assemblyTuningPFKLdiv              = cell(nSessions, 1);
assemblyTuningPFKLdiv_zscore       = cell(nSessions, 1);

assemblyTuningPFKLdiv_time         = cell(nSessions, 1);
% assemblyTuningPFKLdiv_time_zscore  = cell(nSessions, 1);

epochsKLdiv                        = cell(nSessions, 1);
epochsKLdiv_zscore                 = cell(nSessions, 1);




%% calculate matching scores for all sessions

for iSess = 1:nSessions

    
sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name
basePath      = fullfile(parentDir, sessionName);


% % % spikes and behavior data
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'fileInfo', 'spikes_pyr')

behavior = fileInfo.behavior;
behavior.time = behavior.time/3600; % to hours

startT{iSess}.pre  = behavior.time(1,1); endT{iSess}.pre  = behavior.time(2,1); 
startT{iSess}.run  = behavior.time(2,1); endT{iSess}.run  = behavior.time(2,2); 
startT{iSess}.post = behavior.time(2,2); endT{iSess}.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 



spikes = spikes_pyr;

nUnits   = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);



% bidirectional spatial tuning

spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end

spikes_pooled{iSess} = spikes;



% % % SWR LEARNED TUNINGS


% % learned tunings calculated based on the entire epochs
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), ...
    'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore', 'asTuning_PF_KLdiv', 'asTuning_PF_KLdiv_zscore');

assemblyTuningPFcorr{iSess} = s.assemblyTuningPFcorr;
assemblyTuningPFcorr_zscore{iSess} = s.assemblyTuningPFcorr_zscore;
assemblyTuningPFKLdiv{iSess} = s.asTuning_PF_KLdiv;
assemblyTuningPFKLdiv_zscore{iSess} = s.asTuning_PF_KLdiv_zscore;



% % learned tunings in 15-minutes time windows
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', ...
    'assemblyTunings_time', 'assemblyTuningPFcorr_time_zscore', ...
    'assemblyTuningPFcorr_time', 'asTuning_PF_KLdiv_time', 'asTuning_PF_KLdiv_time_zscore', ...
    'activeUnits');

binCenters{iSess}                  = s.binCenters;
assemblyTunings_time{iSess}        = s.assemblyTunings_time;
activeUnits{iSess}                 = s.activeUnits;


assemblyTuningPFcorr_time{iSess}         = s.assemblyTuningPFcorr_time;
% assemblyTuningPFcorr_time_zscore{iSess}  = s.assemblyTuningPFcorr_time_zscore;
assemblyTuningPFKLdiv_time{iSess}        = s.asTuning_PF_KLdiv_time;
% assemblyTuningPFKLdiv_time_zscore{iSess} = s.asTuning_PF_KLdiv_time_zscore;



for iEpoch = 1:numel(epochNames)
    binCenters{iSess}.(epochNames{iEpoch}) = binCenters{iSess}.(epochNames{iEpoch})/3600;
end



% okUnits are the units with at least 100 PBEs duing PRE and POST
% okUnits = intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post);



% truncate the long recording sessions to the first 4 hours of post sleep

idx = find(binCenters{iSess}.post < endT{iSess}.post);

binCenters{iSess}.post = binCenters{iSess}.post(idx);
assemblyTunings_time{iSess}.post.data = assemblyTunings_time{iSess}.post.data(:, :, idx);



% remove learned tunnigs in overlapping time windows
% we used 15minutes-long sliding windows with 5minutes step size, so we
% need to keep every third time windows

    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap{iSess}.(epochName) = binCenters{iSess}.(epochName)(1:3:end);
    assemblyTunings_time_nonOverlap{iSess}.(epochName).data = assemblyTunings_time{iSess}.(epochName).data(:, :, 1:3:end);

end



%% calculate a cross-correlation matrix for the 15minues time window learned tunings
% removed the part related to KL divergence

epochsXCorr{iSess}.data   = nan(3, 3, nUnits);
epochsXCorr{iSess}.ui     = nan(3, 3, nShuffles);
epochsXCorr_zscore{iSess} = nan(3, 3, nUnits); 
   

for iEpoch = 1:3    
    iEpoch
    
    currData1 = assemblyTunings_time_nonOverlap{iSess}.(epochNames{iEpoch}).data;

    for jEpoch = iEpoch:3
        currData2 = assemblyTunings_time_nonOverlap{iSess}.(epochNames{jEpoch}).data;

        for iUnit = 1:nUnits
            
            otherUnits = setdiff(1:nUnits, iUnit);
                        
            % % data
            
            % the real tunings
            tuning1 = squeeze(currData1(iUnit, :, :));
            tuning2 = squeeze(currData2(iUnit, :, :));

            
            % calculate the matching scores: pearson correlation and KL
            % divergence
            corrMat = corr(tuning1, tuning2, 'type', 'pearson');
            
            if jEpoch == iEpoch
                
                offDiagonalIdx = true(size(corrMat));
                for ii = 1:size(corrMat, 1)
                    offDiagonalIdx(ii, ii) = false;
                end
                corrMat = corrMat(offDiagonalIdx);
            end
            
            epochsXCorr{iSess}.data(iEpoch, jEpoch, iUnit) = nanmedian(corrMat(:));

        end
        
        
        % % shuffle 
        for inst = 1:nShuffles

            % generate shuffle tunings
            nTimeBins1 = size(currData1, 3);
            tuning1    = zeros(nPosBins, nTimeBins1);
            for iTimeBin = 1:nTimeBins1
                tuning1(:, iTimeBin) = squeeze(currData1(randperm(nUnits, 1), :, iTimeBin));
            end


            nTimeBins2 = size(currData2, 3);
            tuning2 = zeros(nPosBins, nTimeBins2);
            for iTimeBin = 1:nTimeBins2 
                tuning2(:, iTimeBin) = squeeze(currData2(randperm(nUnits, 1), :, iTimeBin));
            end


            % calculate the matching scores
            corrMat = corr(tuning1, tuning2, 'type', 'pearson');

            if jEpoch == iEpoch
                offDiagonalIdx = true(size(corrMat));
                for ii = 1:size(corrMat, 1)
                    offDiagonalIdx(ii, ii) = false;
                end
                corrMat = corrMat(offDiagonalIdx);
            end

            epochsXCorr{iSess}.ui(iEpoch, jEpoch, inst) = nanmedian(corrMat(:));

        end
    end    
end

epochsXCorr_zscore{iSess} = ((epochsXCorr{iSess}.data) - nanmean(epochsXCorr{iSess}.ui, 3)) ./ nanstd(epochsXCorr{iSess}.ui, [], 3);

end


% %% plot the 15-minutes time windows learned tunings and matching scores for all units in a session
% 
% 
% plotLearnedTunings_allUnits




%% plot the learning tunings across time for example units in an example session

iSess = 1;
selectedUnits = [72 117 79 63 116 90]; % Achilles

spikes = spikes_pooled{iSess};

% selectedUnits = [5 15 20 25 35 55]; % rat U


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

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time{iSess}.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch}./repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters{iSess}.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

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
    
    xlim([0 binCenters{iSess}.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    if ismember(iUnit, [1 4])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters{iSess}.(currEpoch);
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


cnctBinCenters = [binCenters{iSess}.pre binCenters{iSess}.run binCenters{iSess}.post];
tickLabels = 0:2:cnctBinCenters(end);

for iUnit = 1:numel(selectedUnits)
    
    
    currUnit = selectedUnits(iUnit);
    
    
    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time{iSess}.(currEpoch).data(currUnit, :, :, 1));

    end
    
    tmpConcat = cell2mat(learnedTuning');

    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters{iSess}.(epochNames{iEpoch}), binCenters{iSess}.(epochNames{jEpoch}), corrMat)
 
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
           h = text(-4, median(binCenters{iSess}.(currEpoch)), currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           h = text(-2, median(cnctBinCenters), 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           
           text(median(binCenters{iSess}.(currEpoch)), -4 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           
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
    
    
    tt = 0;
    for iEpoch = [1 3] % pre and post within-epoch consistency of learned tunings over time
        tt = tt + 1; 
        % data
        data = epochsXCorr{iSess}.data(iEpoch, iEpoch, currUnit);
        
%         [xData, yData] = myViolin(data);
% 
%         patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')

        % ui
        ui = squeeze(epochsXCorr{iSess}.ui(iEpoch, iEpoch, :));
        
        
        try 
        [xData, yData] = myViolin(ui);
        patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
        scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')

        catch
            continue
        end
        
    end
    
    
    tt = tt+1;
    
    % similarity of learned tunings between pre and post
    
    % data
    data = epochsXCorr{iSess}.data(1, 3, currUnit, 1);
        
%     [xData, yData] = myViolin(data);
% 
%     h1 = patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    h1 = scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none'); 
    
    
    % ui
    ui = squeeze(epochsXCorr{iSess}.ui(1, 3, :));
    
    try
        [xData, yData] = myViolin(ui);    
        h2 = patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none');
    catch
        continue
    end
    
    
    
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
yl = ylim;

ylim([yl(1)-0.1 yl(2)+0.1])




%% plot the distribution of cross-correlation for an example session (panel D)

iSess = 7;

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
    
    
%     currMeds = nan(numel(okUnits), 1);
%     for iUnit = 1:numel(okUnits)
%         currData = squeeze(epochsXCorr{iSess}.data{iEpoch, iEpoch, okUnits(iUnit), 1});  
%         currMeds(iUnit) = nanmedian(currData);
%     end

    currMeds = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits, 1));

%     signIdx  = squeeze(epochsXCorr_ui_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4 & squeeze(epochsXCorr_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4;
%     scatter(tt+0.05*randn(numel(find(~signIdx)), 1), currMeds(~signIdx), 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
%     scatter(tt+0.05*randn(numel(find(signIdx)), 1), currMeds(signIdx), 1, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);

    scatter(tt+0.05*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    
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
thresh_z = 2;

clear ax med LTstability_pooled PFcorrs_pooled


for iEpoch = [1 2 3]
    
    hold on
       
    
    LTstability  = cell(nSessions, 1); % LT consistency
    
    PFcorrs_zscore = cell(nSessions, 1);
    KLdivs_zscore  = cell(nSessions, 1);
    
    PFcorrs = cell(nSessions, 1);
    KLdivs  = cell(nSessions, 1);
    
    
    for iSess = setdiff(1:nSessions, [10 11])
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        
        
        
        PFcorrs{iSess} = assemblyTuningPFcorr{iSess}.(epochNames{iEpoch}).data(okUnits);
        KLdivs{iSess} = assemblyTuningPFKLdiv{iSess}.(epochNames{iEpoch}).data(okUnits);
        
        
        PFcorrs_zscore{iSess} = assemblyTuningPFcorr_zscore{iSess}.(epochNames{iEpoch})(okUnits);
        KLdivs_zscore{iSess}  = assemblyTuningPFKLdiv_zscore{iSess}.(epochNames{iEpoch})(okUnits);

%         PFcorrs{iSess}     = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits));
    end

    
    LTstability_pooled.(epochNames{iEpoch}) = cell2mat(LTstability);
     
    KLdivs_zscore_pooled.(epochNames{iEpoch})  = cell2mat(KLdivs_zscore);
    PFcorrs_zscore_pooled.(epochNames{iEpoch}) = cell2mat(PFcorrs_zscore);
    
    KLdivs_pooled.(epochNames{iEpoch})  = cell2mat(KLdivs);
    PFcorrs_pooled.(epochNames{iEpoch}) = cell2mat(PFcorrs);

end



tt = 1;
for iEpoch = 1:3
    
    currParam = KLdivs_zscore_pooled;
    
    stableOnes    = LTstability_pooled.(epochNames{iEpoch}) > thresh_z;
    nonStableOnes = LTstability_pooled.(epochNames{iEpoch}) <= thresh_z;
    
    
    pooledData = [currParam.pre; currParam.run; currParam.post];
    
    bottomLim = prctile(pooledData, 0);
    topLim    = prctile(pooledData, 100);
    
    
    data2Plot = currParam.(epochNames{iEpoch})(stableOnes);
%     data2Plot = data2Plot(data2Plot > bottomLim & data2Plot < topLim);
    [xData, yData] = myViolin_oneSided(data2Plot, 0.1);
    h1 = patch(tt+8*xData, yData, [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    
    try
        data2Plot = currParam.(epochNames{iEpoch})(nonStableOnes);
%         data2Plot = data2Plot(data2Plot > bottomLim & data2Plot < topLim);
        [xData, yData] = myViolin_oneSided(data2Plot, 0.1); 
        h2 = patch(tt-8*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    catch
    end
    
    med.(epochNames{iEpoch}) = [nanmedian(currParam.(epochNames{iEpoch})(nonStableOnes)) nanmedian(currParam.(epochNames{iEpoch})(stableOnes))];
    
    tt = tt+1;
end


ax = gca;
ax.YGrid = 'on';

xlim([0 4])


% legend([h1, h2], 'stable LTs', 'nonstable LTs', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')

% xlabel('epoch', 'fontsize', fontsize)
ylabel('PF-LT Pearson correlation coeff. (z)', 'fontsize', fontsize)
% ylabel('PF-LT KL divergence (z)', 'fontsize', fontsize)


yl = ylim;

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'run'; 'post'}, 'ytick', [yl(1):1:yl(2)], 'TickDir', 'out','TickLength',[0.01, 0.01] ) 





%% a small analysis on correlation between LT-PF corrs and LT-PF KLdivs


KLdivs_all = [KLdivs_zscore_pooled.pre  KLdivs_zscore_pooled.post]; 
PFcorrs_zscore_all = [PFcorrs_zscore_pooled.pre  PFcorrs_zscore_pooled.post];

idx =  PFcorrs_zscore_all > prctile(PFcorrs_zscore_all, 0) & PFcorrs_zscore_all < prctile(PFcorrs_zscore_all, 100); % KLdivs_all > prctile(KLdivs_all, 1) & KLdivs_all < prctile(KLdivs_all, 99) &

KLdivs_all = KLdivs_all(idx);
PFcorrs_zscore_all = PFcorrs_zscore_all(idx);



plotwidth  = 150;
plotheight = 150;
fontsize   = 9;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


scatter(KLdivs_all, PFcorrs_zscore_all, 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4)

xlabel('LT-PF KL divergence', 'fontsize', fontsize)
ylabel('LT-PF Pearson correlation', 'fontsize', fontsize)

% title('correlation between PF matching scores', 'fontsize', 6, 'fontweight', 'normal')

[coeff, pval] = corr(KLdivs_all, PFcorrs_zscore_all, 'type', 'pearson');

xl = xlim;
yl = ylim;

text(median(xl)-0.2*range(xl), yl(1)+ 0.2*range(yl), {sprintf('Pearson corr coeff. = %.2f', coeff); sprintf('pval = %.2f', pval)}, 'fontsize', 6, 'color', [0.4 0.4 0.4])

set(gca, 'box', 'off', 'linewidth', 1, 'xtick', [xl(1):1:xl(2)], 'ytick', [yl(1):1:yl(2)], 'TickDir', 'out','TickLength',[0.01, 0.01] )

grid on



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
        
    LTstability     = cell(nSessions, 1); % LT consistency
    PFcorrs_zscore         = cell(nSessions, 1);
    PFcorrs_zscored = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_pooled{iSess}.(currEpoch)(okUnits);
        PFcorrs_zscored{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(currEpoch)(okUnits);
    end
    
    LTstability_pooled.(currEpoch)     = cell2mat(LTstability);
    PFcorrs_pooled.(currEpoch)         = cell2mat(PFcorrs);
    PFcorrs_zscored_pooled.(currEpoch) = cell2mat(PFcorrs_zscored);
    
    
    consisUnits = find(abs(LTstability_pooled.(currEpoch)) > 2);
    
    
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

binEdges = linspace(min(data)-eps, max(data)+eps, 50);

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


function [xData, yData] = myViolin_oneSided(data, binSize)

% data = data(data > prctile(data, 5) & data < prctile(data, 97.5));

binEdges = (min(data)-eps):binSize:(max(data)+eps);

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
