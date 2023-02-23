
clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = 1:17;
nSessions    = numel(currSessions);

nShuffles = 1000;

startT = cell(nSessions, 1);
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
% assemblyTuningPFcorr_pscore = cell(nSessions, 1);

epochsXCorr         = cell(nSessions, 1);
epochsXCorr_zscore  = cell(nSessions, 1);

wins                = cell(nSessions, 1);

% KL  diveregenc (the same manner as the Pearson correlation)

% assemblyTuningPFKLdiv        = cell(nSessions, 1);
% assemblyTuningPFKLdiv_pscore = cell(nSessions, 1);
% 
% epochsKLdiv                  = cell(nSessions, 1);
% epochsKLdiv_zscore           = cell(nSessions, 1);

% load(fullfile(parentDir, 'assemblyTunings_plus_PFunitIDshuffles.mat'), 'learnedTuningPFcorr')


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
startT{iSess}.post = behavior.time(2,2); endT{iSess}.post = behavior.time(3,2); 
    

spikes = spikes_pyr;

ntotUnits   = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);



% bidirectional spatial tuning

spatialTunings_merge = zeros(ntotUnits, nPosBins);
for iUnit = 1:ntotUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end

spikes_pooled{iSess} = spikes;



% % % SWR LEARNED TUNINGS


% % learned tunings calculated based on the entire epochs
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), ...
    'activeUnits', ...
    'assemblyTuningPFcorr');%, ...
%     'assemblyTuningPFKLdiv');

activeUnits{iSess} = s.activeUnits;

okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
nUnits = numel(okUnits);

assemblyTuningPFcorr{iSess}  = s.assemblyTuningPFcorr;




% % learned tunings in 15-minutes time windows

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time_Jan2022.mat']), ...
        'binCenters', 'assemblyTunings_time');
    
binCenters{iSess}           = s.binCenters;
assemblyTunings_time{iSess} = s.assemblyTunings_time;

for iEpoch = 1:numel(epochNames)
    binCenters{iSess}.(epochNames{iEpoch}) = binCenters{iSess}.(epochNames{iEpoch})/3600;
end

% remove learned tunnigs in overlapping time windows
% we used 15minutes-long sliding windows with 5minutes step size, so we
% need to keep every third time windows
    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap{iSess}.(epochName) = binCenters{iSess}.(epochName)(1:3:end);
    assemblyTunings_time_nonOverlap{iSess}.(epochName).data = assemblyTunings_time{iSess}.(epochName).data(okUnits, :, 1:3:end);
end



%% calculate a cross-correlation matrix for the 15minues time window learned tunings

clear currData

maxNumWins = 12;
winLength  = 2; % in hours

epochsXCorr{iSess}.data   = nan(3, nUnits, maxNumWins);
epochsXCorr{iSess}.ui     = nan(3, nShuffles, maxNumWins);
epochsXCorr_zscore{iSess} = nan(3, nUnits, maxNumWins); 
   
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    if ismember(iEpoch, 1:2)
        wins{iSess}.(currEpoch) = [startT{iSess}.(currEpoch) endT{iSess}.(currEpoch)];
    elseif iEpoch == 3
        winStarts = startT{iSess}.(currEpoch) + (0:(endT{iSess}.(currEpoch) - startT{iSess}.(currEpoch) - winLength));
        winEnds   = winStarts + winLength;

        wins{iSess}.(currEpoch) = [winStarts' winEnds'];
    end
    
    nWins = size(wins{iSess}.(currEpoch), 1);
    for iWin = 1:nWins

        currWin = wins{iSess}.(currEpoch)(iWin, :);

        idx = find(binCenters_nonOverlap{iSess}.(currEpoch) > currWin(1) & binCenters_nonOverlap{iSess}.(currEpoch) < currWin(2));
    
        currData = assemblyTunings_time_nonOverlap{iSess}.(currEpoch).data(:, :, idx);
    

        for iUnit = 1:nUnits
                       
            % the real tunings
            currTuning_across_time = squeeze(currData(iUnit, :, :));
            
            % calculate the matching scores: pearson correlation and KL
            % divergence
            corrMat = corr(currTuning_across_time, currTuning_across_time, 'type', 'pearson');
                
            offDiagonalIdx = true(size(corrMat));
            for ii = 1:size(corrMat, 1)
                offDiagonalIdx(ii, ii) = false;
            end
            corrMat = corrMat(offDiagonalIdx);
        
            epochsXCorr{iSess}.data(iEpoch, iUnit, iWin) = nanmedian(corrMat(:));
    
        end
            
            
        % % shuffle 
        for inst = 1:nShuffles
    
            % generate shuffle tunings
            nTimeBins = size(currData, 3);
            currTuning_across_time1    = nan(nPosBins, nTimeBins);
            currTuning_across_time2    = nan(nPosBins, nTimeBins);
            for iTimeBin = 1:nTimeBins
                currTuning_across_time1(:, iTimeBin) = squeeze(currData(randperm(nUnits, 1), :, iTimeBin));
                currTuning_across_time2(:, iTimeBin) = squeeze(currData(randperm(nUnits, 1), :, iTimeBin));
            end
    
    
            % calculate the matching scores
            corrMat = corr(currTuning_across_time1, currTuning_across_time2);
            
            offDiagonalIdx = true(size(corrMat));
            for ii = 1:size(corrMat, 1)
                offDiagonalIdx(ii, ii) = false;
            end
            corrMat = corrMat(offDiagonalIdx);
            
    
            epochsXCorr{iSess}.ui(iEpoch, inst, iWin) = nanmedian(corrMat(:));
    
        end
    end

end

epochsXCorr_zscore{iSess} = ((epochsXCorr{iSess}.data) - nanmean(epochsXCorr{iSess}.ui, 2)) ./ nanstd(epochsXCorr{iSess}.ui, [], 2);

end


%% 


clear currData

plotwidth  = 380;
plotheight = 100;
fontsize   = 6;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    currData.(currEpoch) = cell(nSessions, 1);



    % what should be the x-axis coordinate of the data

    if iEpoch == 1 % PRE

        PREwins = nan(nSessions, 2);
        for iSess = 1:nSessions
            PREwins(iSess, :) = wins{iSess}.pre;
        end
        PREdur = diff(PREwins')';
        [~, PREidx] = max(PREdur);
        curr_x = mean(PREwins(PREidx, :));

    elseif iEpoch == 2 % MAZE/RUN

        RUNwins = nan(nSessions, 2);
        for iSess = 1:nSessions
            RUNwins(iSess, :) = wins{iSess}.run;
        end
        RUNdur = diff(RUNwins')';
        [~, RUNidx] = max(RUNdur);
        curr_x = PREwins(PREidx, 2) + 0.5*(RUNdur(RUNidx));

    elseif iEpoch == 3 % POST
        curr_x = PREwins(PREidx, 2) + RUNdur(RUNidx) + (0:maxNumWins-1) + winLength/2;
    end

    for iSess = 1:nSessions
        currData.(currEpoch){iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, :, :)); %

        % plot for individual sessions

        if ismember(iEpoch, 1:2)
            plot(curr_x+0.2, nanmedian(currData.(currEpoch){iSess}, 1), 'color', cl(iSess, :), 'Marker', '.', 'MarkerSize', 5)
        else
            plot(curr_x, nanmedian(currData.(currEpoch){iSess}, 1), 'color', [cl(iSess, :) 0.5])
        end

    end


    % plot pooled data
    pooledData = cell2mat(currData.(currEpoch));
    lq  = prctile(pooledData, 25);
    hq  = prctile(pooledData, 75);
    med = prctile(pooledData, 50);
    
    if iEpoch == 3

        idx = ~isnan(lq);
        lq = lq(idx);
        hq = hq(idx);
        med = med(idx);
        curr_x = curr_x(idx);

        patch([curr_x fliplr(curr_x)], [lq fliplr(hq)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        plot(curr_x, med, 'color', [0 0 0 .8], 'linewidth', 1); 
    else
        
        lq = lq(1);
        hq = hq(1);

        line([curr_x curr_x], [lq hq], 'color', [0 0 0 .8], 'linewidth', 1)
        line([curr_x-0.2 curr_x+0.2], [lq lq], 'color', [0 0 0 0.8], 'linewidth', 1)
        line([curr_x-0.2 curr_x+0.2], [hq hq], 'color', [0 0 0 0.8], 'linewidth', 1)

        plot(curr_x, med, 'color', [0 0 0 .8], 'Marker', '.', 'MarkerSize', 10)

    end

end

xlim([0 17])
xl = xlim;


set(gca,'box', 'off', 'linewidth', 1, 'XTick', 0:2:xl(2), 'XTickLabels', '', 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) % ,  'YTick', [0 25 50 75 100]

xticklabels([])


axis tight
currAx = gca;
currAx.YGrid = 'on';



%% plot the distribution of cross-corrrelation for all sessions

plotwidth  = 170;
plotheight = 120;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on

nSessions = size(epochsXCorr_zscore, 1);
med_pooled = nan(3,1);
p_value    = nan(3,1);
z_value    = nan(3,1);

med= [];

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};

    currData.(currEpoch) = cell(nSessions, 1);
    for iSess = 1:nSessions
        curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        currData.(currEpoch){iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits)); % for plotting the stability
 

        med.(currEpoch)(iSess) = nanmedian(currData.(currEpoch){iSess});
        lq  = nanprctile(currData.(currEpoch){iSess}, 25);
        hq  = nanprctile(currData.(currEpoch){iSess}, 75);
        data_min = min(currData.(currEpoch){iSess});
        data_max = max(currData.(currEpoch){iSess});


        scatter(curr_x, med.(currEpoch)(iSess), 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

        % whiskers
        line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
        line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

        % box
        patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    end
    
    pooledData.(currEpoch) = cell2mat(currData.(currEpoch));
    
    [p_value(tt), ~, stats] = signrank(pooledData.(currEpoch), 0, 'tail', 'right');
    z_value(tt) = stats.zval;


    med_pooled(tt) = nanmedian(pooledData.(currEpoch));
    lq  = nanprctile(pooledData.(currEpoch), 25);
    hq  = nanprctile(pooledData.(currEpoch), 75);
    data_min = min(pooledData.(currEpoch));
    data_max = max(pooledData.(currEpoch));

    curr_x = tt + (iSess+1-ceil(nSessions/2))*0.85/nSessions;

    scatter(curr_x, med_pooled(tt), 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    
    tt = tt + 1;
end 


currData.pre_w_post = cell(nSessions, 1);
for iSess = 1:nSessions
    curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
    okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);

    currData.pre_w_post{iSess} = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits)); % for plotting the stability

    med.pre_w_post(iSess) = nanmedian(currData.pre_w_post{iSess});
    lq  = nanprctile(currData.pre_w_post{iSess}, 25);
    hq  = nanprctile(currData.pre_w_post{iSess}, 75);
    data_min = min(currData.pre_w_post{iSess});
    data_max = max(currData.pre_w_post{iSess});


    scatter(curr_x, med.pre_w_post(iSess), 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
end

pooledData.pre_w_post = cell2mat(currData.pre_w_post);

[p_value(tt), ~, stats] = signrank(pooledData.pre_w_post, 0, 'tail', 'right');
z_value(tt) = stats.zval;


med_pooled(tt) = nanmedian(pooledData.pre_w_post);
lq  = nanprctile(pooledData.pre_w_post, 25);
hq  = nanprctile(pooledData.pre_w_post, 75);
data_min = min(pooledData.pre_w_post);
data_max = max(pooledData.pre_w_post);

curr_x = tt+(iSess+1-ceil(nSessions/2))*0.85/nSessions;

scatter(curr_x, med_pooled(tt), 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

% whiskers
line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

% box
patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)





xlim([0.5 3.5])

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])


xtickangle(gca, 45)
ylabel('LT stability/similarity (z)', 'fontsize', fontsize)

ax = gca;
ax.YGrid = 'on';



%% (added later) violin plot of pooled data

iSess = 2;

plotwidth  = 120;
plotheight = 120;
fontsize   = 6;



f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

med_pooled     = nan(3,1);
p_value = nan(3, 1);
z_value = nan(3, 1);

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
   
    currMeds = pooledData.(currEpoch);
   
    indivSessionMedian = med.(currEpoch);
    
    dots = indivSessionMedian(setdiff(1:nSessions, 2));
    scatter(tt+0.1*randn(numel(dots), 1), dots, 2, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    scatter(tt+0.1*randn(1), indivSessionMedian(2), 2, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);


    [xData, yData] = myViolin(currMeds);
    patch(tt+4*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    
 
    med_pooled(tt) = nanmedian(currMeds);

   
    [p_value(tt), ~, stats] = signrank(indivSessionMedian, 0, 'tail', 'right');
    line([tt-0.5 tt+0.5], [med_pooled(tt) med_pooled(tt)], 'linewidth', 1, 'color', 'k')
    
    
    tt = tt + 1;
end

currMeds = pooledData.pre_w_post;


indivSessionMedian = med.pre_w_post;

dots = indivSessionMedian(setdiff(1:nSessions, 2));
scatter(tt+0.1*randn(numel(dots), 1), dots, 2, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
scatter(tt+0.1*randn(1), indivSessionMedian(2), 2, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
   
    
[xData, yData] = myViolin(currMeds);
patch(tt+4*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)

med_pooled(tt) = nanmedian(currMeds);

[p_value(tt), ~, stats] = signrank(indivSessionMedian, 0, 'tail', 'right');% currMeds


line([tt-0.5 tt+0.5], [med_pooled(tt) med_pooled(tt)], 'linewidth', 1, 'color', 'k')   


xlim([0.5 3.5])

set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])

xtickangle(gca, 45)
ylabel('LT stability/similarity (z)', 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



% comparing between the different epochs/conditions



[pp, ~, stats] = signrank(pooledData.pre, 0, 'tail', 'right');
zz = stats.zval;



%% sub-functions

function prctileScore = calPrctile(variable)


[nUnits, nShuffles] = size(variable.ui);
prctileScore = zeros(nUnits, 1);

for iUnit = 1:nUnits
    
    if isnan(variable.data(iUnit))
        prctileScore(iUnit) = nan;
    else
        prctileScore(iUnit) = numel(find(variable.ui(iUnit, :) <= variable.data(iUnit)))/nShuffles * 100;
    end

end

end


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

gw = gausswindow(6,18);
count = conv(count, gw, 'same');

xData = [count zeros(1, numel(count))];
yData = [binCenters fliplr(binCenters)];

end
