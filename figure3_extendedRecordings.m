


clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = 7:9;
nSessions    = numel(currSessions);


startT = cell(nSessions, 1);
endT   = cell(nSessions, 1);


activeUnits = cell(nSessions, 1);
epochNames = {'pre'; 'run'; 'post'};


% initialize the variables
spikes_pooled         = cell(nSessions, 1);
binCenters            = cell(nSessions, 1);
binCenters_nonOverlap = cell(nSessions, 1);

assemblyTunings_time = cell(nSessions, 1);


% Pearson correlation; between learned tunnings of different units and between learned tunings
% and the correponding place fields of individual units

assemblyTuningPFcorr                   = cell(nSessions, 1);
assemblyTuningPFcorr_pscore            = cell(nSessions, 1);
assemblyTuningPFcorr_time              = cell(nSessions, 1);
assemblyTuningPFcorr_time_shuffle      = cell(nSessions, 1);
assemblyTuningPFcorr_time_prctileScore = cell(nSessions, 1);
    

% KL  diveregenc (the same manner as the Pearson correlation)

assemblyTuningPFKLdiv                   = cell(nSessions, 1);
assemblyTuningPFKLdiv_pscore            = cell(nSessions, 1);
assemblyTuningPFKLdiv_time              = cell(nSessions, 1);
assemblyTuningPFKLdiv_time_shuffle      = cell(nSessions, 1);
assemblyTuningPFKLdiv_time_prctileScore = cell(nSessions, 1);



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
startT{iSess}.post = behavior.time(3,1); endT{iSess}.post = behavior.time(3,2); % in long recordings only the first 4 hours of post sleep is used in this figures 


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
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), ...
    'activeUnits', ...
    'assemblyTuningPFcorr', ...
    'assemblyTuningPFKLdiv');

activeUnits{iSess} = s.activeUnits;
okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);



% % learned tunings in 15-minutes time windows

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time_Jan2022.mat']), ...
        'binCenters', ...
        'assemblyTunings_time', ...
        'assemblyTuningPFcorr_time', ...
        'assemblyTuningPFKLdiv_time');
    
binCenters{iSess}           = s.binCenters;
assemblyTunings_time{iSess} = s.assemblyTunings_time;

for iEpoch = 1:numel(epochNames)
    binCenters{iSess}.(epochNames{iEpoch}) = binCenters{iSess}.(epochNames{iEpoch})/3600;
end

assemblyTuningPFcorr_time{iSess}  = s.assemblyTuningPFcorr_time;
assemblyTuningPFKLdiv_time{iSess} = s.assemblyTuningPFKLdiv_time; 




% truncate the long recording sessions to the first 4 hours of post sleep

idx = find(binCenters{iSess}.post <= endT{iSess}.post);

binCenters{iSess}.post = binCenters{iSess}.post(idx);
assemblyTunings_time{iSess}.post = assemblyTunings_time{iSess}.post.data(:, :, idx);

assemblyTuningPFcorr_time{iSess}.post  = assemblyTuningPFcorr_time{iSess}.post(:, idx);
assemblyTuningPFKLdiv_time{iSess}.post = assemblyTuningPFKLdiv_time{iSess}.post(:, idx);




% calculate surrogate LT - PF matching scores

nShuffles = numel(okUnits) - 1;
for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};
   
    currLTs = s.assemblyTunings_time.(currEpoch).data;
    currSpatialTunings = spatialTunings_merge;
    
    nTimeBins = size(currLTs, 3);
    
    shuffledPFcorrs = nan(nUnits, nShuffles, nTimeBins);
    shuffledKLdivs  = nan(nUnits, nShuffles, nTimeBins); 
    
    PFcorrs_prctileScore = nan(nUnits, nTimeBins);
    KLdiv_prctileScore   = nan(nUnits, nTimeBins);
    
    for iTimeBin = 1:nTimeBins
       
        currTimeBinLTs  = currLTs(:, :, iTimeBin);
        
        % In shuffle matching score keep the spatial tunings the same
        % as in real data

        allCorrs  = corr(currTimeBinLTs', spatialTunings_merge', 'type', 'Pearson');
        allKLdivs = calKLDivergence(currTimeBinLTs', spatialTunings_merge');

        for iUnit = 1:numel(okUnits)
            currUnit = okUnits(iUnit);
            
            shuffledPFcorrs(currUnit, :, iTimeBin) = allCorrs(setdiff(okUnits, currUnit), currUnit);
            shuffledKLdivs(currUnit, :, iTimeBin)  = allKLdivs(setdiff(okUnits, currUnit), currUnit); 
            
            
            currData    = assemblyTuningPFcorr_time{iSess}.(currEpoch)(currUnit, iTimeBin);
            currShuffle =  shuffledPFcorrs(currUnit, :, iTimeBin);
            if ~isnan(currData)
                PFcorrs_prctileScore(currUnit, iTimeBin) = numel(find(currData > currShuffle))/nShuffles * 100;
            end
            
            
            currData    = assemblyTuningPFKLdiv_time{iSess}.(currEpoch)(currUnit, iTimeBin);
            currShuffle = shuffledKLdivs(currUnit, :, iTimeBin);
            if ~isnan(currData)
                KLdiv_prctileScore(currUnit, iTimeBin)   = numel(find(currData < currShuffle))/nShuffles * 100;  
            end
        
        end
        
    end
    
    assemblyTuningPFcorr_time_shuffle{iSess}.(currEpoch)  = shuffledPFcorrs;
    assemblyTuningPFKLdiv_time_shuffle{iSess}.(currEpoch) = shuffledKLdivs;  
    
    assemblyTuningPFcorr_time_prctileScore{iSess}.(currEpoch)  = PFcorrs_prctileScore;
    assemblyTuningPFKLdiv_time_prctileScore{iSess}.(currEpoch) = KLdiv_prctileScore;
    
end

end



%% plot the LT-PF correlation and KL div across time


plotwidth  = 380;
plotheight = 265;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 2;
noc = 1;

leftmargin = 50;  rightmargin = 30;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 15;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);

nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);


unifiedBinsCenters = cell(3,1);


currStart = 0;
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};
    
    allNBins = nan(nSessions, 1);
    for iSess = 1:nSessions    
        allNBins(iSess) = numel(binCenters{iSess}.(currEpoch));
    end
    [maxNBins, sessIdx] = max(allNBins);
    
    unifiedBinsCenters{iEpoch} = currStart + (1:maxNBins)*0.0833;
    currStart = unifiedBinsCenters{iEpoch}(end);
    
    
    % pooling the PF corr across sessions
    
    PFcorrelation.(currEpoch) = cell(nSessions, 1);
    KLdiv.(currEpoch)         = cell(nSessions, 1);
    for iSess = 1:nSessions
        
        nUnits = size(assemblyTuningPFcorr_time{iSess}.(currEpoch), 1);
        
        PFcorrelation.(currEpoch){iSess} = nan(nUnits, maxNBins);
        KLdiv.(currEpoch){iSess}         = nan(nUnits, maxNBins);
        
        switch iEpoch
            case 1
                PFcorrelation.(currEpoch){iSess}(:, end-allNBins(iSess)+1:maxNBins) = assemblyTuningPFcorr_time_prctileScore{iSess}.(currEpoch);
                KLdiv.(currEpoch){iSess}(:, end-allNBins(iSess)+1:maxNBins)         = assemblyTuningPFKLdiv_time_prctileScore{iSess}.(currEpoch);
            case 2
                PFcorrelation.(currEpoch){iSess}(:, 1:allNBins(iSess))              = assemblyTuningPFcorr_time_prctileScore{iSess}.(currEpoch);
                KLdiv.(currEpoch){iSess}(:, 1:allNBins(iSess))                      = assemblyTuningPFKLdiv_time_prctileScore{iSess}.(currEpoch);
            case 3
                PFcorrelation.(currEpoch){iSess}(:, 1:allNBins(iSess))              = assemblyTuningPFcorr_time_prctileScore{iSess}.(currEpoch);
                KLdiv.(currEpoch){iSess}(:, 1:allNBins(iSess))                      = assemblyTuningPFKLdiv_time_prctileScore{iSess}.(currEpoch);
        end
    end

end


% 
% for iSess = 1:nSessions
%     
%     
%     PFcorrelation_refDist = PFcorrelation.run{iSess};
% %     temp = nanmean(PFcorrelation_refDist, 1);
% %     lq = nanprctile(temp(:), 25);
% %     hq = nanprctile(temp(:), 75);
% %     med = nanmedian(temp(:));
% %     
% %     idx = temp > (med+1.5*(hq-lq)) | temp < (med-1.5*(hq-lq));
% %     PFcorrelation_refDist(:, idx) = nan;
%     
% 
%     KLdiv_refDist = KLdiv.run{iSess};
% %     temp = nanmean(KLdiv_refDist, 1);
% %     lq = nanprctile(temp(:), 25);
% %     hq = nanprctile(temp(:), 75);
% %     med = nanmedian(temp(:));
% %     
% %     idx = temp > (med+1.5*(hq-lq)) | temp < (med-1.5*(hq-lq));
% %     KLdiv_refDist(:, idx) = nan;
%     
%     
%     
%     for iEpoch = 1:3
%         currEpoch = epochNames{iEpoch};
%         
%         nTimeBins = size(PFcorrelation.(currEpoch){iSess}, 2);
%         
%         PFcorrelation.(currEpoch){iSess} = (PFcorrelation.(currEpoch){iSess} - repmat(nanmean(PFcorrelation_refDist, 2), [1, nTimeBins])) ...
%             ./repmat(nanstd(PFcorrelation_refDist, [], 2), [1, nTimeBins]);
%         
%         KLdiv.(currEpoch){iSess} = (KLdiv.(currEpoch){iSess} - repmat(nanmean(KLdiv_refDist, 2), [1, nTimeBins])) ...
%             ./repmat(nanstd(KLdiv_refDist, [], 2), [1, nTimeBins]); 
%         
%     end
% end
% 

PFcorrelation_pooled = cell(3,1);
KLdiv_pooled         = cell(3,1);

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    PFcorrelation_pooled{iEpoch} = cell2mat(PFcorrelation.(currEpoch));
    KLdiv_pooled{iEpoch}         = cell2mat(KLdiv.(currEpoch));
end





% LT-PF correlation

% cl = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.75, 0.75, 0]];

cl = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.75, 0.75, 0]; [0.4940, 0.1840, 0.5560]; [0, 0.5, 0]];



ax(1) = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');


hold on
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    currData = PFcorrelation_pooled{iEpoch};
    
    x = unifiedBinsCenters{iEpoch};
    
    if iEpoch == 2
%         currData = currData(:, 1:end-2);
%         x = x(1:end-2);
    end

    
    lq = nanprctile(currData, 25);
    hq = nanprctile(currData, 75);
    
    
    patch([x fliplr(x)], [lq' fliplr(hq')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        
    
    for iSess = 1:nSessions
        
        sessionData = PFcorrelation.(currEpoch){iSess};
        if iEpoch == 2
%             sessionData = sessionData(:, 1:end-2);
        end
        h(iSess) = plot(x, nanmedian(sessionData), 'color', [cl(iSess, :) 0.6], 'linewidth', 0.5); colormap('jet') ;
    end
    
    plot(x, nanmedian(currData), 'color', 'k', 'linewidth', 1); colormap('jet')

    
end

% legend(h, 'RatN', 'RatS', 'RatU', 'box', 'off', 'location', 'northout')
legend(h, 'Achilles-maze1', 'Achilles-maze2', 'Buddy', 'Cicero', 'Gatsby', 'box', 'off', 'location', 'northout')
% legend(h, 'Roy-Day1', 'Ted-Maze1', 'Ted-Maze2', 'Ted-Maze3', 'Kevin', 'box', 'off', 'location', 'northout')
% legend(h, 'Roy-Day1', 'Roy-Day2', 'Roy-Day3', 'Ted-Maze3', 'Kevin', 'box', 'off', 'location', 'northout')

ylabel({'LT-PF'; 'Pearson correlation coeff. %'}, 'fontsize', fontsize)

set(gca, 'position', sub_pos{1},'box', 'off',  'YTick', [0 25 50 75 100], 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 

xticklabels([])



axis tight
currAx = gca;
currAx.YGrid = 'on';


% LT-PF KL divergence

ax(2) = axes('position',sub_pos{2},'XGrid','off','XMinorGrid','off','FontSize', 6,'Box','off','Layer','top');


hold on
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    currData = KLdiv_pooled{iEpoch};
    
    x = unifiedBinsCenters{iEpoch};
    
    if iEpoch == 2
%         currData = currData(:, 1:end-2);
%         x = x(1:end-2);
    end

    
    lq = nanprctile(currData, 25);
    hq = nanprctile(currData, 75);
    
    
    
    patch([x fliplr(x)], [lq' fliplr(hq')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    
    for iSess = 1:nSessions
        
        sessionData = KLdiv.(currEpoch){iSess};
        if iEpoch == 2
%             sessionData = sessionData(:, 1:end-2);
        end
        
        plot(x, nanmedian(sessionData), 'color', [cl(iSess, :) 0.6], 'linewidth', 0.5); colormap('jet');
    end
    
    
    plot(x, nanmedian(currData), 'color', 'k', 'linewidth', 1); colormap('jet')

end


xlabel('time (hour)', 'fontsize', fontsize)
ylabel({'LT-PF'; 'KL divergence %'}, 'fontsize', fontsize)


set(gca, 'box', 'off', 'linewidth', 1, 'YTick', [0 25 50 75 100], 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 

currAx = gca;
currAx.YGrid = 'on';


linkaxes(ax, 'xy')

axis tight


ylim([0 100])

% xlim([0 9.6])
xlim([0 12])



%% sub-functions


function prct = nanprctile(data, percentile)


nTimeBins = size(data, 2);
prct = nan(nTimeBins, 1);

for iTimeBin = 1:nTimeBins
    
    currData = data(:, iTimeBin);
    
    currData = currData(~isnan(currData));
    prct(iTimeBin) = prctile(currData, percentile);

end

end


function KLdiv = calKLDivergence(learnedTunings, spatialTunings)


% the units should be in columns of learned tunings or spatialTunings
% matrix
% 
% if size(spatialTunings, 2) > size(spatialTunings, 1)
%     spatialTunings = spatialTunings';
% end

spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [size(spatialTunings, 1) 1]);
spatialTunings = spatialTunings + eps;

[~, nUnits] = size(spatialTunings);

% if size(learnedTunings, 1) ~= nPosBins
%     learnedTunings = learnedTunings';
% end

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [size(learnedTunings, 1) 1]);
learnedTunings = learnedTunings + eps;


KLdiv = nan(nUnits);

for iUnit = 1:nUnits
    currSpatialTuning = spatialTunings(:, iUnit);
    
    spatialTuningTerm = repmat(currSpatialTuning, [1 nUnits]);
    KLdiv(:, iUnit) = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);
end

end
