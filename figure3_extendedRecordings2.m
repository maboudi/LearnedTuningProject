
clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
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

assemblyTuningPFcorr        = cell(nSessions, 1);
assemblyTuningPFcorr_pscore = cell(nSessions, 1);
assemblyTuningPFcorr_time   = cell(nSessions, 1);


% KL  diveregenc (the same manner as the Pearson correlation)

assemblyTuningPFKLdiv        = cell(nSessions, 1);
assemblyTuningPFKLdiv_pscore = cell(nSessions, 1);
assemblyTuningPFKLdiv_time   = cell(nSessions, 1);


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

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), ...
        'binCenters', ...
        'assemblyTunings_time', ...
        'assemblyTuningPFcorr_time', ...
        'asTuning_PF_KLdiv_time', ...
        'assemblyTuningPFcorr_time_zscore', ...
        'asTuning_PF_KLdiv_time_zscore');
    
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

end



cnctSpatialTunings = cell2mat(spatialTuning);
[~, PFpeakLocs] = max(cnctSpatialTunings, [], 2);
[~, sortIdx]  = sort(PFpeakLocs, 'ascend');

cnctSpatialTunings_1 = cnctSpatialTunings(sortIdx, :);




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

PFcorrelation_pooled = cell(3,1);
KLdiv_pooled         = cell(3,1);


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
    
    PFcorrelation = cell(nSessions, 1);
    KLdiv         = cell(nSessions, 1);
    for iSess = 1:nSessions
        
        nUnits = size(assemblyTuningPFcorr_time{iSess}.(currEpoch), 1);
        
        PFcorrelation{iSess} = nan(nUnits, maxNBins);
        KLdiv{iSess}         = nan(nUnits, maxNBins);
        
        switch iEpoch
            case 1
                PFcorrelation{iSess}(:, end-allNBins(iSess)+1:maxNBins) = assemblyTuningPFcorr_time{iSess}.(currEpoch);
                KLdiv{iSess}(:, end-allNBins(iSess)+1:maxNBins)         = assemblyTuningPFKLdiv_time{iSess}.(currEpoch);
            case 2
                PFcorrelation{iSess}(:, 1:allNBins(iSess))              = assemblyTuningPFcorr_time{iSess}.(currEpoch);
                KLdiv{iSess}(:, 1:allNBins(iSess))                      = assemblyTuningPFKLdiv_time{iSess}.(currEpoch);
            case 3
                PFcorrelation{iSess}(:, 1:allNBins(iSess))              = assemblyTuningPFcorr_time{iSess}.(currEpoch);
                KLdiv{iSess}(:, 1:allNBins(iSess))                      = assemblyTuningPFKLdiv_time{iSess}.(currEpoch);
        end
    end
    
    PFcorrelation_pooled{iEpoch} = cell2mat(PFcorrelation);
    KLdiv_pooled{iEpoch}         = cell2mat(KLdiv);

end


% LT-PF correlation

ax(1) = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');


hold on
for iEpoch = 1:3
    
    currData = PFcorrelation_pooled{iEpoch};
    
    lq = nanprctile(currData, 25);
    hq = nanprctile(currData, 75);
    
    patch([unifiedBinsCenters{iEpoch} fliplr(unifiedBinsCenters{iEpoch})], [lq' fliplr(hq')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    
    plot(unifiedBinsCenters{iEpoch}, nanmedian(currData), 'color', 'k', 'linewidth', 1); colormap('jet')
    
end

ylabel({'LT-PF'; 'Pearson correlation coeff.'}, 'fontsize', fontsize)
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
xticklabels([])

xlim([0 13.5])


% axis tight
currAx = gca;
currAx.YGrid = 'on';


% LT-PF KL divergence

ax(2) = axes('position',sub_pos{2},'XGrid','off','XMinorGrid','off','FontSize', 6,'Box','off','Layer','top');


hold on
for iEpoch = 1:3
    
    currData = KLdiv_pooled{iEpoch};
    
    lq = nanprctile(currData, 25);
    hq = nanprctile(currData, 75);
    
    patch([unifiedBinsCenters{iEpoch} fliplr(unifiedBinsCenters{iEpoch})], [lq' fliplr(hq')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    
    plot(unifiedBinsCenters{iEpoch}, nanmedian(currData), 'color', 'k', 'linewidth', 1); colormap('jet')
    
end

xlabel('time (hour)', 'fontsize', fontsize)
ylabel({'LT-PF'; 'KL divergence'}, 'fontsize', fontsize)
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
currAx = gca;
currAx.YGrid = 'on';


xlim([0 13.5])


% linkaxes(ax, 'x')

% axis tight


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
