
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

assemblyTunings_time            = cell(nSessions, 1);



%% calculate matching scores for all sessions

iSess = 11;
    
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
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), ...
    'activeUnits', ...
    'assemblyTuningPFcorr', ...
    'assemblyTuningPFKLdiv');

activeUnits{iSess} = s.activeUnits;
okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);



% % learned tunings in 15-minutes time windows

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time_Jan2022.mat']), ...
        'binCenters', 'assemblyTunings_time');
    
binCenters{iSess}           = s.binCenters;
assemblyTunings_time{iSess} = s.assemblyTunings_time;

for iEpoch = 1:numel(epochNames)
    binCenters{iSess}.(epochNames{iEpoch}) = binCenters{iSess}.(epochNames{iEpoch})/3600;
end



% truncate the long recording sessions to the first 4 hours of post sleep

idx = find(binCenters{iSess}.post < endT{iSess}.post);

binCenters{iSess}.post = binCenters{iSess}.post(idx);
assemblyTunings_time{iSess}.post.data = assemblyTunings_time{iSess}.post.data(:, :, idx);




%% plot the learning tunings across time for example units in an example session


selectedUnits = 1:nUnits; 
% selectedUnits = [77 78 74 44 6 8 38 54 73 71 64 14 15];% Rat S
% selectedUnits = [7 15 17 37 49 56 68 70 73 16 82]; % Rat N
% selectedUnits = [89 126 163]; % RatU
% selectedUnits = [86 88 111 114 61]; % AChilles_10252013


spikes = spikes_pooled{iSess};


plotwidth  = 500;
plotheight = 800;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 30;
noc = 6;

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
    
    plot(0.5*currSpatialTuning + startT{iSess}.run , 1:nPosBins, 'linewidth', 1, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(ax(iUnit));

    
    if iUnit == 1
        line(0.5*[0 2] + startT{iSess}.run, [yl(2)-10 yl(2)-10], 'linewidth', 1, 'color', [0 0 0 0.6])
        text(startT{iSess}.run, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+12, ['unit ' num2str(iUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters{iSess}.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    if ismember(iUnit, [1])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters{iSess}.(currEpoch);
           text(median(currBinCenters), yl(2)+20, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if iUnit == 1
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, [20])
        xlabel('time (hour)')
    end
    
    if iUnit >= 21
        set(ax(iUnit), 'YTickLabel', {})
    end
    
    if ~ismember(iUnit, [20])
       set(ax(iUnit), 'XTickLabel', {}) 
    end
    
    set(gca, 'color', [0.5 0.5 0.5])
    
    hold off
    
end

