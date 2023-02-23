clear; clc; close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);
nSessions    = numel(sub) - 2; % subtracting two becuase of '.' and '..'
% sessionNames = cell(nSessions, 1);


epochNames = {'pre'; 'run'; 'post'};


replayScoreMethods = {'rt_ui'; 'wc_ui'; 'wc_ts'};


replayScoreMethods_fullName = {'radon Integral - unit ID shuffle'; ...
                               'weighted Corr - unit ID shuffle' ; ...
                               'weighted Corr - wPBE time swap'};

sessionNames = {'Achilles-Maze1'; 'Achilles-Maze2'; 'Buddy'; 'Cicero'; 'Gatsby'; ...
                'Kevin';
                'RatN'; 'RatS'; 'RatU';
                'Roy-D1'; 'Roy-D2'; 'Roy-D3'; 'Ted-M1'; 'Ted-M2'; 'Ted-M3'};
                           
% initialize

nUnits = zeros(nSessions, 1);
peakPFlocation_session = cell(nSessions, 1);
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};

    % track spatial tunings
    spatialTuning_pooled = cell(nSessions, 1);
    
    
    % assembly tunings
    assemblyTuning_pooled.(currEpoch).data = cell(nSessions, 1);
    assemblyTuning_pooled.(currEpoch).ui   = cell(nSessions, 1);
    

    % PF correlation
    assemblyTuningPFcorr_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPFcorr_pooled.(currEpoch).ui     = cell(nSessions, 1);
    assemblyTuningPFcorr_zscore_pooled.(currEpoch) = cell(nSessions, 1);
    %
    assemblyTuningPFcorr_pscore_pooled.(currEpoch) = cell(nSessions, 1);

    
    % preffered location distance
    assemblyTuning_prefPos_pooled.(currEpoch).data = cell(nSessions, 1);
    assemblyTuning_prefPos_pooled.(currEpoch).ui = cell(nSessions, 1);
    
    assemblyTuningPF_prefPosDist_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPF_prefPosDist_pooled.(currEpoch).ui     = cell(nSessions, 1);
    assemblyTuningPF_prefPosDist_zscore_pooled.(currEpoch) = cell(nSessions, 1);
    %
    assemblyTuningPF_prefPosDist_pscore_pooled.(currEpoch) = cell(nSessions , 1);

    
    % KL divergence
    assemblyTuningPFKLdiv_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPFKLdiv_pooled.(currEpoch).ui     = cell(nSessions, 1);
    assemblyTuningPFKLdiv_zscore_pooled.(currEpoch) = cell(nSessions, 1);
    %
    assemblyTuningPFKLdiv_pscore_pooled.(currEpoch) = cell(nSessions, 1);

end

gw = gausswindow(5,15);

flag = 1;

% currSessions = 1:2;
currSessions = [1:5 7:9 10 13:15 6];

for iSession = 1:numel(currSessions)
    
    currSession = currSessions(iSession);
    
    
    sessionName = sub(currSession+2).name
    
    
    load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
    spikes = spikes_pyr;
    
    load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), ...
        'nPBEs', 'activeUnits', 'assemblyTunings', ...
        'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore', ...
        'assemblyTuning_prefPos', 'assemblyTuningPF_prefPosDist', 'assemblyTuningPF_prefPosDist_zscore', ...
        'assemblyTuningPFKLdiv', 'assemblyTuningPFKLdiv_zscore'); 
%     load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), ...
%         'nPBEs', 'activeUnits', 'assemblyTunings', ...
%         'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore', ...
%         'asTuning_prefPos', 'asTuning_PF_prefPosDist', 'asTuning_PF_prefPosDist_zscore', ...
%         'asTuning_PF_KLdiv', 'asTuning_PF_KLdiv_zscore'); 
%     assemblyTuning_prefPos = asTuning_prefPos;
%     assemblyTuningPF_prefPosDist = asTuning_PF_prefPosDist;
%     assemblyTuningPF_prefPosDist_zscore = asTuning_PF_prefPosDist_zscore;
%     assemblyTuningPFKLdiv = asTuning_PF_KLdiv;
%     assemblyTuningPFKLdiv_zscore = asTuning_PF_KLdiv_zscore;
    
   

    okUnits = intersect(activeUnits.pre, activeUnits.post);
    okUnits = intersect(okUnits, activeUnits.run);
    
    spikes = spikes(okUnits);
        
    
    % track spatial tunings
    nUnits(iSession) = numel(spikes); 
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    
    interpPosBins = linspace(0, nPosBins, 200);
    nPosBins_interp = numel(interpPosBins);
    
    
    % non-directional spatial tuning
    spatialTunings_merge = zeros(nUnits(iSession), nPosBins_interp);
    
    peakPFlocation    = zeros(nUnits(iSession), 1);
    peakPFfiring      = zeros(nUnits(iSession), 1);
    
    
    
    for iUnit = 1:nUnits(iSession)
        currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;
        currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
        currSpatialTuning = currSpatialTuning./max(currSpatialTuning);
        
        spatialTunings_merge(iUnit, :) = currSpatialTuning;
        [peakPFfiring(iUnit), peakPFlocation(iUnit)] = max(currSpatialTuning);            
    end
    
    
    [~, sortIdx] = sort(peakPFlocation, 'ascend');
    sortedUnitIdx = okUnits(sortIdx);
    peakPFlocation_session{iSession} = peakPFlocation(sortIdx);
    
    
    temp = sum(spatialTunings_merge, 1);
    if flag == 1
        startBin = find(~isnan(temp), 1, 'first');
        endBin   = find(~isnan(temp), 1, 'last');
        flag = 0;
    end
    
    spatialTuning_pooled{iSession} = spatialTunings_merge(sortIdx, startBin:endBin);
    nPosBins_interp = endBin - startBin + 1;
    
        
    for iEpoch = 1:3
        currEpoch = epochNames{iEpoch};
        
        
        % assembly Tunings
        assemblyTuning_pooled.(currEpoch).data{iSession} = zeros(nUnits(iSession), nPosBins_interp);
        leanedTunings_normPosition = zeros(nUnits(iSession), nPosBins_interp);
        
        for iUnit = 1:nUnits(iSession)
            
            currLearnedTuning = assemblyTunings.(currEpoch).data(sortedUnitIdx(iUnit), :);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;
            currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            currLearnedTuning = currLearnedTuning(startBin:endBin);
            
            leanedTunings_normPosition(iUnit, :) = currLearnedTuning;
            
            assemblyTuning_pooled.(currEpoch).data{iSession}(iUnit, :) = currLearnedTuning/max(currLearnedTuning);
             
        end

        % spatial bin correlation
        assemblyTuningPFcorr_pooled.(currEpoch).data{iSession}   = assemblyTuningPFcorr.(currEpoch).data(sortedUnitIdx); % data        
        assemblyTuningPFcorr_pooled.(currEpoch).ui{iSession}     = assemblyTuningPFcorr.(currEpoch).ui(sortedUnitIdx, :); % shuffle
        assemblyTuningPFcorr_zscore_pooled.(currEpoch){iSession} = assemblyTuningPFcorr_zscore.(currEpoch)(sortedUnitIdx, :); % z-score
        % 
        allPFcorr_prctiles = calPrctile(assemblyTuningPFcorr.(currEpoch));
        assemblyTuningPFcorr_pscore_pooled.(currEpoch){iSession} = allPFcorr_prctiles(sortedUnitIdx);
        
        % preferred position distance
        assemblyTuning_prefPos_pooled.(currEpoch).data{iSession}   = assemblyTuning_prefPos.(currEpoch).data(sortedUnitIdx);
        assemblyTuning_prefPos_pooled.(currEpoch).ui{iSession}     = assemblyTuning_prefPos.(currEpoch).ui(sortedUnitIdx, :);
        
        % preferred position distance betwene the LTs and the place fields
        % of identical units
        assemblyTuningPF_prefPosDist_pooled.(currEpoch).data{iSession}   = assemblyTuningPF_prefPosDist.(currEpoch).data(sortedUnitIdx); % data        
        assemblyTuningPF_prefPosDist_pooled.(currEpoch).ui{iSession}     = assemblyTuningPF_prefPosDist.(currEpoch).ui(sortedUnitIdx, :); % shuffle
        assemblyTuningPF_prefPosDist_zscore_pooled.(currEpoch){iSession} = assemblyTuningPF_prefPosDist_zscore.(currEpoch)(sortedUnitIdx, :); % z-score
        % 
        allPeakDistances_prctiles = calPrctile(assemblyTuningPF_prefPosDist.(currEpoch));
        assemblyTuningPF_prefPosDist_pscore_pooled.(currEpoch){iSession} = allPeakDistances_prctiles(sortedUnitIdx);
        

   
        % KL divergence
        assemblyTuningPFKLdiv_pooled.(currEpoch).data{iSession}   = assemblyTuningPFKLdiv.(currEpoch).data(sortedUnitIdx); % data        
        assemblyTuningPFKLdiv_pooled.(currEpoch).ui{iSession}     = assemblyTuningPFKLdiv.(currEpoch).ui(sortedUnitIdx, :); % shuffle
        assemblyTuningPFKLdiv_zscore_pooled.(currEpoch){iSession} = assemblyTuningPFKLdiv_zscore.(currEpoch)(sortedUnitIdx, :); % z-score
        
        % 
        allKLdiv_prctiles = calPrctile(assemblyTuningPFKLdiv.(currEpoch));
        assemblyTuningPFKLdiv_pscore_pooled.(currEpoch){iSession} = allKLdiv_prctiles(sortedUnitIdx);
        
    end
    
end

%% plot histograms of peak distance, spatial bin correlation, and KL divergence (panel C left column)

plotwidth  = 120;
plotheight = 280;
fontsize   = 6;

nor = 3;
noc = 1;

leftmargin = 20;  rightmargin = 10;    topmargin = 60;    bottommargin = 60;

gapc = 0;
gapr = 20;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% peak position distances bw LTs and PFs

ax1 = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customCDF(assemblyTuningPF_prefPosDist_pscore_pooled, 'right')

xlabel({'peak distance(normalized)'}, 'fontsize', fontsize)
ylabel('fraction of units', 'fontsize', fontsize)
% title('All PBEs in an epoch', 'fontsize', fontsize)


rr = legend(ax1, 'pre', 'run', 'post', 'box', 'off', 'location', 'northout');
rr.FontSize = fontsize;

set(ax1, 'box', 'off', 'position', sub_pos{1},'linewidth', 1)


% testing signifinace difference between the groups

rr = assemblyTuningPFKLdiv_pscore_pooled;
temp = [cell2mat(rr.pre) cell2mat(rr.run) cell2mat(rr.post)];

[a, b, c] = friedman(temp);


p_value = nan(3);
z_value = nan(3);

for iEpoch = 1:3
    for jEpoch = iEpoch+1:3

        [p_value(iEpoch, jEpoch), ~, stats] = signrank(cell2mat(rr.(epochNames{iEpoch})), cell2mat(rr.(epochNames{jEpoch})));
        z_value(iEpoch, jEpoch) = stats.zval;

    end
end



% Pearson correlation
ax2 = axes('position',sub_pos{2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customCDF(assemblyTuningPFcorr_pscore_pooled, 'left')


xlabel({'correlation %'}, 'fontsize', fontsize)
ylabel('fraction of units', 'fontsize', fontsize)

set(ax2, 'box', 'off', 'position', sub_pos{2},'linewidth', 1)



% KL divegence
ax3 = axes('position',sub_pos{3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customCDF(assemblyTuningPFKLdiv_pscore_pooled, 'right')

xlabel({'KL divergence %'}, 'fontsize', fontsize)
ylabel('fraction of units', 'fontsize', fontsize)

set(ax3, 'box', 'off', 'position', sub_pos{3},'linewidth', 1)



%% plot the distribution of the place field matching scores for individual sessions (panel C right column)

plotwidth  = 140;
plotheight = 300;
fontsize   = 6;

nor = 3;
noc = 1;

leftmargin = 30;  rightmargin = 10;    topmargin = 60;    bottommargin = 70;

gapc = 0;
gapr = 20;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% spatial bin correlation
ax1 = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customPlotIndividualSessionDist(assemblyTuningPF_prefPosDist_pscore_pooled)

ylim([0 100])

set(ax1, 'xtick', [1 2 3], 'xticklabel', [], 'box', 'off', 'linewidth', 1)
ylabel('peak distance', 'fontsize', fontsize)


% preffered position distance
ax2 = axes('position',sub_pos{2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customPlotIndividualSessionDist(assemblyTuningPFcorr_pscore_pooled)

ylim([0 100])

set(ax2, 'xtick', [1 2 3], 'xticklabel', [], 'box', 'off', 'linewidth', 1)
ylabel('correlation %', 'fontsize', fontsize)




% KL divegence
ax3 = axes('position',sub_pos{3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customPlotIndividualSessionDist(assemblyTuningPFKLdiv_pscore_pooled)

ylim([0 100])

set(ax3, 'xtick', [1 2 3], 'xticklabel', epochNames, 'box', 'off', 'linewidth', 1)
ylabel('KL divergence %', 'fontsize', fontsize)



%% panel D

plotwidth  = 200;
plotheight = 240;
fontsize   = 6;

nor = 2;
noc = 3;

leftmargin = 20;  rightmargin = 10;    topmargin = 60;    bottommargin = 60;

gapc = 20;
gapr = 20;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% scatter plot of spatial bin correlation

%
ax1 = axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

epoch1 = 'pre';
epoch2 = 'run';
customScatter(assemblyTuningPFcorr_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)
title('correlation', 'fontsize', fontsize, 'fontweight', 'normal')

set(ax1, 'box', 'off', 'position', sub_pos{1,1},'linewidth', 1, 'fontsize', fontsize)
axis square


%
ax2 = axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

epoch1 = 'run';
epoch2 = 'post';
customScatter(assemblyTuningPFcorr_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)

set(ax2, 'box', 'off', 'position', sub_pos{1,2},'linewidth', 1, 'fontsize', fontsize)
axis square

%
ax3 = axes('position',sub_pos{1,3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

epoch1 = 'pre';
epoch2 = 'post';
customScatter(assemblyTuningPFcorr_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)

set(ax3, 'box', 'off', 'position', sub_pos{1,3},'linewidth', 1, 'fontsize', fontsize)
axis square


% scatter plot of KL divergence divergence 
ax4 = axes('position',sub_pos{2,1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
epoch1 = 'pre';
epoch2 = 'run';
customScatter(assemblyTuningPFKLdiv_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)
title('KL divergence', 'fontsize', fontsize, 'fontweight', 'normal')

set(ax4, 'box', 'off', 'position', sub_pos{2,1},'linewidth', 1, 'fontsize', fontsize)
axis square


%
ax5 = axes('position',sub_pos{2,2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
epoch1 = 'run';
epoch2 = 'post';
customScatter(assemblyTuningPFKLdiv_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)

set(ax5, 'box', 'off', 'position', sub_pos{2,2},'linewidth', 1, 'fontsize', fontsize)
axis square


%
ax6 = axes('position',sub_pos{2,3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
epoch1 = 'pre';
epoch2 = 'post';
customScatter(assemblyTuningPFKLdiv_pscore_pooled, epoch1, epoch2)

xlabel(epoch1, 'fontsize', fontsize)
ylabel(epoch2, 'fontsize', fontsize)

set(ax6, 'box', 'off', 'position', sub_pos{2,3},'linewidth', 1, 'fontsize', fontsize)
axis square


%% histogram of peak preferred location (panel B, top row)


% place fields
cnctSpatialTunings = cell2mat(spatialTuning_pooled);
nUnits = size(cnctSpatialTunings, 1);

[~, cnctSpatialTuning_peakLocs] = max(cnctSpatialTunings, [], 2);


% learned tunings

for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};
    cnctAssemblyTunings = cell2mat(assemblyTuning_pooled.(currEpoch).data);

    [~, cnctLearnedTuning_peakLocs.(currEpoch)] = max(cnctAssemblyTunings, [], 2);
end



plotwidth  = 170;
plotheight = 170;
fontsize   = 6;

nor = 1;
noc = 4;

leftmargin = 20;  rightmargin = 10;    topmargin = 60;    bottommargin = 60;

gapc = 10;
gapr = 10;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



bins = linspace(0, 1, 20);

% place fields
ax1 = axes('position',sub_pos{1,1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customBars(cnctSpatialTuning_peakLocs/200, bins)

ylabel('fraction of units', 'fontsize', fontsize)
title('Place fields', 'fontsize', fontsize, 'fontweight', 'normal')

set(ax1, 'box', 'off','linewidth', 1, 'fontsize', fontsize)

 
% learned tunings
ax2 = axes('position',sub_pos{1,2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customBars(cnctLearnedTuning_peakLocs.pre/200, bins)

xlabel('peak position(normalized)', 'fontsize', fontsize)

title('pre', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax2, 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'ytick', [])


ax3 = axes('position',sub_pos{1,3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customBars(cnctLearnedTuning_peakLocs.run/200, bins)

title('run', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax3, 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'ytick', [])


ax4 = axes('position',sub_pos{1,4},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
customBars(cnctLearnedTuning_peakLocs.post/200, bins)

title('post', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax4, 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'ytick', [])


linkaxes([ax1, ax2, ax3, ax4], 'xy')





%% scatter of peak locations of preferred position (Panel B, bottom row)


plotwidth  = 200;
plotheight = 170;
fontsize   = 6;

nor = 1;
noc = 3;

leftmargin = 50;  rightmargin = 10;    topmargin = 60;    bottommargin = 60;

gapc = 10;
gapr = 10;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


ax1 = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customScatter2(cnctSpatialTuning_peakLocs/200, cnctLearnedTuning_peakLocs.pre/200)

ylabel({'learned tuning'; 'peak position(norm)'}, 'fontsize', fontsize)
title('pre', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax1, 'box', 'off','linewidth', 1, 'fontsize', fontsize)




ax2 = axes('position',sub_pos{2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customScatter2(cnctSpatialTuning_peakLocs/200, cnctLearnedTuning_peakLocs.run/200)

xlabel('PF peak position(norm)', 'fontsize', fontsize)
title('run', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax2, 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'yticklabel', [])




ax3 = axes('position',sub_pos{3},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

customScatter2(cnctSpatialTuning_peakLocs/200, cnctLearnedTuning_peakLocs.post/200)
title('run', 'fontsize', fontsize, 'fontweight', 'normal')
set(ax3, 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'yticklabel', [])




% ax8 = axes('position',sub_pos{2,4},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
% 
% customScatter2(cnctLearnedTuning_peakLocs.pre/200, cnctLearnedTuning_peakLocs.post/200)
% title('run', 'fontsize', fontsize, 'fontweight', 'normal')
% set(ax6, 'box', 'off','linewidth', 1, 'fontsize', fontsize)





%% functions


function customHistogram(variable)


epochNames = {'pre'; 'run'; 'post'};

if isfield(variable.pre, 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:3
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2

end 



allValues = [cell2mat(variable.pre); cell2mat(variable.run); cell2mat(variable.post)];
% allValues = allValues(allValues > prctile(allValues, 0.1) & allValues < prctile(allValues, 99.9));

bins = linspace(min(allValues), max(allValues)+0.01, 20);


% colors = {'b'; 'none';'r'};
% edgeColors = {'none'; 'k'; 'none'};
% medianColors = {'b'; 'k'; 'r'};

colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255];

hold on
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};
    pooledVar = cell2mat(variable.(currEpoch));     
    
    h = histc(pooledVar, bins); h(end) = [];
    
%     plot(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'color', [colors(iEpoch, :) 0.5], 'linewidth', 1, 'displayName', epochNames{iEpoch})
    area([0 bins(1:end-1)+diff(bins(1:2))/2], [0; h]/sum(h), 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.3, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    
%     a(iEpoch) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'Edgecolor', edgeColors{iEpoch}, 'facecolor', colors{iEpoch}, 'facealpha', 0.3);
%     hold on

    medianCorr.(currEpoch) = nanmedian(pooledVar);

end

xlim([min(allValues), max(allValues)])


rr = ylim;
for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};
    line([medianCorr.(currEpoch) medianCorr.(currEpoch)], [0 rr(2)], 'color', [colors(iEpoch, :) 0.6], 'linewidth', 2)
end


% significance scores
ranksumPvalue = nan(3);
for iEpoch = 1:3
   var1 = cell2mat(variable.(epochNames{iEpoch}));
   
   for jEpoch = setdiff(1:3, iEpoch)
       var2 = cell2mat(variable.(epochNames{jEpoch}));
       
       ranksumPvalue(iEpoch, jEpoch) = ranksum(var1, var2, 'tail', 'right');
       
   end
end
       
% add the significance signs here ...  
    
    
end
 


function customCDF(variable, tail)


epochNames = {'pre'; 'run'; 'post'};

if isfield(variable.pre, 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:3
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2
end 


colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255];

hold on
for iEpoch = 1:3

    currEpoch = epochNames{iEpoch};
    pooledVar = cell2mat(variable.(currEpoch));     
    
    nUnits = numel(pooledVar);

    [cdf_pooled, x_pooled] = ecdf(pooledVar);
    auc.(currEpoch) = trapz(x_pooled, cdf_pooled);

    [x_pooled, idx] = unique(x_pooled);
    cdf_pooled = cdf_pooled(idx);
    
%     area(x_pooled, cdf_pooled, 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.3, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    plot(x_pooled, cdf_pooled, 'color', colors(iEpoch, :), 'linewidth', 1)

    medianCorr.(currEpoch) = nanmedian(pooledVar);

end



if ~isempty(tail)
    nIters = 10000;
    auc_chance = nan(nIters, 1);
    
    
    for iIter = 1:nIters
        rndPrctls = randi(100, nUnits, 1);
    
        [cdf_chnace, x_chance] = ecdf(rndPrctls);
    
        auc_chance(iIter) = trapz(x_chance, cdf_chnace);
    end
end


ylim([0 1])
yl = ylim;
xl = xlim;


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    line([xl(1) medianCorr.(currEpoch)], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    line([medianCorr.(currEpoch) medianCorr.(currEpoch)], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    
    if ~isempty(tail)
        [~, pval] = ttest2(auc.(currEpoch), auc_chance, 'tail', tail);
        text(medianCorr.(currEpoch), 0.07, sprintf('%.2f(p=%.2e)', medianCorr.(currEpoch), pval), 'fontsize', 6, 'color', colors(iEpoch, :))
    end
end


% significance scores
signrankPvalue = nan(3);
for iEpoch = 1:3
   var1 = cell2mat(variable.(epochNames{iEpoch}));
   
   for jEpoch = setdiff(1:3, iEpoch)
       var2 = cell2mat(variable.(epochNames{jEpoch}));

       signrankPvalue(iEpoch, jEpoch) = signrank(var1, var2);
   end
end
       
% add the significance signs here ...  
    
    
end


function customScatter(variable, epoch1, epoch2)

epochNames = {'pre'; 'run'; 'post'};

if isfield(variable.pre, 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:3
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2

end 

for iEpoch = 1:3
    variable.(epochNames{iEpoch}) = cell2mat(variable.(epochNames{iEpoch}));
end

allValues = [variable.pre; variable.run; variable.post];
minVal = prctile(allValues, 0.1);
maxVal = prctile(allValues, 99.9);

idx = (variable.pre > minVal & variable.pre < maxVal) & (variable.run > minVal & variable.run < maxVal) & (variable.post > minVal & variable.post < maxVal);

for iEpoch = 1:3
    variable.(epochNames{iEpoch}) = variable.(epochNames{iEpoch})(idx);
end


var1 = variable.(epoch1);
var2 = variable.(epoch2);

scatter(var1, var2, 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

hold on
% 
% VAR1 = [ones(size(var1)) var1];
% b = VAR1\var2;
% 
% var2Calc = VAR1*b;
% var2Calc(isnan(var2Calc)) = interp1(var1(~isnan(var2Calc)), var2Calc(~isnan(var2Calc)), var1(isnan(var2Calc)));
% 
% plot(var1, var2Calc, 'r', 'linewidth', 1.5)
% 
% grid on
% 
% axis square
% 
% 
% 
% [corrValue, pvalue] = corr(var1, var2, 'type', 'pearson')
% 


end

function customBars(peakLocs, bins)

[h,b] = hist(peakLocs, bins);
h = h/sum(h);

bar(b, h, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])

axis square

end


function customScatter2(var1, var2)


scatter(var1, var2, 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

hold on

VAR1 = [ones(size(var1)) var1];
b = VAR1\var2;

var2Calc = VAR1*b;
% var2Calc(isnan(var2Calc)) = interp1(var1(~isnan(var2Calc)), var2Calc(~isnan(var2Calc)), var1(isnan(var2Calc)));

plot(var1, var2Calc, 'color',[1 0 0 0.5], 'linewidth', 1.5)

% grid on

axis square

[corrValue, pvalue] = corr(var1, var2, 'type', 'pearson');

text(0.5, 0.8, sprintf('rho=%.2f, p=%.2e', corrValue, pvalue), 'FontSize', 6)



end

function customPlotIndividualSessionDist(variable)


epochNames = {'pre'; 'run'; 'post'};
colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255];


if isfield(variable.pre, 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:3
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2

end 


hold on
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    epochData = variable.(currEpoch);
    
    nSessions = numel(epochData);
    
    for iSess = 1:nSessions
        
        currSessionData = epochData{iSess};
        currMed = nanmedian(currSessionData);
        curr_hq   = nanprctile(currSessionData, 75);
        curr_lq   = nanprctile(currSessionData, 25);
        
        
        curr_x = iEpoch+(iSess-ceil(nSessions/2))*1/nSessions;
        plot([curr_x curr_x], [curr_lq curr_hq], 'linewidth', 0.5, 'color', [colors(iEpoch, :) 0.6])
        
        scatter(curr_x, currMed, 5, 'MarkerFaceColor', colors(iEpoch, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.7)
        
    end
    
end

hold off

% xlim([0.5 3.5])


end

function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end


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