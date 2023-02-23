
% in this script, I pooled the learned tuning across all sessions and
% sorted the units in terms of their peak normalized positions 



clear; clc; close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

sub = dir(parentDir);
nSessions = numel(sub) - 2; % subtracting two because of '.' and '..'
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
    
    
    % track learned tunings
    runLearnedTunings_pooled  = cell(nSessions, 1);
    
    
    % assembly tunings
    assemblyTuning_pooled.(currEpoch).data = cell(nSessions, 1);
    assemblyTuning_pooled.(currEpoch).ui   = cell(nSessions, 1);
    

    % PF correlation
    assemblyTuningPFcorr_pooled.(currEpoch).data   = cell(nSessions, 1);
    assemblyTuningPFcorr_pooled.(currEpoch).ui     = cell(nSessions, 1);

    assemblyTuningPFcorr_zscore_pooled.(currEpoch) = cell(nSessions, 1);

end

gw = gausswindow(5,15);

flag = 1;


currSessions = setdiff(1:15, 11:12);

for iSession = 1:numel(currSessions)
    
    currSession = currSessions(iSession);
    
    
    sessionName = sub(currSession+2).name
    
    
    load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
    spikes = spikes_pyr;
    
    load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'nPBEs', 'activeUnits') % , 'assemblyTunings', 'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore'
    load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_1hour.mat']), 'assemblyTunings', 'assemblyTuningPFcorr')
    load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_activeRun.mat']), 'learnedTunings')
    
    
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
        
        for iUnit = 1:nUnits(iSession)
            
            currLearnedTuning = assemblyTunings.(currEpoch).data(sortedUnitIdx(iUnit), :);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;
            currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            currLearnedTuning = currLearnedTuning(startBin:endBin);
            
            assemblyTuning_pooled.(currEpoch).data{iSession}(iUnit, :) = currLearnedTuning/max(currLearnedTuning);
             
        end

        % PF correlation
        assemblyTuningPFcorr_pooled.(currEpoch).data{iSession}   = assemblyTuningPFcorr.(currEpoch).data(sortedUnitIdx);
    end
    
end

%%


plotheight = 11;
plotwidth  = 11;
fontsize   = 8;

nor = 1;
noc = 3;

leftmargin = 2;  rightmargin = 2;    topmargin = 2;    bottommargin = 2;

gapc = 1;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);




% for each session

plotwidth_epoch = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_epoch = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_epoch = 2;
noc_epoch = 2;

leftmargin_epoch = 0;  rightmargin_epoch = 0;     bottommargin_epoch = 0;    topmargin_epoch = 0;
gapr_epoch = 0.5;    gapc_epoch = 0.25;



cnctSpatialTunings = cell2mat(spatialTuning_pooled);

for iEpoch = 1:3

    epochName = epochNames{iEpoch};
    
        
    cnctAssemblyTunings = cell2mat(assemblyTuning_pooled.(epochName).data);
    cnctPFcorrelation   = cell2mat(assemblyTuningPFcorr_pooled.(epochName).data);

    
    hiLearnedTuning_unitIdx = find(cnctPFcorrelation >= prctile(cnctPFcorrelation, 50));
    hi_nUnits = numel(hiLearnedTuning_unitIdx);
    
    loLearnedTuning_unitIdx = find(cnctPFcorrelation < prctile(cnctPFcorrelation, 50));
    lo_nUnits = numel(loLearnedTuning_unitIdx);

    
    hiLearnedTunings = cnctAssemblyTunings(hiLearnedTuning_unitIdx, :);
    loLearnedTunings = cnctAssemblyTunings(loLearnedTuning_unitIdx, :);

    
    hiSpatialTunings = cnctSpatialTunings(hiLearnedTuning_unitIdx, :);
    loSpatialTunings = cnctSpatialTunings(loLearnedTuning_unitIdx, :);
    


    % sorting again because we are concatentating multiple sessions that
    % were sorted separately
    
    % best learned tunings
    [~, peakLocs] = max(hiSpatialTunings, [], 2);
    [~, sortIdx]  = sort(peakLocs, 'ascend');

    hiLearnedTunings = hiLearnedTunings(sortIdx, :);
    hiSpatialTunings = hiSpatialTunings(sortIdx, :);

    % worst learned tunings
    [~, peakLocs] = max(loSpatialTunings, [], 2);
    [~, sortIdx]  = sort(peakLocs, 'ascend');

    loLearnedTunings = loLearnedTunings(sortIdx, :);
    loSpatialTunings = loSpatialTunings(sortIdx, :);

    
    
    % %
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{iEpoch}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{iEpoch}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    
    ax1 = axes('position',sub_pos_epoch{1,1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    imagesc(hiSpatialTunings)
    xlim([0 nPosBins_interp])
    ylim([1 hi_nUnits])
    set(ax1, 'xTick', [], 'yTick', [1 hi_nUnits], 'YDir', 'normal', 'box', 'off')

    if iEpoch == 1
        h = text(-400, hi_nUnits/4, 'LT-PF corr > 50%', 'fontsize', fontsize);
        set(h, 'rotation', 90)
        
        h = text(-200, hi_nUnits/2, 'units', 'fontsize', fontsize);
        set(h, 'rotation', 90)

    end
    
    dim = [sub_pos_epoch{1,1}(1)-0.03 sub_pos_epoch{1,1}(2)+0.05 0.3 0.3];
    annotation('textbox',dim,'String','PFs' ,'FitBoxToText','on', 'EdgeColor', 'none', 'fontsize', fontsize)
    
    
%     title('place fields', 'fontsize', fontsize, 'fontweight', 'normal')
    
    ax2 = axes('position',sub_pos_epoch{1,2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    imagesc(hiLearnedTunings)
    xlim([0 nPosBins_interp])
    ylim([1 hi_nUnits])
    set(ax2, 'xTick', [], 'yTick', [], 'YDir', 'normal', 'box', 'off')
    
    
    dim = [sub_pos_epoch{1,2}(1)-0.03 sub_pos_epoch{1,1}(2)+0.05 0.3 0.3];
    annotation('textbox',dim,'String',[epochName ' PBEs'] ,'FitBoxToText','on', 'EdgeColor', 'none', 'fontsize', fontsize)
    
    
    ax3 = axes('position',sub_pos_epoch{2,1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    imagesc(loSpatialTunings)
    xlim([0 nPosBins_interp])
    ylim([1 hi_nUnits])
    set(ax3, 'xTick', [0 nPosBins_interp], 'xTickLabel', {'0'; '1'}, 'yTick', [1 hi_nUnits], 'YDir', 'normal', 'box', 'off')

    if iEpoch == 1
        h = text(-400, lo_nUnits/4, 'LT-PF corr. < 50%', 'fontsize', fontsize);
        set(h, 'rotation', 90)
        
        h = text(-200, lo_nUnits/2, 'units', 'fontsize', fontsize);
        set(h, 'rotation', 90)
    end
    xlabel('position', 'fontsize', fontsize)

    
    
    ax4 = axes('position',sub_pos_epoch{2,2},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    imagesc(loLearnedTunings)
    xlim([0 nPosBins_interp])
    ylim([1 hi_nUnits])
    set(ax4, 'xTick', [0 nPosBins_interp], 'xTickLabel', {'0'; '1'}, 'yTick', [], 'YDir', 'normal', 'box', 'off')
    
    xlabel('position', 'fontsize', fontsize)

    
end
dim = [sub_pos{2}(1)-0.05 sub_pos_epoch{1,1}(2)+0.15 0.3 0.3];
annotation('textbox',dim,'String','learned tunings' ,'FitBoxToText','on', 'EdgeColor', 'none')
    


colormap('jet')