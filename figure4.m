
% currently it's figure 5

clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 12 15:17 6];
nSessions    = numel(currSessions);

epochNames = {'pre'; 'run'; 'post'};
curr_rsm = 'wc_ts';




activeUnits = cell(nSessions, 1);
nUnits = nan(nSessions, 1);

spatialTuning = cell(nSessions, 1);


learnedTuning_sub = cell(nSessions, 1);


learnedTuningPFcorr_sub         = cell(nSessions, 1);
learnedTuningPFcorr_sub_pscore  = cell(nSessions, 1);

learnedTuningPFKLdiv_sub        = cell(nSessions, 1);
learnedTuningPFKLdiv_sub_pscore = cell(nSessions, 1);



gw = gausswindow(3,9);

flag = 1;


for iSess = 1:nSessions
    
sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name
basePath      = fullfile(parentDir, sessionName);
    

load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
spikes = spikes_pyr;



% % learned tunings calculated based on the entire epochs
s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), 'activeUnits');

activeUnits{iSess} = s.activeUnits;

okUnits = intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post);
okUnits = intersect(okUnits, activeUnits{iSess}.run);

nUnits(iSess) = numel(okUnits);



% % learned tunings for different subset of PBEs

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores_Lthresh1e_3_wcts_lowMaxJump.mat']), ... %radonUI
    'assemblyTunings_sub', 'asTuningPFcorr_sub', 'asTuningPFKLdiv_sub', ...
    'asTuning_prefPos_sub', 'asTuningPF_prefPosDist_sub');


replayScoreMethods = fieldnames(s.asTuningPFcorr_sub.pre);



spikes = spikes(okUnits);

% track spatial tunings
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

interpPosBins = linspace(0, nPosBins, 200);
nPosBins_interp = numel(interpPosBins);


% non-directional spatial tuning
spatialTunings_merge = zeros(nUnits(iSess), nPosBins_interp);

peakPFlocation    = zeros(nUnits(iSess), 1);
peakPFfiring      = zeros(nUnits(iSess), 1);


for iUnit = 1:nUnits(iSess)

    currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;
    currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
    currSpatialTuning = currSpatialTuning./max(currSpatialTuning);

    spatialTunings_merge(iUnit, :) = currSpatialTuning;
    [peakPFfiring(iUnit), peakPFlocation(iUnit)] = max(currSpatialTuning);

end

[~, sortIdx] = sort(peakPFlocation, 'ascend');
sortedUnitIdx = okUnits(sortIdx);
% peakPFlocation_session{iSess} = peakPFlocation(sortIdx);


temp = sum(spatialTunings_merge, 1);
if flag == 1
    startBin = find(~isnan(temp), 1, 'first');
    endBin   = find(~isnan(temp), 1, 'last');
    flag = 0;
end

spatialTuning{iSess} = spatialTunings_merge(sortIdx, startBin:endBin);

nPosBins_interp = endBin - startBin + 1;
nShuffles = 10000;
        
for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};


    % assembly Tunings
    learnedTuning_sub{iSess}.(currEpoch)          = nan(nUnits(iSess), nPosBins_interp, 4);
    learnedTuningPFcorr_sub{iSess}.(currEpoch).ui = nan(nUnits(iSess), nShuffles, 4);
    
    for iSub = 1:4
        for iUnit = 1:nUnits(iSess)

            currLearnedTuning = s.assemblyTunings_sub.(currEpoch).(curr_rsm).data(sortedUnitIdx(iUnit), :, iSub);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;
            currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            currLearnedTuning = currLearnedTuning(startBin:endBin);

            learnedTuning_sub{iSess}.(currEpoch)(iUnit, :, iSub) = currLearnedTuning/max(currLearnedTuning);
        end

        % was added later
        
        shuffle_PFcorr = nan(nUnits(iSess), nShuffles);
        allCorrs = corr(learnedTuning_sub{iSess}.(currEpoch)(:, :, iSub)', spatialTuning{iSess}'); % already sorted so no need further

        for i_s = 1:nShuffles         
            cidx = randi(nUnits(iSess), nUnits(iSess), 1);
            shuffle_PFcorr(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:nUnits(iSess))', cidx));
        end
        learnedTuningPFcorr_sub{iSess}.(currEpoch).ui(:, :, iSub) = shuffle_PFcorr;
        learnedTuningPFcorr_sub{iSess}.(currEpoch).data(:, iSub)  = diag(allCorrs);
    end
    
    
%     learnedTuningPFcorr_sub{iSess}.(currEpoch).data   = s.asTuningPFcorr_sub.(currEpoch).(curr_rsm).data(sortedUnitIdx, :);
%     learnedTuningPFKLdiv_sub{iSess}.(currEpoch)  = s.asTuningPFKLdiv_sub.(currEpoch).(curr_rsm).data(sortedUnitIdx, :);


    allPFcorr = calPrctile(s.asTuningPFcorr_sub.(currEpoch).(curr_rsm));
    learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch) = allPFcorr(sortedUnitIdx, :);

%     allPFKLdiv = calPrctile(s.asTuningPFKLdiv_sub.(currEpoch).(curr_rsm));
%     learnedTuningPFKLdiv_sub_pscore{iSess}.(currEpoch) = allPFKLdiv(sortedUnitIdx, :);

    
end


end


%% plot the learning tunings for different quartiles of replay score

plotheight = 300;
plotwidth  = 250;
% fontsize   = 8;

nor = 1;
noc = 3;

leftmargin = 25;  rightmargin = 25;    topmargin = 25;    bottommargin = 25;

gapc = 10;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% for each epoch

plotwidth_epoch = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_epoch = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_epoch = 1;
noc_epoch = 2;

leftmargin_epoch = 0;  rightmargin_epoch = 0;     bottommargin_epoch = 0;    topmargin_epoch = 0;
gapr_epoch = 0;    gapc_epoch = 3;



cnctSpatialTunings = cell2mat(spatialTuning);
[~, PFpeakLocs] = max(cnctSpatialTunings, [], 2);
[~, sortIdx]  = sort(PFpeakLocs, 'ascend');

cnctSpatialTunings_1 = cnctSpatialTunings(sortIdx, :);


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    currLearnedTunings = cell(nSessions, 1);
    for iSess = 1:nSessions      
        currLearnedTunings{iSess} = learnedTuning_sub{iSess}.(currEpoch);
    end
    pooledLT.(currEpoch) = cell2mat(currLearnedTunings);
    
    % sorting again because we are concatentating multiple sessions that
    % were sorted separately
    

    pooledLT.(currEpoch) = pooledLT.(currEpoch)(sortIdx, :, :);
%     cnctSpatialTunings_1 = cnctSpatialTunings(sortIdx, :);
    
    tn_units = size(pooledLT.(currEpoch), 1);
    
    
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{iEpoch}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{iEpoch}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    
    
    for ii=1:2
        
        if ii == 1; iSub = 1;
        elseif ii == 2; iSub = 4;
        end
        
        ax1 = axes('position',sub_pos_epoch{ii},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

        imagesc(pooledLT.(currEpoch)(:,:,iSub))
        colormap('jet')
        xlim([0 nPosBins_interp])
        ylim([1 tn_units])
        set(ax1, 'xTick', [], 'yTick', [], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'linewidth', 1)
    end
    
end



%% PV spatial bin correlations

plotheight = 100;
plotwidth  = 250;
fontsize   = 6;

nor = 1;
noc = 3;

leftmargin = 25;  rightmargin = 25;    topmargin = 25;    bottommargin = 45;

gapc = 10;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);



% for each epoch

plotwidth_epoch = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
plotheight_epoch = (plotheight - (topmargin+bottommargin) - (nor-1) * gapr)/nor;

nor_epoch = 1;
noc_epoch = 2;

leftmargin_epoch = 0;  rightmargin_epoch = 0;     bottommargin_epoch = 0;    topmargin_epoch = 0;
gapr_epoch = 0;    gapc_epoch = 3;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



nEpochs = 1;
for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{nEpochs}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{nEpochs}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    

%     nSub = 1;
    for ii= 1:2

        if ii == 1; iSub = 1;
        elseif ii == 2; iSub = 4;
        end
    

        ax1 = axes('position',sub_pos_epoch{ii},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        
        currLTs = pooledLT.(currEpoch)(:, :, iSub);
        currLTs(isnan(currLTs)) = 0;

        inclusionIdx = sum(currLTs, 2) ~= 0;
        
        corrMat = corr(cnctSpatialTunings_1(inclusionIdx, :), currLTs(inclusionIdx, :));

        imagesc(corrMat)
        caxis([0 0.5])

        colormap('jet')
        xlim([0 nPosBins_interp])
        ylim([0 nPosBins_interp])
        
        
        if iSub == 1 && iEpoch == 1
            
            xlabel('learned tuning', 'fontsize', fontsize)
            ylabel('place field', 'fontsize', fontsize)

        end
    
        
        if iSub == 2 && iEpoch == 3

            h = colorbar('location', 'southoutside', 'TickDirection', 'out', 'FontSize', fontsize, 'LineWidth', 1);
            h.Label.String = 'PV correlation';
        end

        
        set(ax1, 'xTick', [0 nPosBins_interp], 'xTickLabel', {'0'; '1'}, 'yTick', [0 nPosBins_interp], ...
            'YTickLabel', {'0'; '1'}, 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.05, 0.05], 'fontsize', fontsize, 'linewidth', 1, 'position', sub_pos_epoch{ii})

        if iSub ~= 1 || iEpoch ~= 1
    
            set(ax1, 'YTickLabel', [])
        end

%         if iSub == 1
%             tl = 'low score';
%         elseif iSub == 4
%             tl = 'high score';
%         end

%         title(tl, 'fontsize', fontsize, 'FontWeight', 'normal')

        title(sprintf('%s-Q%d', currEpoch, iSub), 'fontsize', fontsize, 'FontWeight', 'normal')
        
        axis square
        ax1.XGrid = 'on';

%         nSub = nSub + 1;
   
    end

    nEpochs = nEpochs + 1;
end





%% Distribution of peak locations 
% and the relationship between the peak locations' distance from the track
% ends and PF correlation across units


for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    
    currLearnedTunings = cell(nSessions, 1);
    for iSess = 1:nSessions      
        currLearnedTunings{iSess} = learnedTuning_sub{iSess}.(currEpoch);
    end
    cnctLearnedTunings = cell2mat(currLearnedTunings);
    
    
    peakLocs = nan(size(cnctLearnedTunings, 1), 4);
    for iSub = 1:4
        [~, peakLocs(:,iSub)] = max(cnctLearnedTunings(:, :, iSub), [], 2);
    end
    
    learnedTuningPeakLocations.(currEpoch) = peakLocs;
    clear peakLocs
    
    
    % the PF correlation
    PFcorr = cell(nSessions, 1);
    for iSess = 1:nSessions        
        PFcorr{iSess} = learnedTuningPFcorr_sub{iSess}.(currEpoch);
    end
    PFcorr_data.(currEpoch) = cell2mat(PFcorr);
    
end



% plot the distribution of learned tuning's peak locations

plotwidth  = 250;
plotheight = 170;
fontsize   = 6;

nor = 1;
noc = 6;

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
quartiles2plot = [1 4];

ax = nan(6,1);

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    for ii = 1:2
        currQuartile = quartiles2plot(ii);
        
        isubplot = (iEpoch-1)*2 + ii;
        ax(isubplot) = axes('position',sub_pos{isubplot},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

        customBars(learnedTuningPeakLocations.(currEpoch)(:, currQuartile)/197, bins)

        if isubplot == 1
            
            xlabel('peak position(normalized)', 'fontsize', fontsize)
            ylabel('fraction of units', 'fontsize', fontsize)
        end
        
        set(ax(isubplot), 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01])
        
        if isubplot > 1
            set(ax(isubplot), 'YTickLabel', [])
        end
    end
end

linkaxes(ax, 'xy')
      

%% plot the distribution of learned tuning correlation versus LT peak location distance from
% the track ends

plotwidth  = 250;
plotheight = 170;
fontsize   = 6;

nor = 1;
noc = 6;

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


quartiles2plot = [1 4];

ax = nan(6,1);

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    for ii = 1:2
        currQuartile = quartiles2plot(ii);
        
        isubplot = (iEpoch-1)*2 + ii;
        ax(isubplot) = axes('position',sub_pos{isubplot},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        
%         currPeakLocs = learnedTuningPeakLocations.(currEpoch)(:, currQuartile)/197;
        currPeakLocs = PFpeakLocs/197;

        distanceFromTrackEnd = min([currPeakLocs 1-currPeakLocs], [], 2);
        
        currPFcorrs = PFcorr_data.(currEpoch)(:, currQuartile);
        scatter(distanceFromTrackEnd, currPFcorrs, 1, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
        
        
        if isubplot == 1
            xlabel('peak position(normalized)', 'fontsize', fontsize)
        end
        
        if isubplot == 1
            ylabel('PFcorrelation', 'fontsize', fontsize)
        end
        
        set(ax(isubplot), 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01])
        
        if isubplot > 1
            set(ax(isubplot), 'YTickLabel', [])
        end
        
        axis square
    end
end

linkaxes(ax, 'xy')




%% a barplot of PF correlation for the first 1/3rd (track end locations) and remaining 2/3rd as track middle locations (sould combine this with the previous plot)

plotwidth  = 270;
plotheight = 170;
fontsize   = 6;


nor = 1;
noc = 6;

leftmargin = 40;  rightmargin = 10;    topmargin = 60;    bottommargin = 60;


gapc = 10;
gapr = 10;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


quartiles2plot = [1 4];

ax = nan(6,1);

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};
    PFcorrs_segements.(currEpoch) = cell(2,2);
    

    signrank_pval_trackEnds.(currEpoch)   = nan(2,1);
    signrank_pval_trackMiddle.(currEpoch) = nan(2,1);
    trackEndCorrs.(currEpoch)             = cell(2,1);
    trackMidCorrs.(currEpoch)             = cell(2,1);

    for ii = 1:2
        currQuartile = quartiles2plot(ii);
        
        isubplot = (iEpoch-1)*2 + ii;
        ax(isubplot) = axes('position', sub_pos{isubplot},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
        
%         currPeakLocs = learnedTuningPeakLocations.(currEpoch)(:, currQuartile)/197;
        currPeakLocs = PFpeakLocs/197;
        distanceFromTrackEnd = min([currPeakLocs 1-currPeakLocs], [], 2);
        
        
        trackLocationIndex = zeros(size(distanceFromTrackEnd));
        trackLocationIndex(distanceFromTrackEnd > 0.167) = 1;
        
        
        
        currPFcorrs = PFcorr_data.(currEpoch)(:, currQuartile);
        boxplot(currPFcorrs, trackLocationIndex, 'boxstyle', 'filled', 'colors', [0.3 0.3 0.3]);
        

        % statistical difference from zero


        trackEndCorrs.(currEpoch){ii} = currPFcorrs(trackLocationIndex == 0);
        trackMidCorrs.(currEpoch){ii} = currPFcorrs(trackLocationIndex == 1);
        

        signrank_pval_trackEnds.(currEpoch)(ii) = signrank(trackEndCorrs.(currEpoch){ii}, 0, 'tail', 'right');
        signrank_pval_trackMiddle.(currEpoch)(ii) = signrank(trackMidCorrs.(currEpoch){ii}, 0, 'tail', 'right');


        
        
        if isubplot == 1
            xlabel('peak position(normalized)', 'fontsize', fontsize)
        end
        
        if isubplot == 1
            ylabel({'LT-PF'; 'Pearson correlation coeff.'}, 'fontsize', fontsize)
        end
        
        set(ax(isubplot), 'box', 'off','linewidth', 1, 'fontsize', fontsize, 'XTickLabel', {'end'; 'middle'}, ...
            'position', sub_pos{isubplot}, 'TickDir', 'out','TickLength',[0.01, 0.01])
        
        xtickangle(ax(isubplot), 45)
        
        set(ax(isubplot), 'YTick', -1:0.5:1)
        
        if isubplot > 1
            set(ax(isubplot), 'YTickLabel', [])
        end
        
        ylim([-1 1])
        
%         axis square
        PFcorrs_segements.(currEpoch){ii, 1} = currPFcorrs(trackLocationIndex == 0);
        PFcorrs_segements.(currEpoch){ii, 2} = currPFcorrs(trackLocationIndex == 1);
        
        h = gca;
        h.YGrid = 'on';
        
    end
    
    trackEnd_bw_PBEsubsets.(currEpoch) = signrank(trackEndCorrs.(currEpoch){2}, trackEndCorrs.(currEpoch){1}, 'tail', 'right');
    trackMid_bw_PBEsubsets.(currEpoch) = signrank(trackMidCorrs.(currEpoch){2}, trackMidCorrs.(currEpoch){1}, 'tail', 'right');
    

end

linkaxes(ax, 'xy')



%%

% plot the PF-matching scores for learned tunings calculated using
% PBEs in different quartiles of replay score

plotheight = 180;
plotwidth  = 210;
fontsize   = 6;

nor = 1;
noc = 3;

leftmargin = 40;  rightmargin = 40;    topmargin = 40;    bottommargin = 80;

gapc = 7;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



for iEpoch = 2
    
    
    ax(iEpoch) = axes('position',sub_pos{iEpoch},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    
    currEpoch = epochNames{iEpoch};
    
    PFcorr.data = cell(nSessions, 1);
    PFcorr.ui   = cell(nSessions, 1);


    KLdiv  = cell(nSessions, 1);
    
    PFcorr_pscore = cell(nSessions, 1);
    KLdiv_pscore  = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        PFcorr.data{iSess} = learnedTuningPFcorr_sub{iSess}.(currEpoch).data;
        PFcorr.ui{iSess}   = learnedTuningPFcorr_sub{iSess}.(currEpoch).ui;

%         KLdiv{iSess}  = learnedTuningPFKLdiv_sub{iSess}.(currEpoch);
        
        PFcorr_pscore{iSess} = learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch);
%         KLdiv_pscore{iSess}  = learnedTuningPFKLdiv_sub_pscore{iSess}.(currEpoch);
         
    end
    
    PFcorr_data  = cell2mat(PFcorr.data);
    PFcorr_ui    = cell2mat(PFcorr.ui); 

%     KLdiv_pooled  = cell2mat(KLdiv);
    
    PFcorr_pscore_pooled = cell2mat(PFcorr_pscore);
%     KLdiv_pscore_pooled  = cell2mat(KLdiv_pscore);
    
    
    currentVariable.(currEpoch) = PFcorr_data;

    pvalue.(currEpoch)   = nan(4,1);
    med_data.(currEpoch) = nan(4,1);

    for iSub = 1:4
        
        med_data.(currEpoch)(iSub) = nanmedian(PFcorr_data(:, iSub));
        med_ui   = nanmedian(PFcorr_ui(:,:, iSub));

        pvalue.(currEpoch)(iSub) = numel(find(med_ui >= med_data.(currEpoch)(iSub)))/size(med_ui, 2);
    end        
    
    
    h = violinplot(currentVariable.(currEpoch));

    for ii = 1:4

        h(ii).ViolinColor = [0 0 0];
        h(ii).ViolinAlpha = 0.7;
        h(ii).EdgeColor   = 'none';
        h(ii).ShowData    = 0;
        h(ii).MedianPlot.SizeData = 2.5;
%         h(ii).MedianPlot.Marker = 'square';
        h(ii).MedianPlot.MarkerEdgeColor = 'none';
        h(ii).MedianPlot.MarkerFaceColor = [1 1 1];
        h(ii).WhiskerPlot.LineWidth = 1;
        h(ii).WhiskerPlot.Color = [1 1 1];

    end
    

    xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})

    if iEpoch == 1
        ylabel('Pearson correlation', 'fontsize', fontsize)
    end

    if iEpoch == 2
        xlabel('replay scores', 'fontsize', fontsize)
    end
    
    xtickangle(45)
    

    
%     ylim([-1 1])
    
%     ylim([-0.1 2])
    currAx = gca;
    currAx.YGrid = 'on';
    
    if iEpoch >= 2
        currAx.YAxis.Visible = 'off';
    end
    

    currentVariable.(currEpoch) = currentVariable.(currEpoch)(~isnan(sum(currentVariable.(currEpoch), 2)), :);

    [pval(iEpoch), AnovaTab, stats] = friedman(currentVariable.(currEpoch));

    

%     [c, m, h, gnames] = multcompare(stats);
    

%     yl = ylim;
%     text(0, yl(2), sprintf('friedman p = %.0e', pval), 'fontsize', 6)
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
%     set(gca, 'YTick', 0:0.5:2)
    
end

 
linkaxes(ax, 'y')


% comparing the high score group between pre and post

[curr_p, ~,stats] = signrank(currentVariable.pre(:, 4), currentVariable.post(:, 4));



% comparing the magnitude of increase between low and high replay score groups 

postDiff = diff(currentVariable.post(:, [1 4]), [], 2);
preDiff  = diff(currentVariable.pre(:, [1 4]), [], 2);

statsPostDiff = [prctile(postDiff, 25) median(postDiff) prctile(postDiff, 75)];
statsPreDiff  = [prctile(preDiff, 25) median(preDiff) prctile(preDiff, 75)];

[curr_p, ~,stats] = signrank(postDiff, preDiff);





%% plot CDFs of LTPF percentile scores instead of pdf of actual LTPF matching scores

plotheight = 165;
plotwidth  = 210;
fontsize   = 6;

nor = 1;
noc = 3;

leftmargin = 40;  rightmargin = 40;    topmargin = 80;    bottommargin = 40;

gapc = 7;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



for iEpoch = 1:3
    
    
    ax(iEpoch) = axes('position',sub_pos{iEpoch},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    
    currEpoch = epochNames{iEpoch};
    
    PFcorr = cell(nSessions, 1);
    KLdiv  = cell(nSessions, 1);
    
    PFcorr_pscore = cell(nSessions, 1);
    KLdiv_pscore  = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        PFcorr_pscore{iSess} = learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch);
        KLdiv_pscore{iSess}  = learnedTuningPFKLdiv_sub_pscore{iSess}.(currEpoch);
         
    end
    
    PFcorr_pscore_pooled = cell2mat(PFcorr_pscore);
    KLdiv_pscore_pooled  = cell2mat(KLdiv_pscore);
    
    
    currentVariable = KLdiv_pscore_pooled;
    
    h = customCDF(currentVariable);


    if iEpoch == 3
        legend(h, {'0-25%'; '25-50%'; '50-75%'; '75-100%'}, 'Location', 'northoutside')
    end

    if iEpoch == 1
        ylabel('fraction of units', 'fontsize', fontsize)
    end

    if iEpoch == 2
        xlabel('LT-PF KL divergence', 'fontsize', fontsize)
    end
    

    
%     ylim([-1 1])
    
%     ylim([-0.1 2])
    currAx = gca;
%     currAx.YGrid = 'on';
    
    if iEpoch >= 2
        currAx.YAxis.Visible = 'off';
    end
    
    pval = kruskalwallis(currentVariable,[],'off');
    
    yl = ylim;
    text(0, yl(2), sprintf('kw pval = %.0e', pval), 'fontsize', 4)
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01], 'position', sub_pos{iEpoch}) 
%     set(gca, 'YTick', 0:0.5:2)
    
end
   
linkaxes(ax, 'y')




%% plot the distribution of LTPF percentile scores for individual sessions

plotheight = 220;
plotwidth  = 160;
fontsize   = 6;

nor = 3;
noc = 1;

leftmargin = 40;  rightmargin = 40;    topmargin = 80;    bottommargin = 40;

gapc = 0;
gapr = 7;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


cl = [0 0 1; 0 1 0; 1 0 0];
for iEpoch = 1:3
    
    
    ax(iEpoch) = axes('position',sub_pos{iEpoch},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    
    currEpoch = epochNames{iEpoch};
    
    PFcorr = cell(nSessions, 1);
    KLdiv  = cell(nSessions, 1);
    
    PFcorr_pscore = cell(nSessions, 1);
    KLdiv_pscore  = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        PFcorr_pscore{iSess} = learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch);
        KLdiv_pscore{iSess}  = learnedTuningPFKLdiv_sub_pscore{iSess}.(currEpoch);
         
    end
    
    
    
    currentVariable = KLdiv_pscore;
    
    customPlotIndividualSessionDist(currentVariable, cl(iEpoch, :))


    if iEpoch == 2
        ylabel('LT-PF KL divergence', 'fontsize', fontsize)
    end

    
%     ylim([-1 1])
    
%     ylim([-0.1 2])
    currAx = gca;
%     currAx.YGrid = 'on';
    

    currAx.XAxis.Visible = 'off';

        
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01], 'position', sub_pos{iEpoch}) 
%     set(gca, 'YTick', 0:0.5:2)
    
end
   

linkaxes(ax, 'y')





%% sub-functions


function prctileScore = calPrctile(variable)

[nUnits, nSubs, nShuffles] = size(variable.ui);
prctileScore = zeros(nUnits, nSubs);

for iUnit = 1:nUnits
    
    for iSub = 1:nSubs

        if isnan(variable.data(iUnit, iSub))
            prctileScore(iUnit, iSub) = nan;
        else
            prctileScore(iUnit, iSub) = numel(find(variable.ui(iUnit, iSub,  :) <= variable.data(iUnit, iSub)))/nShuffles * 100;
        end
    end
    
end

end



function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end




function customBars(peakLocs, bins)

[h,b] = hist(peakLocs, bins);
h = h/sum(h);

bar(b, h, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5])

axis square

end




function h = customCDF(variable)

[nUnits, nGrps] = size(variable);

medianCorr = nan(nGrps, 1);
auc        = nan(nGrps, 1);

hold on
for iGrp = 1:nGrps

    currScores = variable(:, iGrp);     
    
    [cdf_pooled, x_pooled] = ecdf(currScores);
    auc(iGrp) = trapz(x_pooled, cdf_pooled);
    
    [x_pooled, idx] = unique(x_pooled);
    cdf_pooled = cdf_pooled(idx);
    
%     area(x_pooled, cdf_pooled, 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.3, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    h(iGrp) = plot(x_pooled, cdf_pooled, 'color', [0 0 0 0.4+(iGrp-1)*0.15], 'linewidth', 1);

    medianCorr(iGrp) = nanmedian(currScores);
end

ylim([0 1])
yl = ylim;
xl = xlim;


%plotting the median and also comparing the area under curve to
% chance/uniform CDFs

% random number generator simulating uniform CDFs


nIters = 10000;
auc_chance = nan(nIters, 1);


for iIter = 1:nIters
    rndPrctls = randi(100, nUnits, 1);

    [cdf_chnace, x_chance] = ecdf(rndPrctls);

    auc_chance(iIter) = trapz(x_chance, cdf_chnace);
end




for iGrp = 1:nGrps
    
    [~, pval] = ttest2(auc(iGrp), auc_chance, 'tail', 'right');


    line([xl(1) medianCorr(iGrp)], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [0 0 0 0.4+(iGrp-1)*0.15])
    line([medianCorr(iGrp) medianCorr(iGrp)], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [0 0 0 0.4+(iGrp-1)*0.15])
    text(medianCorr(iGrp), 0.07+(iGrp-1)*0.1, sprintf('%d(p=%.2e)', floor(medianCorr(iGrp)), pval), 'fontsize', 4, 'color', [0 0 0 0.4+(iGrp-1)*0.15])

end



% significance scores
ranksumPvalue = nan(nGrps);

for iGrp = 1:nGrps
   var1 = variable(:, iGrp);
   
   for jGrp = setdiff(1:nGrps, iGrp)
       var2 = variable(:, jGrp);
       ranksumPvalue(iGrp, jGrp) = ranksum(var1, var2, 'tail', 'right');
   end
end
       
% add the significance signs here ...  
    
    
end


function customPlotIndividualSessionDist(variable, cl)


nGrps = 4; % quartiles of replay score

hold on

nSessions = numel(variable);

for iSess = 1:nSessions
    
    currSessVar = variable{iSess};

    for iGrp = 1:nGrps
        currGroup = currSessVar(:, iGrp);

        currMed = nanmedian(currGroup);
        curr_hq = nanprctile(currGroup, 75);
        curr_lq = nanprctile(currGroup, 25);
        
        curr_x = iSess+(iGrp-ceil(nGrps))*1/(1.2*nGrps);
      
        h= line([curr_x curr_x], [curr_lq curr_hq], 'linewidth', 0.5);
        h.Color = [cl 0.4+(iGrp-1)*0.15];

        scatter(curr_x, currMed, 5, 'MarkerFaceColor', cl, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4+(iGrp-1)*0.15)
    end
end

hold off


end