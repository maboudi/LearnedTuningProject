clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 12 15:17 6]; % the recording sessions
nSessions    = numel(currSessions);

epochNames = {'pre'; 'run'; 'post'};
curr_rsm = 'wc_ts'; % the sequence metric and surrogate distribution used for calculating the replay scores



activeUnits = cell(nSessions, 1);
nUnits      = nan(nSessions, 1);

spatialTuning = cell(nSessions, 1);


learnedTuning_sub = cell(nSessions, 1);


learnedTuningPFcorr_sub         = cell(nSessions, 1);
learnedTuningPFcorr_sub_pscore  = cell(nSessions, 1);

learnedTuningPFKLdiv_sub        = cell(nSessions, 1);
learnedTuningPFKLdiv_sub_pscore = cell(nSessions, 1);



gw = gausswindow(5,15); % Gaussian Kernel to smooth the spatial tunnig curves
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

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores_Lthresh1e_3_wcts_lowMaxJump.mat']), ...
    'assemblyTunings_sub', 'asTuningPFcorr_sub', 'asTuningPFKLdiv_sub', ...
    'asTuning_prefPos_sub', 'asTuningPF_prefPosDist_sub');


replayScoreMethods = fieldnames(s.asTuningPFcorr_sub.pre);



spikes = spikes(okUnits);

% track spatial tunings
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

interpPosBins   = linspace(0, nPosBins, 200); % to use the same number of bins across sessions with different length of the track, hence number of position bins
nPosBins_interp = numel(interpPosBins);


% non-directional spatial tunings
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


temp = sum(spatialTunings_merge, 1);
if flag == 1
    startBin = find(~isnan(temp), 1, 'first');
    endBin   = find(~isnan(temp), 1, 'last');
    flag = 0;
end

spatialTuning{iSess} = spatialTunings_merge(sortIdx, startBin:endBin);

nPosBins_interp = endBin - startBin + 1;
nShuffles       = 10000;
        
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

        
        shuffle_PFcorr = nan(nUnits(iSess), nShuffles);
        allCorrs = corr(learnedTuning_sub{iSess}.(currEpoch)(:, :, iSub)', spatialTuning{iSess}'); 

        for i_s = 1:nShuffles % surrogate distributions by shuffling matching the PFs of non-identical units       
            cidx = randi(nUnits(iSess), nUnits(iSess), 1);
            shuffle_PFcorr(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:nUnits(iSess))', cidx));
        end
        learnedTuningPFcorr_sub{iSess}.(currEpoch).ui(:, :, iSub) = shuffle_PFcorr;
        learnedTuningPFcorr_sub{iSess}.(currEpoch).data(:, iSub)  = diag(allCorrs);
    end
    

    allPFcorr = calPrctile(s.asTuningPFcorr_sub.(currEpoch).(curr_rsm));
    learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch) = allPFcorr(sortedUnitIdx, :);

    
end


end


%% plot the learning tunings for different quartiles of replay score

plotheight = 300;
plotwidth  = 250;

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

plotwidth_epoch  = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
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
    

    pooledLT.(currEpoch) = pooledLT.(currEpoch)(sortIdx, :, :);
    
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

plotwidth_epoch  = (plotwidth - (rightmargin+leftmargin) - (noc-1) * gapc)/noc;
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
    
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, ...
        bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{nEpochs}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{nEpochs}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    

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

        title(sprintf('%s-Q%d', currEpoch, iSub), 'fontsize', fontsize, 'FontWeight', 'normal')
        
        axis square
        ax1.XGrid = 'on';
   
    end

    nEpochs = nEpochs + 1;
end


%% plot the PF fidleity scores for learned tunings calculated using
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



for iEpoch = 1:3
    
    
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
        
        PFcorr_pscore{iSess} = learnedTuningPFcorr_sub_pscore{iSess}.(currEpoch);
         
    end
    
    PFcorr_data  = cell2mat(PFcorr.data);
    PFcorr_ui    = cell2mat(PFcorr.ui); 

    
    PFcorr_pscore_pooled = cell2mat(PFcorr_pscore);
    
    
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
    

    currAx = gca;
    currAx.YGrid = 'on';
    
    if iEpoch >= 2
        currAx.YAxis.Visible = 'off';
    end
    

    currentVariable.(currEpoch) = currentVariable.(currEpoch)(~isnan(sum(currentVariable.(currEpoch), 2)), :);
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
    
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


