clear
clc
% close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/';
rr = dir(fullfile(parentDir, 'assemblyTuning_finalResults'));


currSessions = [9 10 11];
nSessions    = numel(currSessions);


epochNames = {'pre';'maze'; 'post'; 'remaze'};
% epochNames = {'pre'; 'run'; 'post'};


activeUnits = cell(nSessions, 1);
nUnits      = nan(nSessions, 1);

spatialTuning    = cell(nSessions, 1);
spatialTuning_re = cell(nSessions, 1);

learnedTuning_sub = cell(nSessions, 1);

learnedTuningPFcorr       = cell(nSessions, 1);
remaze_post_LTcorrelation = cell(nSessions, 1);
remaze_maze_LTcorrelation = cell(nSessions, 1);


gw = gausswindow(3,9); % for smoothing the tunings 
flag = 1;


for iSess = 1:nSessions
    

    sessionNumber = currSessions(iSess);
    sessionName   = rr(sessionNumber+2).name
    basePath      = fullfile(parentDir, 'assemblyTuning_finalResults', sessionName);
        
    load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
    spikes = spikes_pyr;

    if sessionNumber == 9
        load(fullfile(parentDir, 'sessions_calculated_PBEinfo4', sessionName, 'spikes', [sessionName '.spikes2.mat']), 'spikes_pyr2')%%%%
        spikes_re = spikes_pyr2;%%%%
    end

%     % measuring the firing rate during re-maze
%     load(fullfile(parentDir, 'Bapun_NSD_datasets', sessionName, [sessionName '.epochs.mat']))
%     
%     if exist('re_maze', 'var')
%         remaze = double(re_maze);
%     elseif exist('maze2', 'var')
%         remaze = double(maze2);
%     else
%         error('There is no remaze in this session!')
%     end 
% 
%     fr_remaze = nan(numel(spikes), 1);
%     for iUnit = 1:numel(spikes)
%         fr_remaze(iUnit) = numel(find(spikes(iUnit).time > remaze(1) & spikes(iUnit).time < remaze(2)));
%     end
%     
%     activeUnits_remaze = find(fr_remaze > 100);


    % % learned tunings and PF-matching scores


    % POST

    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTunings', 'assemblyTuningPFcorr', 'activeUnits');
   
    okUnits = s.activeUnits.post; 
%     okUnits = intersect(s.activeUnits.post, activeUnits_remaze);

%     okUnits = intersect(s.activeUnits.pre, s.activeUnits.post);
%     okUnits = intersect(okUnits, s.activeUnits.run);
    nUnits(iSess) = numel(okUnits);


%     for iEpoch = 1:numel(epochNames)
% 
%         currEpoch = epochNames{iEpoch};
%         LTs.(currEpoch)           = s.assemblyTunings.(currEpoch).data;
%         LTPFcorr.(currEpoch).data = s.assemblyTuningPFcorr.(currEpoch).data;
%         LTPFcorr.(currEpoch).ui   = s.assemblyTuningPFcorr.(currEpoch).ui;
%     end
        
    LTs.pre           = s.assemblyTunings.pre.data;
    LTPFcorr.pre.data = s.assemblyTuningPFcorr.pre.data;
    LTPFcorr.pre.ui   = []; %s.assemblyTuningPFcorr.post.ui;
    

    LTs.post           = s.assemblyTunings.post.data 
    LTPFcorr.post.data = s.assemblyTuningPFcorr.post.data;
    LTPFcorr.post.ui   = []; %s.assemblyTuningPFcorr.post.ui;

    clear s


%     s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Jan2022_first2hours.mat']), 'assemblyTunings', 'assemblyTuningPFcorr');
    
%     LTs.post           = s.assemblyTunings.post.data;
%     LTPFcorr.post.data = s.assemblyTuningPFcorr.post.data;
%     LTPFcorr.post.ui   = s.assemblyTuningPFcorr.post.ui;

%     clear s

    

%     s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores_Lthresh1e_3_radonUI.mat']), 'assemblyTunings_sub', 'asTuningPFcorr_sub');
    
%     LTs.post           = s.assemblyTunings_sub.post.rt_ui.data(:, :, 4);
%     LTPFcorr.post.data = s.asTuningPFcorr_sub.post.rt_ui.data(:, 4);
%     LTPFcorr.remaze.ui   = squeeze(s.asTuningPFcorr_sub.post.rt_ui.ui(:, 4, :));

    
%     clear s

    

    
    % Maze
    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_maze_stricterThetaDetection.mat']), 'learnedTunings', 'learnedTuningPFcorr');
         
    LTs.maze      = s.learnedTunings;
    LTPFcorr.maze = s.learnedTuningPFcorr;
  
    clear s



%     s = load(fullfile(parentDir, sessionName, 'assemblyTunings', [sessionName '.assemblyTunings_activeRun.mat']), 'learnedTunings', 'learnedTuningPFcorr');
%      
%     LTs.maze      = s.learnedTunings;
%     LTPFcorr.maze = s.learnedTuningPFcorr;
%     
%     for iUnit = 1:nUnits(iSess)
%         LTs.maze(iUnit, :) = conv(LTs.maze(iUnit, :), gw, 'same');
%     end
%     
%     clear s

    
    
    
    % remaze
    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_remaze_stricterThetaDetection.mat']), 'learnedTunings', 'learnedTuningPFcorr');
    
    LTs.remaze      = s.learnedTunings;
    LTPFcorr.remaze = s.learnedTuningPFcorr;
    
    clear s




    spikes = spikes(okUnits);

    if sessionNumber == 9
        spikes_re = spikes_re(okUnits);
    end

    
    % track spatial tunings
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    interpPosBins   = linspace(0, nPosBins, 200);
    nPosBins_interp = numel(interpPosBins);
    
    
    % non-directional spatial tuning
    spatialTunings_merge  = zeros(nUnits(iSess), nPosBins_interp);
    spatialTunings_merge_re = zeros(nUnits(iSess), nPosBins_interp); %%%

    peakPFlocation    = zeros(nUnits(iSess), 1);
    peakPFfiring      = zeros(nUnits(iSess), 1);
    
    for iUnit = 1:nUnits(iSess)
    
        currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;
        currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
        currSpatialTuning = currSpatialTuning./max(currSpatialTuning);
    
        spatialTunings_merge(iUnit, :) = currSpatialTuning;
        [peakPFfiring(iUnit), peakPFlocation(iUnit)] = max(currSpatialTuning);

        % place fields during REMAZE

        if sessionNumber == 9
            currSpatialTuning = spikes_re(iUnit).spatialTuning_smoothed.uni;
        else
            currSpatialTuning = spikes(iUnit).spatialTuning_smoothed_re.uni;
        end

        currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
        currSpatialTuning = currSpatialTuning./max(currSpatialTuning);

        spatialTunings_merge_re(iUnit, :) = currSpatialTuning;
    
    end
    

    [~, sortIdx] = sort(peakPFlocation, 'ascend');
    sortedUnitIdx = okUnits(sortIdx);
    
    temp = nanmean(spatialTunings_merge, 1); %#ok<NANMEAN> 
    
    if flag == 1
        startBin = find(~isnan(temp), 1, 'first');
        endBin   = find(~isnan(temp), 1, 'last');
        flag = 0;
    end
    
    nPosBins_interp = endBin - startBin + 1;

    spatialTuning{iSess} = spatialTunings_merge(sortIdx, startBin:endBin);
    spatialTuning_re{iSess} = spatialTunings_merge_re(sortIdx, startBin:endBin);
    
    for iEpoch = 1:4
        currEpoch = epochNames{iEpoch};
    
        % assembly Tunings
        learnedTuning_sub{iSess}.(currEpoch) = nan(nUnits(iSess), nPosBins_interp);
    
        for iUnit = 1:nUnits(iSess)
    
            currLearnedTuning = LTs.(currEpoch)(sortedUnitIdx(iUnit), :);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;

            if ismember(iEpoch , [1 3])
                currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            end

            currLearnedTuning = currLearnedTuning(startBin:endBin);
    
            learnedTuning_sub{iSess}.(currEpoch)(iUnit, :) = currLearnedTuning/max(currLearnedTuning);
        end
        
        if ismember(iEpoch , [1 3])

            allCorrs = corr(learnedTuning_sub{iSess}.(currEpoch)', spatialTuning{iSess}');
            learnedTuningPFcorr{iSess}.(currEpoch).data   = diag(allCorrs);
%             learnedTuningPFcorr{iSess}.(currEpoch).ui     = LTPFcorr.(currEpoch).ui(sortedUnitIdx, :);

        else
            learnedTuningPFcorr{iSess}.(currEpoch)   = LTPFcorr.(currEpoch)(sortedUnitIdx, :);
        end
        
    end
    

    % correlation between LTs in different period (not individually with the place field)

    allCorrs = corr(LTs.remaze', LTs.post');
    temp     = diag(allCorrs);
    
    remaze_post_LTcorrelation{iSess} = temp(sortedUnitIdx);
    
    

    allCorrs = corr(LTs.remaze', LTs.maze');
    temp     = diag(allCorrs);

    remaze_maze_LTcorrelation{iSess} = temp(sortedUnitIdx);

end



% sorting the units when all sessions considered together
cnctSpatialTunings = cell2mat(spatialTuning);
[~, PFpeakLocs]    = max(cnctSpatialTunings, [], 2);
[~, sortIdx]       = sort(PFpeakLocs, 'ascend');

cnctSpatialTunings_re = cell2mat(spatialTuning_re);
[~, PFpeakLocs_re]    = max(cnctSpatialTunings_re, [], 2);




PFpeakLocs = PFpeakLocs(sortIdx);


temp = mean(cnctSpatialTunings, 1); 

flag = 1;
if flag == 1
%     startBin = find(temp > 0, 1, 'first');
%     endBin   = find(temp > 0, 1, 'last');

    startBin = 1;
    endBin = nPosBins_interp;
    flag = 0;
end

nPosBins_interp = endBin - startBin + 1;
cnctSpatialTunings    = cnctSpatialTunings(:, startBin:endBin);
cnctSpatialTunings_re = cnctSpatialTunings_re(:, startBin:endBin);



for iEpoch = 1:4
    
    currEpoch = epochNames{iEpoch};
    
    currLearnedTunings = cell(nSessions, 1);
    currLTPFcorrs.data = cell(nSessions, 1);
%     currLTPFcorrs.ui   = cell(nSessions, 1);

    for iSess = 1:nSessions      
        currLearnedTunings{iSess} = learnedTuning_sub{iSess}.(currEpoch)(:, startBin:endBin);

        if ismember(iEpoch , [1 3])
            currLTPFcorrs.data{iSess} = learnedTuningPFcorr{iSess}.(currEpoch).data;
%             currLTPFcorrs.ui{iSess}   = learnedTuningPFcorr{iSess}.(currEpoch).ui;
        else
            currLTPFcorrs.data{iSess} = learnedTuningPFcorr{iSess}.(currEpoch);
        end
    end

    cnctLearnedTunings = cell2mat(currLearnedTunings);
    cnctLTPFcorrs.data = cell2mat(currLTPFcorrs.data);

    if ismember(iEpoch , [1 3])
%         cnctLTPFcorrs.ui = cell2mat(currLTPFcorrs.ui);
    end

    % sorting again because we are concatentating multiple sessions that
    % were sorted separately
    
    LTs_pooled.(currEpoch)       = cnctLearnedTunings(sortIdx, :);
    
    LTPFcorrs_pooled.(currEpoch).data = cnctLTPFcorrs.data(sortIdx);

    if ismember(iEpoch , [1 3])
%         LTPFcorrs_pooled.(currEpoch).ui   = cnctLTPFcorrs.ui(sortIdx, :);
    end

end


remaze_post_LTcorrelation = cell2mat(remaze_post_LTcorrelation);
remaze_post_LTcorrelation = remaze_post_LTcorrelation(sortIdx);

remaze_maze_LTcorrelation = cell2mat(remaze_maze_LTcorrelation);
remaze_maze_LTcorrelation = remaze_maze_LTcorrelation(sortIdx);



%% plot the LTs separately for low and high post LT-PF correlation

hiUnits = find(LTPFcorrs_pooled.post.data > nanmedian(LTPFcorrs_pooled.post.data));
loUnits = setdiff(1:numel(LTPFcorrs_pooled.post.data), hiUnits);




%% plot the learned tunings for maze theta, post PBE, and remaze theta 


plotheight = 200;
plotwidth  = 220;%180;
fontsize   = 6;


nor = 1;
noc = 5;

leftmargin = 25;  rightmargin = 25;    topmargin = 25;    bottommargin = 25;

gapc = 3;
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
noc_epoch = 1;

leftmargin_epoch = 0;  rightmargin_epoch = 0;     bottommargin_epoch = 0;    topmargin_epoch = 0;
gapr_epoch = 0;    gapc_epoch = 3;



% plot the place fields

ax1 = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

cnctSpatialTunings_1 = cnctSpatialTunings(sortIdx, :);
imagesc(cnctSpatialTunings_1) %(hiUnits, :)

tn_units = size(cnctSpatialTunings_1, 1);

colormap('jet')
xlim([0 nPosBins_interp])
ylim([1 tn_units])
% ylim([0.5 numel(hiUnits)+0.5])
set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [1 tn_units], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', 6)


tt = 2;

for iEpoch = [2:4]
    
    currEpoch = epochNames{iEpoch};

    tn_units = size(LTs_pooled.(currEpoch) , 1);
    
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{tt}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{tt}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    

    for ii=1
        
        ax1 = axes('position',sub_pos_epoch{ii},'XGrid','off','XMinorGrid','off','FontSize', fontsize,'Box','off','Layer','top');
    
        hold on

        imagesc(LTs_pooled.(currEpoch))  %(hiUnits, :)
%         plot(PFpeakLocs, 1:numel(PFpeakLocs), 'color', 'w', 'linewidth', 0.5)

        
        colormap('jet')
        xlim([0 nPosBins_interp])
        ylim([1 tn_units])
%         ylim([0.5 numel(hiUnits)+0.5])

        set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [1 tn_units], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', fontsize)


%         if iEpoch > 1
           set(ax1, 'yTick', [])
%         end

    end
    
    tt = tt+1;
    
end


% remaze place fields

ax5 = axes('position',sub_pos{5},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

cnctSpatialTunings_2 = cnctSpatialTunings_re(sortIdx, :);
imagesc(cnctSpatialTunings_2) %(hiUnits, :)

tn_units = size(cnctSpatialTunings_2, 1);

colormap('jet')
xlim([0 nPosBins_interp])
ylim([1 tn_units])
% ylim([0.5 numel(loUnits)+0.5])
set(ax5, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', 6)


%% To determine the dependence of PF_rePF correlation on Post LT-PF correlation


allCorrs = corr(cnctSpatialTunings_2', cnctSpatialTunings_1');
PFstability = diag(allCorrs);

plotheight = 100;
plotwidth  = 100;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on


scatter(LTPFcorrs_pooled.post.data, PFstability, 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

ylabel('PF stability', 'fontsize', fontsize)
xlabel('post LT-PF correlation', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(LTPFcorrs_pooled.post.data)) LTPFcorrs_pooled.post.data];
[b, ~,~,~, stats] = regress(PFstability, X);
pval = stats(3);
r2   = stats(1);

estimated_PFstability = b(1)+b(2)*LTPFcorrs_pooled.post.data;
h1 = plot(LTPFcorrs_pooled.post.data, estimated_PFstability, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(LTPFcorrs_pooled.post.data), nanmedian(estimated_PFstability), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 min(1, xl(2)+0.1)])
ylim([yl(1)- 0.1 min(1, yl(2)+0.1)])




%% PF correlation and KL divergence


% plot the distribution of LT-PFcorrelation

plotheight = 81;
plotwidth  = 101;
fontsize   = 6;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


customCDF(LTPFcorrs_pooled)


yl = ylim;

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'XTick', -1:0.5:1)

xlabel({'LT-PF'; 'Pearson correlation coeff.'}, 'fontsize', fontsize)

ylabel('fraction of units', 'fontsize', fontsize)



%% relationship between LTPF correlation in maze, post ad remaze


% % %  post and remaze


% [corrCoeff, pval, inclusionIdx] = nancorr(LTPFcorrs_pooled.post, LTPFcorrs_pooled.remaze, []);


plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(LTPFcorrs_pooled.post.data, LTPFcorrs_pooled.remaze.data, 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('post LT-PF correlation', 'fontsize', fontsize)
ylabel('remaze LT-PF correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(LTPFcorrs_pooled.post.data)) LTPFcorrs_pooled.post.data];
[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.remaze.data, X);
pval = stats(3);
r2   = stats(1);

estimated_remaze = b(1)+b(2)*LTPFcorrs_pooled.post.data;
h1 = plot(LTPFcorrs_pooled.post.data, estimated_remaze, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(LTPFcorrs_pooled.post.data), nanmedian(estimated_remaze), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 xl(2)+0.1])
ylim([yl(1)- 0.1 yl(2)+0.1])

grid on




% % % remaze and maze

%% [corrCoeff, pval, inclusionIdx] = nancorr(LTPFcorrs_pooled.remaze, LTPFcorrs_pooled.maze, []);

plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(LTPFcorrs_pooled.maze.data, LTPFcorrs_pooled.remaze.data, 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('maze LT-PF correlation', 'fontsize', fontsize)
ylabel('remaze LT-PF correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(LTPFcorrs_pooled.remaze.data)) LTPFcorrs_pooled.maze.data];
[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.remaze.data, X);
pval = stats(3);
r2   = stats(1);

estimated_remaze = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_remaze_maze = LTPFcorrs_pooled.remaze.data - estimated_remaze;


h1 = plot(LTPFcorrs_pooled.maze.data, estimated_remaze, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(LTPFcorrs_pooled.maze.data), nanmedian(estimated_remaze), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


axis square



%% multiple regression analysis

X = [LTPFcorrs_pooled.pre.data LTPFcorrs_pooled.maze.data LTPFcorrs_pooled.post.data];

lm = fitlm(X, PFstability);
% lm = fitlm(X, LTPFcorrs_pooled.remaze.data);



%% post residual

X = [ones(size(LTPFcorrs_pooled.maze.data)) LTPFcorrs_pooled.maze.data];
[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.post.data, X);
pval = stats(3);
r2   = stats(1);

estimated_post = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_post_maze = LTPFcorrs_pooled.post.data - estimated_post;



%% PF stability residual

X = [ones(size(LTPFcorrs_pooled.maze.data)) LTPFcorrs_pooled.maze.data];
[b, ~,~,~, stats] = regress(PFstability, X);
pval = stats(3);
r2   = stats(1);

estimated_PFstability = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_PFstability_maze = PFstability - estimated_PFstability;




%% to plot PF stability residual vs POST LT-PFcorrs residual both regressed against MAZE LT-PFcorrs
plotheight = 100;
plotwidth  = 100;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on


scatter(residual_post_maze, residual_PFstability_maze, 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

ylabel('residual PF stability', 'fontsize', fontsize)
xlabel('residual PRE LT-PF correlation', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(residual_post_maze)) residual_post_maze];
[b, ~,~,~, stats] = regress(residual_PFstability_maze, X);
pval = stats(3);
r2   = stats(1);

estimated_resPFstability = b(1)+b(2)*residual_post_maze;
h1 = plot(residual_post_maze, estimated_resPFstability, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(residual_post_maze), nanmedian(estimated_resPFstability), {sprintf('R2 = %.2f', r2); sprintf('p = %.1e', pval)}, 'fontsize', 4)


xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 min(1, xl(2)+0.1)])
ylim([yl(1)- 0.1 min(1, yl(2)+0.1)])



%% post-rePF correlation in different tirtiles of the residual PF stability


% allCorrs = corr(cnctSpatialTunings_2', LTs_pooled.post');
allCorrs = corr(LTs_pooled.remaze', LTs_pooled.post');

post_rePFcorr = diag(allCorrs);

xvariable = residual_remaze_maze;

quartiles = [0 33.3;33.3 66.7; 66.7 100];
% quartiles = [0 50;50 100];

nq = size(quartiles, 1);


med = nan(nq,1);
lq  = nan(nq,1);
hq  = nan(nq,1);
xCenter = nan(nq,1);

med_ui = nan(nq,1);
lq_ui  = nan(nq,1);
hq_ui  = nan(nq,1);
xCenter_ui = nan(nq,1);



currQuartileUnits = cell(nq,1);

for iq = 1:nq

    currQuartileUnits{iq} = find(xvariable > prctile(xvariable, quartiles(iq, 1)) & xvariable <= prctile(xvariable, quartiles(iq, 2)));
    
%     allCorrs = corr(cnctSpatialTunings_2(currQuartileUnits{iq}, :)', LTs_pooled.post(currQuartileUnits{iq}, :)');
    allCorrs = corr(LTs_pooled.remaze(currQuartileUnits{iq}, :)', LTs_pooled.post(currQuartileUnits{iq}, :)');


    % data    
    post_rePF_quartile.data{iq} = diag(allCorrs);

    med(iq) = nanmedian(post_rePF_quartile.data{iq});
    lq(iq) = prctile(post_rePF_quartile.data{iq}, 25);
    hq(iq) = prctile(post_rePF_quartile.data{iq}, 75);
    xCenter(iq) = nanmean(xvariable(currQuartileUnits{iq}));

    

    % shuffled LTs

    nonDiagonals = setdiff(allCorrs, diag(allCorrs));
    post_rePF_quartile.ltui{iq} = nonDiagonals;

    med_ui(iq) = nanmedian(post_rePF_quartile.ltui{iq});
    lq_ui(iq) = prctile(post_rePF_quartile.ltui{iq}, 25);
    hq_ui(iq) = prctile(post_rePF_quartile.ltui{iq}, 75);


end


plotheight = 100;
plotwidth  = 100;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



hold on
scatter(xvariable, post_rePFcorr, 2 , 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 

plot(xCenter, med, 'color', [1 0 0 0.5], 'linewidth', 1)
plot(xCenter+0.05, med_ui, 'color', [0 0 0 0.5], 'linewidth', 1)


pval = nan(nq, 1);
zval = nan(nq, 1);

for iq = 1:nq
    line([xCenter(iq) xCenter(iq)], [lq(iq) hq(iq)], 'color', [1 0 0 0.5], 'linewidth', 1)
    line([xCenter(iq) xCenter(iq)]+0.05, [lq_ui(iq) hq_ui(iq)], 'color', [0 0 0 0.5], 'linewidth', 1)

    [pval(iq), ~, stats] = ranksum(post_rePF_quartile.data{iq}, post_rePF_quartile.ltui{iq}, 'tail', 'right');
    zval(iq) = stats.zval;

    text(xCenter(iq), hq(iq)+0.1, sprintf('p=%.2e', pval(iq)), 'fontsize', 4)

end


ylabel('PRE LT-REMAZE PF correlation', 'fontsize', fontsize)
xlabel('residual PF stability', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])


axis square





%% temporary



plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


% unitsSet1  = find(remaze_post_LTcorrelation > 0.5 & temp < -0.5);
% unitsSet2  = find(remaze_post_LTcorrelation > 0.5 & temp > -0.5);

remaze_maze_LTPFCorrdiff = LTPFcorrs_pooled.remaze.data - LTPFcorrs_pooled.maze.data;

% unitsSet1  = find(residual_remaze_maze < -0.3 & residual_post_maze < -0.3);
% % unitsSet2  = find(residual_remaze_maze > 0.2  & residual_post_maze >  0.2);
% unitsSet2 = currUnits;

% unitsSet2  = unitsSet2(randperm(numel(unitsSet2), numel(unitsSet1)));

% unitsSet3  = find(residual_remaze_maze < -0.3  & residual_post_maze >  0);

% unitsSet3  = find(remaze_maze_LTPFCorrdiff < -0.5 & remaze_post_LTcorrelation > 0.5);



nUnits = numel(residual_remaze_maze);
otherUnits = 1:nUnits; %setdiff(1:nUnits, [unitsSet1; unitsSet2]);

hold on

% scatter(residual_post_maze(unitsSet1), residual_remaze_maze(unitsSet1), 5, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(residual_post_maze(unitsSet2), residual_remaze_maze(unitsSet2), 5, 'r', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(residual_post_maze(unitsSet3), residual_remaze_maze(unitsSet3), 5, 'b', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 


scatter(residual_post_maze(otherUnits), residual_remaze_maze(otherUnits), 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 



ylabel('residual LT-PF correlation_{(remaze, maze)}', 'fontsize', fontsize)
xlabel('residual LT-PF correlation_{(post, maze)}', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(residual_post_maze)) residual_post_maze];
[b, ~,~,~, stats] = regress(residual_remaze_maze, X); % residual
pval = stats(3);

r2  = stats(1);

residualCalc = b(1)+b(2)*residual_post_maze;
h1 = plot(residual_post_maze, residualCalc, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(residual_post_maze), nanmedian(residualCalc), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


axis square

%%

plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


% unitsSet1  = find(remaze_post_LTcorrelation > 0.5 & temp < -0.5);
% unitsSet2  = find(remaze_post_LTcorrelation > 0.5 & temp > -0.5);

remaze_maze_LTPFCorrdiff = LTPFcorrs_pooled.remaze.data - LTPFcorrs_pooled.maze.data;

% unitsSet1  = find(residual_remaze_maze < -0.3 & residual_post_maze < -0.3);
% unitsSet2  = find(residual_remaze_maze > 0.2  & residual_post_maze >  0.2);

unitsSet1  = find(remaze_maze_LTPFCorrdiff < -0.5 & remaze_post_LTcorrelation > 0.5);
unitsSet2  = find(remaze_maze_LTPFCorrdiff > -0.5 & remaze_post_LTcorrelation > 0.5);
unitsSet3  = find(remaze_maze_LTPFCorrdiff > -0.5 & remaze_post_LTcorrelation < 0.2);
unitsSet4  = find(remaze_maze_LTPFCorrdiff < -0.5 & remaze_post_LTcorrelation < 0.2);


nUnits = numel(residual_remaze_maze);
otherUnits = 1:nUnits; %setdiff(1:nUnits, [unitsSet1; unitsSet2; unitsSet3; unitsSet4]);

hold on

% scatter(remaze_post_LTcorrelation(unitsSet1), remaze_maze_LTPFCorrdiff(unitsSet1), 2, 'y', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(remaze_post_LTcorrelation(unitsSet2), remaze_maze_LTPFCorrdiff(unitsSet2), 2, [1 0 1], 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(remaze_post_LTcorrelation(unitsSet3), remaze_maze_LTPFCorrdiff(unitsSet3), 2, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(remaze_post_LTcorrelation(unitsSet4), remaze_maze_LTPFCorrdiff(unitsSet4), 2, 'b', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 


scatter(remaze_maze_LTPFCorrdiff(otherUnits), remaze_post_LTcorrelation(otherUnits), 2 , 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 



ylabel('remaze-post LT correlation', 'fontsize', fontsize)
xlabel('\Delta(remaze,maze) LT-PF correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

% X = [ones(size(remaze_post_LTcorrelation)) remaze_maze_LTPFCorrdiff];
% [b, ~,~,~, stats] = regress(remaze_post_LTcorrelation, X); % residual
% pval = stats(3);
% 
% r2  = stats(1);
% 
% residualCalc = b(1)+b(2)*remaze_maze_LTPFCorrdiff;
% h1 = plot(remaze_maze_LTPFCorrdiff, residualCalc, 'color', [0 0 0 0.5], 'linewidth', 1);
% text(max(remaze_maze_LTPFCorrdiff), nanmedian(residualCalc), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


axis square








%%
tt = 1;

currMed = nanmedian(residual_post_maze);


group1 = residual_remaze_maze(residual_post_maze > currMed);
group2 = residual_remaze_maze(residual_post_maze < currMed);

plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


[xData, yData] = myViolin_oneSided(group1, 0.01);
h1 = patch(tt+10*xData, yData, [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

lq  = prctile(group1, 25);
hq  = prctile(group1, 75);
med = nanmedian(group1);

line([tt tt+0.5], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
line([tt tt+0.5], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
line([tt tt+0.5], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 1)




[xData, yData] = myViolin_oneSided(group2, 0.01); 
h2 = patch(tt-10*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

lq  = prctile(group2, 25);
hq  = prctile(group2, 75);
med = nanmedian(group2);

line([tt-0.5 tt], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
line([tt-0.5 tt], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
line([tt-0.5 tt], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 1)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])



%%



% residual_post_ui = nan(100, 1);
% for i_ui = 1:100
%     residual_post_ui(i_ui)  = calr2(residual_remaze_maze, LTPFcorrs_pooled.post.ui(:, i_ui));
% end



%%




%% controling for the maze LTPF, what is the correlation between remaze LTPF and post LTPF?


plotheight = 80;
plotwidth  = 80;
fontsize   = 6;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);





remaze_maze_LTPFCorrdiff = LTPFcorrs_pooled.remaze.data - LTPFcorrs_pooled.maze.data;


unitsSet1  = find(remaze_maze_LTPFCorrdiff < -0.5 & LTPFcorrs_pooled.post.data < 0);
unitsSet2  = find(remaze_maze_LTPFCorrdiff > -0.4  & LTPFcorrs_pooled.post.data > 0.7);



nUnits = numel(residual_remaze_maze);
otherUnits = setdiff(1:nUnits, [unitsSet1; unitsSet2]);


hold on

scatter(LTPFcorrs_pooled.post.data(unitsSet1), remaze_maze_LTPFCorrdiff(unitsSet1), 5, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
scatter(LTPFcorrs_pooled.post.data(unitsSet2), remaze_maze_LTPFCorrdiff(unitsSet2), 5, 'r', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
scatter(LTPFcorrs_pooled.post.data(otherUnits), remaze_maze_LTPFCorrdiff(otherUnits), 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 




xlabel('post LT-PF correlation', 'fontsize', fontsize)
ylabel('\Delta_{(remaze,maze)} LT-PF correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(LTPFcorrs_pooled.post.data)) LTPFcorrs_pooled.post.data];
[b, ~,~,~, stats] = regress(remaze_maze_LTPFCorrdiff, X); % residual
pval = stats(3);
r2_resid_post  = stats(1);

residualCalc = b(1)+b(2)*LTPFcorrs_pooled.post.data;
h1 = plot(LTPFcorrs_pooled.post.data, residualCalc, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(LTPFcorrs_pooled.post.data), nanmedian(residualCalc), {sprintf('R2=%.2f', r2_resid_post); sprintf('p=%.1e', pval)}, 'fontsize', 4)


axis square




%% ui surrogates of post

r2_resid_pre = 0.0045;
r2_resid_run = 0.0011;


residual_post_ui = nan(100, 1);
for i_ui = 1:100
%     residual_post_ui(i_ui)  = calr2(remaze_maze_LTPFCorrdiff, LTPFcorrs_pooled.post.ui(:, i_ui));
end


plotheight = 80;
plotwidth  = 150;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

% variable.ui = residual_post_ui;

% boxplot(variable.ui, 'boxstyle', 'filled', 'colors', [0.3 0.3 0.3])

scatter(2, r2_resid_pre, 5, 'b', 'filled', 'Marker', 'square')
scatter(2, r2_resid_run, 5, 'g', 'filled', 'Marker', 'square')
scatter(2, r2_resid_post, 5, 'r', 'filled', 'Marker', 'square')


% customHistogram(variable)

axis tight

% yl = ylim;
% 
% line([r2_resid_post r2_resid_post], [0 yl(2)], 'linewidth', 1, 'color', 'r')
% line([r2_resid_pre r2_resid_pre], [0 yl(2)], 'linewidth', 1, 'color', 'b')
% line([r2_resid_run r2_resid_run], [0 yl(2)], 'linewidth', 1, 'color', 'g')

xlim([0 3])
ylim([0 0.06])

ylabel('R^2', 'fontsize', fontsize)
% ylabel('probability', 'fontsize', fontsize)

set(gca, 'ytick', 0:0.02:0.06, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])




%% how remaze LTs are correlated with corresponding post and maze LTs


plotheight = 80;
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(remaze_post_LTcorrelation, remaze_maze_LTcorrelation, 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('remaze and post LT correlation', 'fontsize', fontsize)
ylabel('remaze and maze LT correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

aa = remaze_post_LTcorrelation;
bb = remaze_maze_LTcorrelation;

minn = min([aa; bb]);
maxx = max([aa; bb]);

h1 = line([minn maxx], [minn maxx], 'color', [0 0 0 0.5], 'linewidth', 1, 'linestyle', ':', 'DisplayName', 'identity line');

axis square

set(gca, 'xtick', -1:0.5:1, 'ytick', -1:0.5:1)



%% plotting the place fields and Learend tunings for the units with high remaze-post LT correlation but low remaze-maze LT correlations


% currUnits = find(remaze_post_LTcorrelation > 0.5 & remaze_maze_LTcorrelation < 0.2);
% currUnits = find(remaze_post_LTcorrelation > 0.6 & remaze_maze_LTcorrelation > 0.6);

% currUnits = find(remaze_post_LTcorrelation > 0.5 & temp < -0.5);
% currUnits = find(remaze_post_LTcorrelation > 0.5 & remaze_maze_LTPFCorrdiff < -0.5);
% currUnits = find(LTPFcorrs_pooled.post.data > 0.8 & residual > 0.3);

% idx = sort(randperm(numel(unitsSet2), numel(unitsSet1)));

% currUnits = unitsSet2(idx);
% currUnits = unitsSet1;

plotheight = 100;
plotwidth  = 180;
fontsize   = 6;


nor = 1;
noc = 4;

leftmargin = 25;  rightmargin = 25;    topmargin = 25;    bottommargin = 25;

gapc = 3;
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
noc_epoch = 1;

leftmargin_epoch = 0;  rightmargin_epoch = 0;     bottommargin_epoch = 0;    topmargin_epoch = 0;
gapr_epoch = 0;    gapc_epoch = 3;



% plot the place fields

ax1 = axes('position',sub_pos{1},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

cnctSpatialTunings_1 = cnctSpatialTunings(sortIdx(currUnits), :);
imagesc(cnctSpatialTunings_1)

tn_units = size(cnctSpatialTunings_1, 1);

colormap('jet')
xlim([0 nPosBins_interp])
ylim([0.5 tn_units+0.5])
set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [1 tn_units], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', 6)


tt = 2;

for iEpoch = 1:3
    
    currEpoch = epochNames{iEpoch};

    tn_units = size(LTs_pooled.(currEpoch)(currUnits, :) , 1);
    
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); %% this gives us the positions for the subplots, each corresponded to an event

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{tt}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{tt}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    

    for ii=1
        
        ax1 = axes('position',sub_pos_epoch{ii},'XGrid','off','XMinorGrid','off','FontSize', fontsize,'Box','off','Layer','top');
    
        hold on

        imagesc(LTs_pooled.(currEpoch)(currUnits, :)) 
%         plot(PFpeakLocs, 1:numel(PFpeakLocs), 'color', 'w', 'linewidth', 0.5)

        
        colormap('jet')
        xlim([0 nPosBins_interp])
        ylim([1 tn_units])

        set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [1 tn_units], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', fontsize)


%         if iEpoch > 1
           set(ax1, 'yTick', [])
%         end

    end
    ylim([0.5 tn_units+0.5])


    tt = tt+1;
    
end

%%%




%% 


plotheight = 100;
plotwidth  = 100;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

unitsSet1  = find(remaze_post_LTcorrelation > 0.3670 & remaze_maze_LTcorrelation < 0.5271);
% unitsSet2  = find(remaze_post_LTcorrelation > 0.5 & temp > -0.5);


nUnits = numel(temp);
otherUnits = setdiff(1:nUnits, unitsSet1);

scatter(remaze_maze_LTcorrelation(unitsSet1), remaze_post_LTcorrelation(unitsSet1), 5, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
% scatter(remaze_post_LTcorrelation(unitsSet2), temp(unitsSet2), 5, 'r', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)


scatter(remaze_maze_LTcorrelation(otherUnits), remaze_post_LTcorrelation(otherUnits), 5, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)



xlabel('remaze and maze LT correlation', 'fontsize', fontsize)
ylabel('remaze and post LT correlation', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

% aa = remaze_post_LTcorrelation;
% bb = remaze_maze_LTcorrelation;
% 
% minn = min([aa; bb]);
% maxx = max([aa; bb]);
% 
% h1 = line([minn maxx], [minn maxx], 'color', [0 0 0 0.5], 'linewidth', 1, 'linestyle', ':', 'DisplayName', 'identity line');

axis square

set(gca, 'xtick', -1:0.5:1, 'ytick', -1:0.5:1)




%% plotting the regresssion R2 for different thresh on maze LT-PFcorrelation

% thresholds = [-inf 0.3 0.5 0.7 0.9];
% n = numel(thresholds);
% 
% remaze_post_LTPF_corr = nan(n, 1);
% post_maze_LTPF_corr   = nan(n, 1);
% remaze_maze_LTPF_corr = nan(n, 1);
% 
% 
% for iThresh = 1:n
%         
%     currThresh = thresholds(iThresh);
% 
%     mazehighPFcorrUnits = find(LTPFcorrs_pooled.maze >= currThresh);
% 
%     remaze_correspond_mazehigh = LTPFcorrs_pooled.remaze(mazehighPFcorrUnits);
%     post_correspond_mazehigh = LTPFcorrs_pooled.post(mazehighPFcorrUnits);
%     maze_correspond_mazehigh = LTPFcorrs_pooled.maze(mazehighPFcorrUnits);
% 
%    
%     remaze_post_LTPF_corr(iThresh)  = calr2(post_correspond_mazehigh, remaze_correspond_mazehigh);
%     post_maze_LTPF_corr(iThresh)    = calr2(post_correspond_mazehigh, maze_correspond_mazehigh);
%     remaze_maze_LTPF_corr(iThresh)  = calr2(remaze_correspond_mazehigh, maze_correspond_mazehigh);
% 
% end
% 
% plotheight = 150;
% plotwidth  = 120;
% fontsize   = 6;
% 
% f= figure;
% clf(f);
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% 
% hold on
% 
% plot(remaze_post_LTPF_corr, 'marker', 'square', 'color', [0 0 0 1], 'linewidth', 1, 'displayName', 'remaze post LTPFcorrs')
% plot(post_maze_LTPF_corr, 'marker', 'diamond', 'color', [0 0 0 0.7], 'linewidth', 1, 'displayName', 'post maze LTPFcorrs')
% plot(remaze_maze_LTPF_corr, 'marker', 'o', 'color', [0 0 0 0.4], 'linewidth', 1, 'displayName', 'remaze maze LTPFcorrs')
% 
% 
% xlabel('threshold on maze LT-PF correlation', 'fontsize', fontsize)
% ylabel('R2', 'fontsize', fontsize)
% 
% legend
% 
% set(gca, 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:n, 'xticklabel', {'no thresh.'; '0.3'; '0.5'; '0.7'; '0.9'}, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])
% xtickangle(gca, 45)





%% sub-functions


function customHistogram(variable)

epochNames = fieldnames(variable);

if isfield(variable.(epochNames{1}), 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:numel(epochNames)
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2

end 

allValues = cell(numel(epochNames), 1);
for iEpoch = 1:numel(epochNames)
    allValues{iEpoch} = variable.(epochNames{iEpoch});
end

allValues = cell2mat(allValues);

% allValues = allValues(allValues > prctile(allValues, 0.1) & allValues < prctile(allValues, 99.9));
% allValues(allValues > (nanmedian(allValues) + 2.5*(nanprctile(allValues, 75) - nanprctile(allValues, 25)))) = [];


bins = linspace(min(allValues), max(allValues)+0.01, 20);

% colors = {'b'; 'none';'r'};
% edgeColors = {'none'; 'k'; 'none'};
% medianColors = {'b'; 'k'; 'r'};

colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255];

hold on
for iEpoch = 1:(numel(epochNames))

    currEpoch = epochNames{iEpoch};
    pooledVar = variable.(currEpoch);     
    
    h = histc(pooledVar, bins); h(end) = [];
    
%     plot(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'color', [colors(iEpoch, :) 0.5], 'linewidth', 1, 'displayName', epochNames{iEpoch})
    area([0 bins(1:end-1)]+diff(bins(1:2))/2, [0; h]/sum(h), 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.2, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    
%     a(iEpoch) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/sum(h), 'Edgecolor', 'none', 'facecolor', colors(iEpoch, :), 'facealpha', 0.4);
%     hold on

%     medianCorr.(currEpoch) = nanmedian(pooledVar);

end

xlim([min(allValues), max(allValues)])

% rr = ylim;
% for iEpoch = 1:numel(epochNames)
%     currEpoch = epochNames{iEpoch};
%     line([medianCorr.(currEpoch) medianCorr.(currEpoch)], [0 rr(2)], 'color', [colors(iEpoch, :) 0.6], 'linewidth', 1)
% end


% significance scores
% ranksumPvalue = nan(3);
% for iEpoch = 1:3
%    var1 = variable.(epochNames{iEpoch});
%    
%    for jEpoch = setdiff(1:3, iEpoch)
%        var2 = variable.(epochNames{jEpoch});
%        
%        ranksumPvalue(iEpoch, jEpoch) = ranksum(var1, var2, 'tail', 'right');
%        
%    end
% end
       
% add the significance signs here ...  
    
    
end


function customCDF(variable)


epochNames = fieldnames(variable);

if isfield(variable.(epochNames{1}), 'data') % in case we are plotting the absolute values and not z-scores
    
    for iEpoch = 1:4
        variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
    end
    
    variable = variable2;
    clear variabale2
end 


colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255; [0 0 0]];

hold on

signedRank_pvalue = nan(4,1);
signedRank_zvalue = nan(4,1);

for iEpoch = 2:4

    currEpoch = epochNames{iEpoch};
    pooledVar = variable.(currEpoch);  

    [curr_p, ~, stats] = signrank(pooledVar, 0, 'tail', 'right');
    signedRank_pvalue(iEpoch) = curr_p;
    signedRank_zvalue(iEpoch) = stats.zval;

    
    [cdf_pooled, x_pooled] = ecdf(pooledVar);
    [x_pooled, idx] = unique(x_pooled);
    cdf_pooled = cdf_pooled(idx);
    
%     area(x_pooled, cdf_pooled, 'FaceColor', colors(iEpoch, :), 'FaceAlpha', 0.3, 'EdgeColor', colors(iEpoch, :), 'linewidth', 0.5, 'EdgeAlpha', 0.6)
    ax(iEpoch) = plot(x_pooled, cdf_pooled, 'color', colors(iEpoch, :), 'linewidth', 1, 'DisplayName', currEpoch)

    medianCorr.(currEpoch) = nanmedian(pooledVar);

end

ylim([0 1])
yl = ylim;
xl = xlim;


for iEpoch = 2:4
    currEpoch = epochNames{iEpoch};
    line([xl(1) medianCorr.(currEpoch)], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    line([medianCorr.(currEpoch) medianCorr.(currEpoch)], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors(iEpoch, :))
    text(medianCorr.(currEpoch), 0.07, sprintf('%.2f', medianCorr.(currEpoch)), 'fontsize', 6, 'color', colors(iEpoch, :))

end



% Friedman test

temp = [variable.maze variable.post variable.remaze];
temp = temp(~isnan(sum(temp, 2)), :);

[p_friedman, anovaTab, stats] = friedman(temp);



% significance scores



p_value = nan(4);
z_value  = nan(4);
for iEpoch = 2:4
   var1 = variable.(epochNames{iEpoch});
   
   for jEpoch = iEpoch+1:4
       var2 = variable.(epochNames{jEpoch});

       [p_value(iEpoch, jEpoch),~, stats] = signrank(var1, var2);
       z_value(iEpoch, jEpoch) = stats.zval;
   end
end
       
% add the significance signs here ...  
    

legend(ax, 'location', 'best')

    
end






function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end


function [corrCoeff, pval, inclusionIdx] = nancorr(var1, var2, var3)

if ~isempty(var3)
    inclusionIdx = ~isnan(var1) & ~isnan(var2) & ~isnan(var3);  
    [corrCoeff, pval] = partialcorr(var1(inclusionIdx), var2(inclusionIdx), var3(inclusionIdx));
else
    inclusionIdx = ~isnan(var1) & ~isnan(var2);
    [corrCoeff, pval] = corr(var1(inclusionIdx), var2(inclusionIdx));
end

end


function [r2, pval] = calr2(var1, var2)

X = [ones(size(var2)) var2];
[~, ~,~,~, stats] = regress(var1, X);
pval = stats(3);
r2   = stats(1);

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