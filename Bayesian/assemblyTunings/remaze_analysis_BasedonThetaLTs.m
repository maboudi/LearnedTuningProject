clear
clc
% close all

parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/';
rr = dir(fullfile(parentDir, 'assemblyTuning_finalResults'));


currSessions = 8:11;
nSessions    = numel(currSessions);


epochNames = {'pre';'maze'; 'post'; 'post_late'; 'remaze'; 'maze_theta'}; 
dataTypes  = {'data'; 'ui'};



activeUnits = cell(nSessions, 1);
nUnits      = nan(nSessions, 1);
unitSessNum = cell(nSessions, 1);

spatialTuning    = cell(nSessions, 1);
spatialTuning_re = cell(nSessions, 1);

spatialInfo        = cell(nSessions, 1);
POST_LTstabilities = cell(nSessions, 1);
LTspatialInfo      = cell(nSessions, 1);

learnedTuning_sub = cell(nSessions, 1);

% learnedTuningPFcorr       = cell(nSessions, 1);
remaze_post_LTcorrelation = cell(nSessions, 1);
remaze_maze_LTcorrelation = cell(nSessions, 1);


gw   = gausswindow(3,9); % for smoothing the tunings 
flag = 1;


for iSess = 1:nSessions
    

    sessionNumber = currSessions(iSess);
    sessionName   = rr(sessionNumber+2).name
    basePath      = fullfile(parentDir, 'assemblyTuning_finalResults', sessionName);
       


    % % PFs

    load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr')
    spikes = spikes_pyr;

    if sessionNumber == 9
        load(fullfile(parentDir, 'sessions_calculated_PBEinfo4', sessionName, 'spikes', [sessionName '.spikes2.mat']), 'spikes_pyr2')%%%%
        spikes_re = spikes_pyr2;%%%%
    end




    % % sleep/quiet wake LTs

    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), 'assemblyTunings', 'assemblyTuningPFcorr', 'activeUnits'); % 
   
    stableUnits = s.activeUnits.post; 

    nUnits(iSess) = numel(stableUnits);

    LTs.pre = s.assemblyTunings.pre.data; %%%
    learnedTuningPFcorr.pre.data{iSess} = s.assemblyTuningPFcorr.pre.data; %%%




    %%% added temporary to include MAZE PBEs instead of MAZE theta
    LTs.maze = s.assemblyTunings.run.data; 
    learnedTuningPFcorr.maze.data{iSess} = s.assemblyTuningPFcorr.run.data;
    

    
    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_Lthresh1e_3.mat']), 'assemblyTunings', 'assemblyTuningPFcorr', 'assemblyTuningSpatialInfo'); % _after6hoursPost

    LTs.post = s.assemblyTunings.post.data;
    learnedTuningPFcorr.post.data{iSess} = s.assemblyTuningPFcorr.post.data;
    
%     LTspatialInfo{iSess} = s.assemblyTuningSpatialInfo.post.data; 
%     LTspatialInfo{iSess} = sum(LTs.post.^2, 2)./(sum(LTs.post, 2).^2);
   
%     LTspatialInfo{iSess} = nan(size(LTs.post, 1), 1);
%     for iUnit = 1:size(LTs.post, 1)
%         LTspatialInfo{iSess}(iUnit) = ginicoeff(LTs.post(iUnit, :));
%     end

    

    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs_last4hoursPOST.mat']), 'assemblyTunings', 'assemblyTuningPFcorr'); % _1st_2hoursPOST

    LTs.post_late = s.assemblyTunings.post.data;
    learnedTuningPFcorr.post_late.data{iSess} = s.assemblyTuningPFcorr.post.data;
    



    % % MAZE theta LTs
    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_activeRun_binDur0.02.mat']), 'learnedTunings'); %'.assemblyTunings_maze_stricterThetaDetection.mat'
         
    LTs.maze_theta = s.learnedTunings;
%     learnedTuningPFcorr.maze_theta.data{iSess} = s.learnedTuningPFcorr;


    
    
    % % REMAZE theta LTs
    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_remaze_stricterThetaDetection.mat']), 'learnedTunings', 'learnedTuningPFcorr');
    
    LTs.remaze      = s.learnedTunings;
    learnedTuningPFcorr.remaze.data{iSess} = s.learnedTuningPFcorr;
    

    

    % % load stabilities

    s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.POSTstability.mat']), 'POSTstabilities');
    POSTstability = s.POSTstabilities; 




    % % 

    spikes = spikes(stableUnits);

    if sessionNumber == 9 % rat U
        spikes_re = spikes_re(stableUnits);
    end

    
    % track spatial tunings
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    interpPosBins   = linspace(0, nPosBins, 200);
    nPosBins_interp = numel(interpPosBins);
    
    
    % non-directional spatial tuning
    spatialTunings_merge    = zeros(nUnits(iSess), nPosBins_interp);
    spatialTunings_merge_re = zeros(nUnits(iSess), nPosBins_interp); %%
    spatialInfo_sess        = nan(nUnits(iSess), 1);

    peakPFlocation    = zeros(nUnits(iSess), 1);
    peakPFfiring      = zeros(nUnits(iSess), 1);
    
    for iUnit = 1:nUnits(iSess)
    
        currSpatialTuning = spikes(iUnit).spatialTuning_smoothed.uni;
        currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);
%         currSpatialTuning = currSpatialTuning./max(currSpatialTuning);
        currSpatialTuning = (currSpatialTuning - nanmean(currSpatialTuning))/nanstd(currSpatialTuning);

        % added later
        
%         spatialInfo_sess(iUnit) = spikes(iUnit).peakFR.uni;
        temp = spikes(iUnit).peakPosBin.uni;

        if temp == 1
            spatialInfo_sess(iUnit) = nan;
        else
            spatialInfo_sess(iUnit) = min(temp - 0, nPosBins - temp)/nPosBins;
        end
    
        spatialTunings_merge(iUnit, :) = currSpatialTuning;
        [peakPFfiring(iUnit), peakPFlocation(iUnit)] = max(currSpatialTuning);
        


        % place fields during REMAZE (all session except RatS, for which the position during remaze isn't available)
        
        if ismember(sessionNumber, 9:11) % RatU_Day2, RatV_Day1, RatV_day3

            if sessionNumber == 9
                currSpatialTuning = spikes_re(iUnit).spatialTuning_smoothed.uni;
            else
                currSpatialTuning = spikes(iUnit).spatialTuning_smoothed_re.uni;
            end
    
            currSpatialTuning = interp1(1:nPosBins, currSpatialTuning, interpPosBins);

            currSpatialTuning = (currSpatialTuning - nanmean(currSpatialTuning))/nanstd(currSpatialTuning);

%             currSpatialTuning = currSpatialTuning./max(currSpatialTuning);
    
            spatialTunings_merge_re(iUnit, :) = currSpatialTuning;
        end
    
    end


    [~, sortIdx] = sort(peakPFlocation, 'ascend');
    sortedUnitIdx = stableUnits(sortIdx);
    
    temp = nanmean(spatialTunings_merge, 1); %#ok<NANMEAN> 
    
    if flag == 1
        startBin = find(~isnan(temp), 1, 'first');
        endBin   = find(~isnan(temp), 1, 'last');
        flag = 0;
    end
    
    nPosBins_interp = endBin - startBin + 1;

    spatialTuning{iSess}    = spatialTunings_merge(sortIdx, startBin:endBin);
    spatialTuning_re{iSess} = spatialTunings_merge_re(sortIdx, startBin:endBin);

    % added later

    spatialInfo{iSess}        = spatialInfo_sess(sortIdx);
    POST_LTstabilities{iSess} = POSTstability(sortIdx);
%     LTspatialInfo{iSess}      = LTspatialInfo{iSess}(sortIdx); 
    

    %
    
    for iEpoch = 1:numel(epochNames)
        currEpoch = epochNames{iEpoch};
    
        % assembly Tunings
        learnedTuning_sub{iSess}.(currEpoch) = nan(nUnits(iSess), nPosBins_interp);
    
        for iUnit = 1:nUnits(iSess)
    
            currLearnedTuning = LTs.(currEpoch)(sortedUnitIdx(iUnit), :);
            currLearnedTuning = interp1(1:nPosBins, currLearnedTuning, interpPosBins);
            currLearnedTuning(isnan(currLearnedTuning)) = 0;

            if ismember(iEpoch , 1:4) % ripple LTs
                currLearnedTuning = conv(currLearnedTuning, gw, 'same');
            end

            currLearnedTuning = currLearnedTuning(startBin:endBin);
    
%             learnedTuning_sub{iSess}.(currEpoch)(iUnit, :) = currLearnedTuning; 
            % learnedTuning_sub{iSess}.(currEpoch)(iUnit, :) = currLearnedTuning/max(currLearnedTuning);
            learnedTuning_sub{iSess}.(currEpoch)(iUnit, :) = (currLearnedTuning - nanmean(currLearnedTuning))/nanstd(currLearnedTuning);
        end
        
        allCorrs = corr(learnedTuning_sub{iSess}.(currEpoch)', spatialTuning{iSess}');

%         if ismember(iEpoch, [2 4])
            learnedTuningPFcorr.(currEpoch).data{iSess} = diag(allCorrs);
%         end

        nShuffles      = 10000;
        shuffle_PFcorr = nan(nUnits(iSess), nShuffles);
        for i_s = 1:nShuffles
            
            cidx = randi(nUnits(iSess), nUnits(iSess), 1);
            shuffle_PFcorr(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:nUnits(iSess))', cidx));
        end
        learnedTuningPFcorr.(currEpoch).ui{iSess} = shuffle_PFcorr;   
    end
    

    % correlation between LTs in different period (not individually with the place field)

    allCorrs = corr(LTs.remaze', LTs.post');
    temp     = diag(allCorrs);
    remaze_post_LTcorrelation{iSess} = temp(sortedUnitIdx);

    unitSessNum{iSess} = [sessionNumber * ones(numel(sortedUnitIdx), 1) sortedUnitIdx];
  
end



% sorting the units again in pooled data

cnctSpatialTunings = cell2mat(spatialTuning);
[~, PFpeakLocs]    = max(cnctSpatialTunings, [], 2);
[~, sortIdx]       = sort(PFpeakLocs, 'ascend');

cnctSpatialTunings = cnctSpatialTunings(sortIdx, :);
PFpeakLocs         = PFpeakLocs(sortIdx);



if ismember(sessionNumber, 9:11) % except RatS
    cnctSpatialTunings_re  = cell2mat(spatialTuning_re);
    cnctSpatialTunings_re  = cnctSpatialTunings_re(sortIdx, :);
end

cnctSpatialInfo = cell2mat(spatialInfo);
cnctSpatialInfo = cnctSpatialInfo(sortIdx);

cnctPOST_LTstabilities = cell2mat(POST_LTstabilities);
cnctPOST_LTstabilities = cnctPOST_LTstabilities(sortIdx);


% cnctLTspatialInfo = cell2mat(LTspatialInfo);
% cnctLTspatialInfo = cnctLTspatialInfo(sortIdx);


cnctUnitSessNum = cell2mat(unitSessNum);
cnctUnitSessNum = cnctUnitSessNum(sortIdx, :);


for iEpoch = 1:numel(epochNames)
    
    currEpoch = epochNames{iEpoch};
    
    currLearnedTunings = cell(nSessions, 1);

    for dt = 1:2    
        currLTPFcorrs.(dataTypes{dt}) = cell(nSessions, 1);
    end

    for iSess = 1:nSessions      
        currLearnedTunings{iSess} = learnedTuning_sub{iSess}.(currEpoch);
        

        for dt = 1:2
            currLTPFcorrs.(dataTypes{dt}){iSess} = learnedTuningPFcorr.(currEpoch).(dataTypes{dt}){iSess};
        end
    end

    cnctLearnedTunings = cell2mat(currLearnedTunings);
    LTs_pooled.(currEpoch) = cnctLearnedTunings(sortIdx, :);

    for dt = 1:2
        cnctLTPFcorrs = cell2mat(currLTPFcorrs.(dataTypes{dt}));
        LTPFcorrs_pooled.(currEpoch).(dataTypes{dt}) = cnctLTPFcorrs(sortIdx);
    end

    % calculate the Gini coefficient of the LTs

    cnctLTspatialInfo.(currEpoch) = nan(size(LTs_pooled.(currEpoch), 1), 1);
    for iUnit = 1:size(LTs_pooled.(currEpoch), 1)
        cnctLTspatialInfo.(currEpoch)(iUnit) = ginicoeff(LTs_pooled.(currEpoch)(iUnit, :));
    end

end

% LT similarity between POST and REMAZE
remaze_post_LTcorrelation = cell2mat(remaze_post_LTcorrelation);
remaze_post_LTcorrelation = remaze_post_LTcorrelation(sortIdx);




%% Main Figure 4


%% plot the LTs separately for low and high POST LT-PF correlation


allCorrs         = corr(cnctSpatialTunings_re', cnctSpatialTunings');
PFstability.data = diag(allCorrs);

nShuffles      = 10000;
currNunits     = size(allCorrs, 1);

shuffle_corr = nan(nUnits(iSess), nShuffles);
for i_s = 1:nShuffles  
    cidx = randi(nUnits(iSess), nUnits(iSess), 1);
    shuffle_corr(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:nUnits(iSess))', cidx));
end
PFstability.ui = shuffle_corr; 


% testing the similarity between the reMAZE PFs and MAZE LTs

allCorrs  = corr(cnctSpatialTunings_re', LTs_pooled.maze_theta');
reMAZE_PF_MAZE_theta_corr = diag(allCorrs);



plotheight = 90;
plotwidth  = 90;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(PFstability.data, reMAZE_PF_MAZE_theta_corr, 1, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5)


xlabel('r(reMAZE PF, MAZE PF)', 'fontsize', fontsize)
ylabel('r(reMAZE PF, MAZE LT)', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 min(1, xl(2)+0.1)])
ylim([yl(1)- 0.1 min(1, yl(2)+0.1)])


line(xl, xl, 'LineStyle', ':', 'linewidth', 1, 'color', [0.7 0 0 0.5])




% okUnits = find(~isnan(PFstability));

thresh   = nanmedian(LTPFcorrs_pooled.post.data); % 


hiUnits = find(LTPFcorrs_pooled.post.data > thresh); %find(LTPFcorrs_pooled.post > thresh);
loUnits = setdiff(1:numel(LTPFcorrs_pooled.post.data), hiUnits);




%% Figure 4a: plot the learned tunings for maze theta, post PBE, and remaze theta 
% 
% selectUnits3 = selectUnits2(setdiff(1:numel(selectUnits2), [1 2 5 6 16 17
% 22 27 28 30 31])); % these aer related to the residual analysis
% selectUnits4 = selectUnits1(setdiff(1:numel(selectUnits1), [2 4 8 10 13 22 23 26]));


% selectUnits3 = selectUnits1([8 9 11 14 17 20 21 24 25 26 27 28 31 33 42 44]);
% selectUnits4 = selectUnits2([2 5 7 10 11 12 14 17 19 23 25 28 31 33 34 36]);

selectUnits = [14 35 42 53 68 84 93 116 157 149 128 168 173 118 154 15];
selectUnits_sortIDs = nan(16,1);
for ii = 1:numel(selectUnits)
    selectUnits_sortIDs(ii) = sortIdx(sortedUnitIdx == selectUnits(ii));
end


units2plot  = selectUnits_sortIDs; %hiUnits; %% determine here the units that we want to plot their learned tunings
% units2plot  = selectUnits4;

nUnits2plot = numel(units2plot);

plotheight = 100;
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


imagesc(cnctSpatialTunings(units2plot, :)) 


colormap('jet')
xlim([0 nPosBins_interp])
ylim([0.5 nUnits2plot+0.5])
set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [1 nUnits2plot], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', 6)



tt = 4;
for iEpoch = 3% 2:4 % maze LTs, post LTs, remaze LTs
        
    sub_pos_epoch = subplot_pos(plotwidth_epoch, plotheight_epoch, leftmargin_epoch, rightmargin_epoch, bottommargin_epoch, topmargin_epoch, noc_epoch, nor_epoch, gapr_epoch, gapc_epoch); 

    for ii = 1 : nor_epoch
        for jj = 1 : noc_epoch

            sub_pos_epoch{ii,jj}(1) = (sub_pos_epoch{ii,jj}(1) * plotwidth_epoch + sub_pos{tt}(1) * plotwidth) / plotwidth;
            sub_pos_epoch{ii,jj}(2) = (sub_pos_epoch{ii,jj}(2) * plotheight_epoch + sub_pos{tt}(2) * plotheight) / plotheight;
            sub_pos_epoch{ii,jj}(3) = sub_pos_epoch{ii,jj}(3) * plotwidth_epoch / plotwidth;
            sub_pos_epoch{ii,jj}(4) = sub_pos_epoch{ii,jj}(4) * plotheight_epoch / plotheight;
        end
    end
    

    for ii = 1
        
        ax1 = axes('position',sub_pos_epoch{ii},'XGrid','off','XMinorGrid','off','FontSize', fontsize,'Box','off','Layer','top');
    
        hold on

        imagesc(LTs_pooled.(epochNames{iEpoch})(units2plot, :)) 

        
        colormap('jet')
        xlim([0 nPosBins_interp])
        ylim([0.5 nUnits2plot+0.5])

        set(ax1, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', fontsize)


    end
    
    tt = tt+1;
    
end


% remaze place fields

ax5 = axes('position',sub_pos{5},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

imagesc(cnctSpatialTunings_re(units2plot, :)) 

colormap('jet')
xlim([0 nPosBins_interp])
ylim([0.5 nUnits2plot+0.5])
set(ax5, 'xTick', [0 nPosBins_interp], 'XTickLabel', {'0'; '1'}, 'yTick', [], 'YDir', 'normal', 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'fontsize', 6)




%% Figure 4b, top inset: PF correlation and KL divergence


% plot the distribution of LT-PFcorrelation

plotheight = 60;
plotwidth  = 130;
fontsize   = 6;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


customCDF(LTPFcorrs_pooled)


% adding the ReMAZE place fields


med_data = nanmedian(PFstability.data);
med_ui   = nanmedian(PFstability.ui);

signedRank_pvalue = numel(find(med_ui >= med_data))/size(med_ui, 2);


% hold on
% [curr_p, ~, stats] = signrank(PFstability, 0, 'tail', 'right');
% signedRank_pvalue_rePF = curr_p;
% signedRank_zvalue_rePF = stats.zval;

    
[cdf_pooled, x_pooled] = ecdf(PFstability.data);
[x_pooled, idx] = unique(x_pooled);
cdf_pooled = cdf_pooled(idx);

ax(iEpoch) = plot(x_pooled, cdf_pooled, 'color', 'c', 'linewidth', 1, 'DisplayName', currEpoch);

medianCorr_re = nanmedian(PFstability.data);


ylim([0 1])
yl = ylim;
xl = xlim;


line([xl(1) medianCorr_re], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', 'c')
line([medianCorr_re medianCorr_re], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', 'c')
text(medianCorr_re, 0.07, sprintf('%.2f', medianCorr_re), 'fontsize', 6, 'color', 'c')



set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])
set(gca, 'XTick', -1:0.5:1)

xlabel('LT-PF Pearson correlation coeff.', 'fontsize', fontsize)

ylabel('fraction of units', 'fontsize', fontsize)



%% Distribution of corr(reMAZE PF, POST LT) after regressing out corr(POST LT, MAZE PF) 

allCorrs = corr(cnctSpatialTunings', LTs_pooled.post');
postLT_mazePF_corr = diag(allCorrs); % this is the regressor for both postLT_remazePF_corr and remazePF_mazePF_corr


allCorrs = corr(cnctSpatialTunings_re', LTs_pooled.post');
postLT_remazePF_corr = diag(allCorrs);



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

scatter(res_postLT_remazePF_corr, res_remazePF_mazePF_corr, 2, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5)

xlabel('res-postLT-remazePF-corr', 'fontsize', fontsize)
ylabel('res-remazeP-mazePF-corr', 'fontsize', fontsize)



X = [ones(size(postLT_remazePF_corr)) postLT_mazePF_corr];
[b, ~,~,~, stats] = regress(postLT_remazePF_corr, X);
pval = stats(3);
r2   = stats(1);

estimated_postLT_remazePF_corr = b(1)+b(2)*postLT_mazePF_corr;

res_postLT_remazePF_corr = postLT_remazePF_corr - estimated_postLT_remazePF_corr;




%% Distribution of corr(reMAZE PF, MAZE PF) after regressing out corr(POST LT, MAZE PF) 

allCorrs = corr(cnctSpatialTunings', cnctSpatialTunings_re');
remazePF_mazePF_corr = diag(allCorrs);

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

scatter(res_postLT_remazePF_corr, res_remazePF_mazePF_corr, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('res-postLT-remazePF-corr', 'fontsize', fontsize)
ylabel('res-remazeP-mazePF-corr', 'fontsize', fontsize)


X = [ones(size(remazePF_mazePF_corr)) postLT_mazePF_corr];
[b, ~,~,~, stats] = regress(remazePF_mazePF_corr, X);
pval = stats(3);
r2   = stats(1);

estimated_remazePF_mazePF_corr = b(1)+b(2)*postLT_mazePF_corr;

res_remazePF_mazePF_corr = remazePF_mazePF_corr - estimated_remazePF_mazePF_corr;



%% plotting the residual reMAZE PF - MAZE PF correlation versus residual reMAZE PF - POST LT correlation
% Regressor for variables is POST LT - MAZE PF correlation


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

scatter(res_postLT_remazePF_corr, res_remazePF_mazePF_corr, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('res-postLT-remazePF-corr', 'fontsize', fontsize)
ylabel('res-remazeP-mazePF-corr', 'fontsize', fontsize)


X = [ones(size(res_postLT_remazePF_corr)) res_postLT_remazePF_corr];
[b, ~,~,~, stats] = regress(res_remazePF_mazePF_corr, X);
pval = stats(3);
r2   = stats(1);

estimated_res_remazePF_mazePF_corr = b(1)+b(2)*res_postLT_remazePF_corr;
res_res_remazePF_mazePF_corr = res_remazePF_mazePF_corr - estimated_res_remazePF_mazePF_corr;


h1 = plot(res_postLT_remazePF_corr, estimated_res_remazePF_mazePF_corr, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(res_postLT_remazePF_corr), nanmedian(estimated_res_remazePF_mazePF_corr), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)

xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 xl(2)])
ylim([yl(1)- 0.1 yl(2)]) 

axis square




%% Figure 4b, bottom insets: relationship between LTPF correlation in MAZE, POST and REMAZE

% % Figure 4b, bottom right inset: POST and REMAZE

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

scatter(LTPFcorrs_pooled.post.data, LTPFcorrs_pooled.remaze.data, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

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

xlim([xl(1)- 0.1 xl(2)])
ylim([yl(1)- 0.1 yl(2)])

axis square



% % Figure 4b, bottom left inset: REMAZE and MAZE

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

scatter(LTPFcorrs_pooled.maze.data, LTPFcorrs_pooled.remaze.data, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)

xlabel('maze LT-PF correlation', 'fontsize', fontsize)
ylabel('remaze LT-PF correlation', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(LTPFcorrs_pooled.remaze.data)) LTPFcorrs_pooled.maze.data];
[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.remaze.data, X);
pval = stats(3);
r2   = stats(1);LTs_pooled.post(:)

estimated_remaze = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_remaze_maze = LTPFcorrs_pooled.remaze.data - estimated_remaze;


h1 = plot(LTPFcorrs_pooled.maze.data, estimated_remaze, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(LTPFcorrs_pooled.maze.data), nanmedian(estimated_remaze), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)

xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 xl(2)])
ylim([yl(1)- 0.1 yl(2)]) 

axis square



%% Multiple regression analysis (partial correlation of reMAZE PFs or LTs with MAZE PFs and LTs in preceding epochs)


DV = cnctSpatialTunings_re; % dependent variable 
% DV = LTs_pooled.remaze; % dependent variable 


average_Y  = mean([LTs_pooled.pre ; LTs_pooled.post], 'omitnan');
average_Y  = (average_Y - mean(average_Y, 'omitnan'))/std(average_Y, 'omitnan');
background = repmat(average_Y, [size(DV, 1) 1]);


predictors       = [background(:) LTs_pooled.pre(:) cnctSpatialTunings(:) LTs_pooled.maze_theta(:) LTs_pooled.maze(:)  LTs_pooled.post(:)]; 
predictor_labels = {
    'intercept';...
    'average LT'; ...
    '\color[rgb]{0.2431 0.3882 0.6824}PRE rplLTs';...
    '\color[rgb]{0.1490 0.1490 0.1490}MAZE PFs';...
    '\color[rgb]{0.2039 0.3843 0.1882}MAZE thtLTs';...
    '\color[rgb]{0.0627 0.5451 0.2667}MAZE rplLTs';...
    '\color[rgb]{0.8627 0.2039 0.3647}POST rplLTs'...
    }; % the intercep was included here too for simplicity 
% 



Npredictors      = numel(predictor_labels);   


lm = fitlm(predictors, DV(:));
    
R2     = lm.Rsquared.Ordinary;
coeffs = lm.Coefficients.Estimate;


% comapring the R^2 and regression coefficients to shuffles

nShuffles = 10000;

R2_shuffle     = nan(nShuffles, 1);
coeffs_shuffle = nan(Npredictors, nShuffles);


for ii = 1:nShuffles

    DV_randomized = DV(randperm(size(DV, 1)), :);
    
    shuffle_lm = fitlm(predictors, DV_randomized(:));

    R2_shuffle(ii)        = shuffle_lm.Rsquared.Ordinary;
    coeffs_shuffle(:, ii) = shuffle_lm.Coefficients.Estimate;
    
end

p_r2 = numel(find(R2 < R2_shuffle)) / nShuffles;


coeff_cmpr2shuffle = repmat(coeffs, [1, nShuffles]) < coeffs_shuffle;
p_coeffs = sum(coeff_cmpr2shuffle, 2) / nShuffles;




%% plot


plotheight = 150;
plotwidth  = 40*Npredictors;
fontsize   = 8;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

h = bar(diag(coeffs), 'stacked');
for ih = 1:2; h(ih); h(ih).FaceColor = 'none'; end 
h(3).FaceColor = [0.2431 0.3882 0.6824];
h(4).FaceColor = [0.1490 0.149  0.1490];   
h(5).FaceColor = [0.2039 0.3843 0.1882]; 
h(6).FaceColor = [0.0627 0.5451 0.2667];   
h(7).FaceColor = [0.8627 0.2039 0.3647]; 

for ih = 1:2; h(ih).EdgeColor = 'k'; h(ih).EdgeAlpha = 0.8; h(ih).LineWidth = 1; end 
for ih = 3:7; h(ih).EdgeColor = 'none'; h(ih).FaceAlpha = .8; end


for ip = 1:Npredictors
    text(ip, coeffs(ip)+0.03, calAnnotPval(p_coeffs(ip)), 'fontsize', 5, 'HorizontalAlignment', 'center')
end

text(0.75, 0.2, {sprintf('R^2=%.3f', R2) ;calAnnotPval(p_r2)}, 'fontsize', 5, 'HorizontalAlignment', 'center')


set(gca, 'xtick', 1:Npredictors, 'xticklabel', predictor_labels, 'fontsize', fontsize, 'YGrid', 'on')

ylabel({'reMAZE thtLTs'; 'regression coeffs'}, 'fontsize', 8)

ylim([-0.1 0.3])



%% POST LTs noisiness (measured by staibility across time, gini coefficient, sparsity, etc) vs their correlation with reMAZE PFs (or theta LTs)


reMAZE_tuning = cnctSpatialTunings_re; % the tuning could be PF or LT
% reMAZE_Tuning = LTs_pooled.remaze;


% xVar = cnctPOST_LTstabilities; % x variable
xVar = cnctLTspatialInfo.post;  


% correlation bw post LTs and reMAZE PF (or theta LTs)

allCorrs = corr(reMAZE_tuning', LTs_pooled.post'); 
POST_reMAZE_tuning_corr = diag(allCorrs); 


idx = ~isnan(POST_reMAZE_tuning_corr) & ~isnan(xVar);
[rho, pval] = corr(POST_reMAZE_tuning_corr(idx), xVar(idx));



plotheight = 120;
plotwidth  = 120;
fontsize   = 8;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

scatter(xVar, POST_reMAZE_tuning_corr, 2, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5);

xlim('padded')
ylim('padded')

xlabel('\color[rgb]{0.8627 0.2000 0.3647}POST \color[rgb]{0 0 0}rplLTs Gini coeff.', 'fontsize', fontsize)
ylabel('r(\color[rgb]{0.8627 0.2000 0.3647}POST \color[rgb]{0 0 0}rplLTs, reMAZE PFs)', 'fontsize', fontsize)

text(min(xVar), max(POST_reMAZE_tuning_corr), {sprintf('r = %.2f', rho); calAnnotPval(pval)}, 'fontsize', 6)



%% how differneces in Gini coefficient between POST and MAZE theta LTs is related to their differneces in correlation with the reMAZE PFs (or theta LTs)


reMAZE_tuning = cnctSpatialTunings_re; % the tuning could be PF or LT

xVar = cnctLTspatialInfo.post - cnctLTspatialInfo.maze_theta;


% Caluclate to what extent the POST LTs predict reMAZE beyond MAZE theta
% LTs

% just take the difference (not sure if this makes sense)

% allCorrs = corr(reMAZE_tuning', LTs_pooled.post'); 
% POST_reMAZE_tuning_corr = diag(allCorrs); 
% 
% allCorrs = corr(reMAZE_tuning', LTs_pooled.maze_theta'); 
% MAZE_reMAZE_tuning_corr = diag(allCorrs); 
% 
% yVar = POST_reMAZE_tuning_corr - MAZE_reMAZE_tuning_corr;


% Another way might be to calculate a partial correlation between POST LT
% and reMAZE tuning controlling for MAZE theta LT (or ripple LT)


curr_n_units = size(LTs_pooled.post, 1);
yVar = nan(curr_n_units, 1);

for iUnit = 1:curr_n_units
    yVar(iUnit) = partialcorr(reMAZE_tuning(iUnit, :)', LTs_pooled.post(iUnit, :)', LTs_pooled.maze(iUnit, :)'); % _theta
end


idx = ~isnan(xVar) & ~isnan(yVar);
[rho, pval] = corr(xVar(idx), yVar(idx));


plotheight = 150;
plotwidth  = 150;
fontsize   = 8;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


scatter(xVar, yVar, 2, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5);

xlim('padded')
ylim('padded')

xlabel({'Gini coeff.\color[rgb]{0.8627 0.2000 0.3647}POST\color[rgb]{0 0 0}rplLT'; '-Gini coeff. MAZE thtLT'}, 'fontsize', fontsize)
ylabel({'r(\color[rgb]{0.8627 0.2000 0.3647}POST\color[rgb]{0 0 0}rplLT, reMAZE PF)'; '-r(MAZE thtLT, reMAZE PF)'}, 'fontsize', fontsize)

text(min(xVar), max(yVar), {sprintf('r = %.2f', rho); calAnnotPval(pval)}, 'fontsize', 6)




%%


%%


% a0_vector = sum(repmat(lm.Coefficients.Estimate(5:end), [1, nPosBins_interp]) .* intercept1');

rr = lm.Fitted;
rr = reshape(rr, size(A2));
rr = (rr - repmat(nanmean(rr , 2), [1 nPosBins_interp]))./repmat(nanstd(rr , [], 2), [1 nPosBins_interp]);

figure; 
subplot(151)
imagesc(A2); colormap('jet'); set(gca, 'YDir', 'normal'); caxis([-2 2]); title('remaze PFs')
hold on
plot(PFpeakLocs, 1:size(A2,1), 'linewidth', 1, 'color', 'y')

subplot(152)
imagesc(rr); colormap('jet'); set(gca, 'YDir', 'normal'); caxis([-2 2]); title('regression fit')
hold on
plot(PFpeakLocs, 1:size(A2,1), 'linewidth', 1, 'color', 'y')

subplot(153)
imagesc(cnctSpatialTunings); colormap('jet'); set(gca, 'YDir', 'normal'); caxis([-2 2]); title('MAZE PF')
hold on
plot(PFpeakLocs, 1:size(A2,1), 'linewidth', 1, 'color', 'y')

subplot(154)
imagesc(LTs_pooled.post); colormap('jet'); set(gca, 'YDir', 'normal'); caxis([-2 2]); title('POST LT')
hold on
plot(PFpeakLocs, 1:size(A2,1), 'linewidth', 1, 'color', 'y')

subplot(155)
imagesc(LTs_pooled.pre); colormap('jet'); set(gca, 'YDir', 'normal'); caxis([-2 2]); title('PRE LT')
hold on
plot(PFpeakLocs, 1:size(A2,1), 'linewidth', 1, 'color', 'y')





% doing the shuffle to see if the coefficients for each variable stand out
a0 = nan(nPosBins_interp, 1);
a1 = nan(nPosBins_interp, 1);
a2 = nan(nPosBins_interp, 1);

for ibin= 1:nPosBins_interp

    X = [learendTuning_pooled.PBEs.pre(:, ibin) cnctSpatialTunings_1(:, ibin) learendTuning_pooled.PBEs.pre(:, ibin).*cnctSpatialTunings_1(:, ibin)]; %

    lm = fitlm(X, learendTuning_pooled.PBEs.post(:, ibin));
    a0(ibin) = lm.Coefficients.Estimate(1);
    
    a1(ibin) = lm.Coefficients.Estimate(2);

    a2(ibin) = lm.Coefficients.Estimate(3);

end


% % % test

X = [LTs_pooled.pre LTPFcorrs_pooled.maze.data LTPFcorrs_pooled.post.data];

% lm = fitlm(X, PFstability);
lm = fitlm(X, LTPFcorrs_pooled.remaze.data);


%% multiple regression analysis for individual neurons

% A2 = cnctSpatialTunings_re; 
A2 = LTs_pooled.remaze;

average_Y = nanmean([LTs_pooled.pre ; LTs_pooled.post]);
average_Y = (average_Y - nanmean(average_Y))/nanstd(average_Y);

background = repmat(average_Y, [size(A2, 1) 1]);


nCells     = size(A2, 1);
bg_indiv   = nan(nCells, 1);
pre_indiv  = nan(nCells, 1);
pf_indiv   = nan(nCells, 1);
post_indiv = nan(nCells, 1);

for cell = 1:nCells
    X = [background(cell, :)' LTs_pooled.pre(cell, :)'  LTs_pooled.maze_theta(cell, :)'  LTs_pooled.post(cell, :)']; % cnctSpatialTunings 
    lm = fitlm(X, A2(cell, :));
    
    currCoeffs = lm.Coefficients.Estimate;
    bg_indiv(cell)   = currCoeffs(2);
    pre_indiv(cell)  = currCoeffs(3);
    pf_indiv(cell)   = currCoeffs(4);
    post_indiv(cell) = currCoeffs(5);
end



%%
cellNum = 162;
figure; 
hold on; plot(LTs_pooled.pre(cellNum, :), 'color', 'b', 'displayname', 'PRE LT', 'linewidth', 3)
plot(LTs_pooled.maze_theta(cellNum, :), 'color', 'k', 'displayname', 'MAZE PF', 'linewidth', 3) % cnctSpatialTunings
hold on; plot(A2(cellNum, :), 'color', [0.5 1 0], 'DisplayName', 'reMAZE PF', 'linewidth', 3)
hold on; plot(LTs_pooled.post(cellNum, :), 'color', 'r', 'displayname', 'POST LT', 'linewidth', 3)

xlabel('Maze position')
ylabel('Tuning')
set(gca, 'fontsize', 12, 'linewidth', 2, 'box', 'off')

legend('box', 'off')

%%

figure;
scatter(pre_indiv, post_indiv, '.', 'DisplayName', 'Unit')
set(gca, 'fontsize', 12, 'linewidth', 2, 'box', 'off')
xlabel('MAZE PF \beta')
ylabel('POST LT \beta')
legend('box', 'off')



%% POST LT vs maze PF, maze theta LT, maze replay LT, PRE LT

A2 = LTs_pooled.post;

average_Y = nanmean([LTs_pooled.pre ; LTs_pooled.post]);
average_Y = (average_Y - nanmean(average_Y))/nanstd(average_Y);

background = repmat(average_Y, [size(A2, 1) 1]);


nCells     = size(A2, 1);
bg_indiv   = nan(nCells, 1);
pre_indiv  = nan(nCells, 1);
pf_indiv   = nan(nCells, 1);
mazeTht_indiv = nan(nCells, 1);
maze_indiv = nan(nCells, 1);

for cell = 1:nCells
    X = [background(cell, :)' LTs_pooled.pre(cell, :)'  cnctSpatialTunings(cell, :)' LTs_pooled.maze_theta(cell, :)'  LTs_pooled.maze(cell, :)']; % cnctSpatialTunings 
    lm = fitlm(X, A2(cell, :));
    
    currCoeffs = lm.Coefficients.Estimate;
    bg_indiv(cell)   = currCoeffs(2);
    pre_indiv(cell)  = currCoeffs(3);
    pf_indiv(cell)   = currCoeffs(4);
    mazeTht_indiv(cell) = currCoeffs(5);
    maze_indiv(cell) = currCoeffs(6);
end



%%
figure;
scatter(pre_indiv, maze_indiv, '.', 'DisplayName', 'Unit')
set(gca, 'fontsize', 12, 'linewidth', 2, 'box', 'off')
xlabel('MAZE PF \beta')
ylabel('POST LT \beta')
legend('box', 'off')


%% plot the distribution of r(POST LT, MAZE PF) and r(POST LT, reMAZE PF)

curr_nUnits = size(cnctSpatialTunings, 1);

nShuffles      = 10000;
shuffle_corrs  = nan(curr_nUnits, nShuffles);
allCorrs       = corr(cnctSpatialTunings', LTs_pooled.pre');

POSTLT_MAZEPF_corr.data = diag(allCorrs);

for i_s = 1:nShuffles
    i_s
    cidx = randi(curr_nUnits, curr_nUnits, 1);
    shuffle_corrs(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:curr_nUnits)', cidx));
end

POSTLT_MAZEPF_corr.ui = shuffle_corrs;



plotheight = 80;    
plotwidth  = 80;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


Pval_POSTLT_MAZEPF = customCDF1(POSTLT_MAZEPF_corr); % POSTLT_MAZEPF_corr

xlabel('r(MAZE PF, reMAZE PF)', 'fontsize', fontsize)
ylabel('fraction of units', 'fontsize', fontsize)

axis square

yl = ylim;
ylim([yl(1)-0.1 yl(2)])

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01], 'YTick', [0 0.5 1], 'XTick', -1:0.5:1, 'YGrid', 'on')





%% plot r(reMAZE PF, POST LT) vs. r(reMAZE PF, MAZE PF)

allCorr = corr(cnctSpatialTunings_re', cnctSpatialTunings');
reMAZE_MAZE_PFcorr = diag(allCorr);


allCorr = corr(cnctSpatialTunings', LTs_pooled.post');
MAZE_PF_POST_LTcorr = diag(allCorr);



plotheight = 150;
plotwidth  = 150;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(MAZE_PF_POST_LTcorr, reMAZE_MAZE_PFcorr, 2, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5)





%% Figure 4C: relationship between REMAZE residual and POST residual


% % % post residual


X = [ones(size(LTPFcorrs_pooled.maze.data)) LTPFcorrs_pooled.maze.data];

[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.post.data, X);
pval = stats(3);
r2   = stats(1);

estimated_post = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_post_maze = LTPFcorrs_pooled.post.data - estimated_post;



% % % plot - the blue and yellow units should be determined from the next steps

plotheight = 120;
plotwidth  = 120;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



nUnits = numel(residual_remaze_maze);
otherUnits = 1:nUnits; 

hold on


% 
% 
% unitsSet1 = find(residual_post_maze < -0.2 & residual_remaze_maze < -0.2);
% unitsSet2 = find(residual_post_maze > 0.2 & residual_remaze_maze > 0.2);
% 
% randgen = randperm(numel(unitsSet1));
% randgen = sort(randgen(1:20));
% unitsSet1 = unitsSet1(randgen);
% 
% randgen = randperm(numel(unitsSet1));
% randgen = sort(randgen(1:20));
% unitsSet2 = unitsSet2(randgen);

% otherUnits = setdiff(otherUnits, [unitsSet1; unitsSet3]);



% scatter(residual_post_maze(unitsSet1), residual_remaze_maze(unitsSet1), 3, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% scatter(residual_post_maze(unitsSet3), residual_remaze_maze(unitsSet3), 3, 'r', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 

scatter(residual_post_maze(otherUnits), residual_remaze_maze(otherUnits), 3, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 



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


%% plot the correlation values for units blue and yellow

% [allCorrs, allCorr_Pval] = corr(LTs_pooled.post', LTs_pooled.remaze','tail', 'right');
[allCorrs, allCorr_Pval] = corr(LTs_pooled.post', cnctSpatialTunings_re','tail', 'right');



remaze_post_LTcorrelations = diag(allCorrs);
remaze_post_LTcorrPvals    = diag(allCorr_Pval);



plotheight = 100;
plotwidth  = 100;

cl = redblue;


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


subplot(1,2,1)
imagesc(remaze_post_LTcorrelations(selectUnits3))
caxis([-1 1])

set(gca, 'YTick', [], 'XTick', [], 'YDir', 'normal')


subplot(1,2,2)
imagesc(remaze_post_LTcorrelations(selectUnits4))
caxis([-1 1])

colormap(cl)

set(gca, 'YTick', [], 'XTick', [], 'YDir', 'normal')





%% Figure 4e: relationship between REMAZE residual and PRE residual


% % % pre residual

X = [ones(size(LTPFcorrs_pooled.maze.data)) LTPFcorrs_pooled.maze.data];
[b, ~,~,~, stats] = regress(LTPFcorrs_pooled.pre.data, X);
pval = stats(3);
r2   = stats(1);

estimated_pre = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
residual_pre_maze = LTPFcorrs_pooled.pre.data - estimated_pre;



% plot 

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

scatter(residual_pre_maze, residual_remaze_maze, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 



ylabel('residual LT-PF correlation_{(REMAZE, MAZE)}', 'fontsize', fontsize)
xlabel('residual LT-PF correlation_{(PRE, MAZE)}', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(residual_pre_maze)) residual_pre_maze];
[b, ~,~,~, stats] = regress(residual_remaze_maze, X); % residual
pval = stats(3);

r2  = stats(1);

residualCalc = b(1)+b(2)*residual_pre_maze;
plot(residual_pre_maze, residualCalc, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(residual_pre_maze), nanmedian(residualCalc), {sprintf('R2=%.2f', r2); sprintf('p=%.1e', pval)}, 'fontsize', 4)


axis square




%% Figure 4d: post-remaze LT correlation in different tertiles of the PF stability


currEpoch = 'post';

remazeVar = cnctSpatialTunings_re;
xvariable = reMAZE_MAZE_PFcorr;%PFstability.data; 

% remazeVar = LTs_pooled.remaze;
% xvariable = residual_remaze_maze;

% subs = [0 33.3;33.3 66.7; 66.7 100];
subs = [0 50; 50 100];

nSub = size(subs, 1);


allCorrs = corr(remazeVar', LTs_pooled.(currEpoch)');
currEpoch_remazeCorr = diag(allCorrs);


med        = nan(nSub,1);
lq         = nan(nSub,1);
hq         = nan(nSub,1);
% xCenter    = nan(nSub,1);
xBeg       = nan(nSub, 1);
xEnd       = nan(nSub, 1);

% med_ui     = nan(nSub,1);
% lq_ui      = nan(nSub,1);
% hq_ui      = nan(nSub,1);
% xCenter_ui = nan(nSub,1);
p_value      = nan(nSub, 1);

nCurrUnits     = nan(nSub, 1);
currUnitSubset = cell(nSub, 1);
shuffle_corrs  = cell(nSub, 1);

for iSub = 1:nSub

    currUnitSubset{iSub} = find(xvariable > prctile(xvariable, subs(iSub, 1)) & xvariable <= prctile(xvariable, subs(iSub, 2)));
    nCurrUnits(iSub) = numel(currUnitSubset{iSub});
    
    allCorrs = corr(remazeVar(currUnitSubset{iSub}, :)', LTs_pooled.(currEpoch)(currUnitSubset{iSub}, :)');


%     for iUnit = 1:size(allCorrs, 1)
%         
%         otherUnitsSessNum = cnctUnitSessNum(1:size(allCorrs, 2));
% 
%         allCorrs(iUnit, otherUnitsSessNum ~= cnctUnitSessNum(iUnit)) = nan;
% 
%     end


    % data

    post_remazeCorr_sub.data{iSub} = diag(allCorrs);

    med(iSub) = nanmedian(post_remazeCorr_sub.data{iSub});
    lq(iSub)  = prctile(post_remazeCorr_sub.data{iSub}, 25);
    hq(iSub)  = prctile(post_remazeCorr_sub.data{iSub}, 75);

%     xCenter(iSub) = nanmean(xvariable(currUnitSubset{iSub}));
    xBeg(iSub) = min(xvariable(currUnitSubset{iSub}));
    xEnd(iSub) = max(xvariable(currUnitSubset{iSub}));
    
    

    nShuffles = 10000;
    med_ui    = nan(nShuffles, 1);
    shuffle_corrs{iSub} = nan(nCurrUnits(iSub), nShuffles);

    for i_s = 1:nShuffles
        cidx = randi(nCurrUnits(iSub), nCurrUnits(iSub), 1);
        shuffle_corrs{iSub}(:, i_s) = allCorrs(sub2ind(size(allCorrs), (1:nCurrUnits(iSub))', cidx)); % 

        med_ui(i_s) = nanmedian(shuffle_corrs{iSub}(:, i_s));
    end

    p_value(iSub) = numel(find(med_ui >= med(iSub)))/nShuffles;



%     % shuffled LTs: resulted from matching the LTs of non-oidentical units
% 
%     nonDiagonals = setdiff(allCorrs, diag(allCorrs));
% 
%     currEpoch_remazeCorr_sub.ltui{iSub} = nonDiagonals;
% 
%     med_ui(iSub) = nanmedian(currEpoch_remazeCorr_sub.ltui{iSub});
%     lq_ui(iSub) = prctile(currEpoch_remazeCorr_sub.ltui{iSub}, 25);
%     hq_ui(iSub) = prctile(currEpoch_remazeCorr_sub.ltui{iSub}, 75);


end





% % % test (the scatter plot of reMAZE - POST correlatio versus reMAZE-MAZE correlation)

%     
% 
% shuffleIdx = randperm(size(LTs_pooled.post, 1));
% 
% allCorrs = corr(cnctSpatialTunings', cnctSpatialTunings_re');
% remazePF_mazePF_corr = diag(allCorrs);
% 
% allCorrs = corr(LTs_pooled.post(shuffleIdx, :)', cnctSpatialTunings_re');
% remazePF_postLT_corr = diag(allCorrs);
% 
% plotScatter_plus_diffHistogram(remazePF_mazePF_corr, remazePF_postLT_corr);
% xlabel('r(reMAZE PF, MAZE PF)', 'fontsize', fontsize)
% ylabel('r(reMAZE PF, POST LT)', 'fontsize', fontsize)



plotheight = 80;
plotwidth  = 80;
fontsize   = 8;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on

scatter(remazePF_mazePF_corr, remazePF_postLT_corr, 2, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 


xlabel('r(reMAZE PF, MAZE PF)', 'fontsize', fontsize)
ylabel('r(reMAZE PF, POST LT)', 'fontsize', fontsize)

set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

axis square
xl = xlim;
yl = ylim;

xl = [xl(1)-0.2 xl(2)+0.2];
yl = [yl(1)-0.2 yl(2)+0.2];

xlim(xl)
ylim(yl)

minn = nanmin(remazePF_mazePF_corr); 
maxx = nanmax(remazePF_mazePF_corr);
line([minn maxx], [minn maxx], 'linewidth', 1, 'linestyle', ':', 'color', [[212 0 53]/255 0.9])

% 
% selectedUnits1 = find(remazePF_mazePF_corr < 0.1495 & remazePF_postLT_corr > 0.4);
% selectedUnits2 = find(residual_remaze_maze < prctile(residual_remaze_maze, 33.3) & currEpoch_remazeCorr > 0.3);



%% alternative plot 1

plotheight = 140;
plotwidth  = 120;
fontsize   = 6;
    
f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

for iSub = 1:nSub
    

    currData = post_remazeCorr_sub.data{iSub};
    [xData, yData] = myViolin_oneSided(currData, 0.01);
    h1 = patch(iSub+15*xData, yData, [212 0 53]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    lq  = prctile(currData, 25);
    hq  = prctile(currData, 75);
    med = median(currData);
    
    line([iSub iSub+0.5], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([iSub iSub+0.5], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
    line([iSub iSub+0.5], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 1)
    
    

    currSurrogate = shuffle_corrs{iSub}(:, 2);
    [xData, yData] = myViolin_oneSided(currSurrogate, 0.01); 
    h2 = patch(iSub-15*xData, yData, [1 1 1], 'EdgeColor', [212 0 53]/255, 'FaceAlpha', 0.4);
    
    lq  = prctile(currSurrogate, 25);
    hq  = prctile(currSurrogate, 75);
    med = median(currSurrogate);
    
    line([iSub-0.5 iSub], [lq lq], 'linestyle', '-', 'color', [212 0 53]/255, 'linewidth', 0.5)
    line([iSub-0.5 iSub], [hq hq], 'linestyle', '-', 'color', [212 0 53]/255, 'linewidth', 0.5)
    line([iSub-0.5 iSub], [med med], 'linestyle', '-', 'color', [212 0 53]/255, 'linewidth', 1)

end

ylim([-0.65 1])

xticks([1 2 3])

xticklabels({sprintf('%.2f - %.2f', xBeg(1), xEnd(1)); ...
    sprintf('%.2f - %.2f', xBeg(2), xEnd(2)); ...
    sprintf('%.2f - %.2f', xBeg(3), xEnd(3))})
xtickangle(45)           

% % alternative plot 2
% 
% plotheight = 140;
% plotwidth  = 100;
% fontsize   = 6;
%     
% f= figure;
% clf(f);
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% 
% 
% 
% hold on
% 
% currData = nan(numel(post_remazeCorr_sub.data{nSub}), nSub);
% 
% for ii = 1:nSub
%     currData(1:size(post_remazeCorr_sub.data{ii}), ii) = currEpoch_remazeCorr_sub.data{ii};
% end
% 
% h = violinplot(currData);
% 
% for ii = 1:nSub
% 
%     h(ii).ViolinColor = 'none';
%     h(ii).ViolinAlpha = 0.7;
%     h(ii).EdgeColor   = [212 0 53]/255;
%     h(ii).ShowData    = 0;
%     h(ii).MedianPlot.SizeData = 10;
% %         h(ii).MedianPlot.Marker = 'square';
%     h(ii).MedianPlot.MarkerEdgeColor = 'none';
%     h(ii).MedianPlot.MarkerFaceColor = [212 0 53]/255;       
%     h(ii).WhiskerPlot.LineWidth = 1;
%     h(ii).WhiskerPlot.Color = [212 0 53]/255;
% 
% end
% 
% 
% xticklabels({'1st tertile'; '2nd tertile'; '3rd tertile'})
% 


% % old plot
% 
% plotheight = 100;
% plotwidth  = 100;
% fontsize   = 6;
% 
% f= figure;
% clf(f);
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% 
% 
% 
% hold on
% scatter(xvariable, currEpoch_remazeCorr, 2 , 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)% rresidual 
% 
% plot(xCenter, med, 'color', [1 0 0 0.5], 'linewidth', 1)
% plot(xCenter+0.05, med_ui, 'color', [0 0 0 0.5], 'linewidth', 1)
% 
% 
% pval = nan(nSub, 1);
% zval = nan(nSub, 1);
% 
% for iSub = 1:nSub
%     line([xCenter(iSub) xCenter(iSub)], [lq(iSub) hq(iSub)], 'color', [1 0 0 0.5], 'linewidth', 1)
%     line([xCenter(iSub) xCenter(iSub)]+0.05, [lq_ui(iSub) hq_ui(iSub)], 'color', [0 0 0 0.5], 'linewidth', 1)
% 
%     [pval(iSub), ~, stats] = ranksum(currEpoch_remazeCorr_sub.data{iSub}, currEpoch_remazeCorr_sub.ltui{iSub}, 'tail', 'right');
%     zval(iSub) = stats.zval;
% 
%     text(xCenter(iSub), hq(iSub)+0.1, sprintf('p=%.2e', pval(iSub)), 'fontsize', 4)
% 
% end
% 
% 
% xlabel('residual LT-PF correlation_{(REMAZE, MAZE)}', 'fontsize', fontsize)
% ylabel(['REMAZE-' upper(currEpoch) 'LT correlation'], 'fontsize', fontsize)
% 
% set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])
% 
% 
% axis square






%% Supplementary Figure 4

% these series of analyses are only applicable to the sessions from Rat U
% and Rat V


%% to determine the dependence of PF stability on Post LT-PF correlation


% X = [ones(size(PFstability.data)) LTPFcorrs_pooled.maze.data];
% 
% [b, ~,~,~, stats] = regress(PFstability.data, X);
% pval = stats(3);
% r2   = stats(1);
% 
% estimatedd = b(1)+b(2)*LTPFcorrs_pooled.maze.data;
% residd = PFstability.data - estimatedd;



% y-axis variable

% yVar = PFstability.data;

allCorrs  = corr(LTs_pooled.remaze', LTs_pooled.maze_theta'); % correlation between MAZE and REMAZE theta LTs
yVar = diag(allCorrs);

    
% x-axis variable

currEpoch = 'post';

% xVar = LTPFcorrs_pooled.(currEpoch).data; 

allCorrs  = corr(LTs_pooled.(currEpoch)', LTs_pooled.maze_theta');
xVar = diag(allCorrs);


% Temporary 
% allCorrs  = corr(LTs_pooled.(currEpoch)', LTs_pooled.remaze');
% xVar = diag(allCorrs);




plotheight = 140;
plotwidth  = 140;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

scatter(xVar, yVar, 4, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5) % LTPFcorrs_pooled.(currEpoch).data, PFstability.data


ylabel('r(reMAZE thtLT, MAZE thtLT)', 'fontsize', fontsize)
xlabel('r(MAZE rplLT, MAZE thtLT)', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(xVar)) xVar];
[b, ~,~,~, stats] = regress(yVar, X); % PFstability.data
pval = stats(3);
r2   = stats(1);

yEst = b(1)+b(2)*xVar;
h1 = plot(xVar, yEst, 'color', [0 0 0 0.5], 'linewidth', 1);

if pval >= 0.001
   pval_txt = sprintf('p=%.3f', pval); 
else 
   pval_txt = sprintf('p=%.1e', pval);
end

% text(max(xVar), nanmedian(yEst), {sprintf('R2=%.2f', r2); pval_txt}, 'fontsize', 6)
text(-0.5, 1, {sprintf('R^2=%.2f', r2); pval_txt}, 'fontsize', 6)


xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 min(1, xl(2)+0.1)])
ylim([yl(1)- 0.1 min(1, yl(2)+0.1)])

axis square



selectUnits1 = find(PFstability.data > 0.5 & LTPFcorrs_pooled.(currEpoch).data > 0.5);

selectUnits2 = find(PFstability.data < 0.2 & LTPFcorrs_pooled.(currEpoch).data < 0.2);





%% multiple regression but this time regressing the PF stability against the LTs in different epochs

X = [LTPFcorrs_pooled.pre.data LTPFcorrs_pooled.maze.data LTPFcorrs_pooled.post.data];

lm = fitlm(X, PFstability);



%% to plot PF stability vs POST LT-PFcorrs residual 


plotheight = 120;
plotwidth  = 120;
fontsize   = 8;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on


scatter(residual_post_maze, PFstability.data, 2.66, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 0.5)

% scatter(residual_post_maze(selectUnits4), PFstability.data(selectUnits4), 2.66, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 1)
% 
% scatter(residual_post_maze(selectUnits3), PFstability.data(selectUnits3), 2.66, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.05, 'MarkerFaceAlpha', 1)


ylabel('r(reMAZE PF, MAZE PF)', 'fontsize', fontsize)
xlabel('re[r(POST LT, MAZE PF) | r(MAZE LT, MAZE PF)]', 'fontsize', fontsize)


set(gca, 'linewidth', 1, 'fontsize', fontsize, 'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.01])

X = [ones(size(residual_post_maze)) residual_post_maze];
[b, ~,~,~, stats] = regress(PFstability.data, X);
pval = stats(3);
r2   = stats(1);

estimated_PFstability = b(1)+b(2)*residual_post_maze;
h1 = plot(residual_post_maze, estimated_PFstability, 'color', [0 0 0 0.5], 'linewidth', 1);
text(max(residual_post_maze), nanmedian(estimated_PFstability), {sprintf('R2 = %.2f', r2); sprintf('p = %.1e', pval)}, 'fontsize', 4)


xl = xlim;
yl = ylim;

xlim([xl(1)- 0.1 min(1, xl(2)+0.1)])
ylim([yl(1)- 0.1 min(1, yl(2)+0.1)])

axis square



selectUnits1 = find(PFstability.data > 0.6 & residual_post_maze > 0.2);

selectUnits2 = find(PFstability.data < 0 & residual_post_maze < 0);



%% sub-functions

function customCDF(variable)


epochNames = fieldnames(variable);

% if isfield(variable.(epochNames{1}), 'data') % in case we are plotting the absolute values and not z-scores
%     
%     for iEpoch = 1:5
%         variable2.(epochNames{iEpoch}) = variable.(epochNames{iEpoch}).data;
%     end
%     
%     variable = variable2;
%     clear variabale2
% end 


colors = [[0 0 200]/255; [0 153 0]/255; [204 0 0]/255; [0 0 0]];

hold on

signedRank_pvalue = nan(4,1);
signedRank_zvalue = nan(4,1);

for iEpoch = 2:4

    currEpoch = epochNames{iEpoch};
    pooledVar = variable.(currEpoch).data;
    pooledVar_ui = variable.(currEpoch).ui;

%     [curr_p, ~, stats] = signrank(pooledVar, 0, 'tail', 'right');
%     signedRank_pvalue(iEpoch) = curr_p;
%     signedRank_zvalue(iEpoch) = stats.zval;

    med_data = nanmedian(pooledVar);
    med_ui   = nanmedian(pooledVar_ui);

    signedRank_pvalue(iEpoch) = numel(find(med_ui >= med_data))/size(med_ui, 2);


    
    [cdf_pooled, x_pooled] = ecdf(pooledVar);
    [x_pooled, idx] = unique(x_pooled);
    cdf_pooled = cdf_pooled(idx);
    
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
% 
% temp = [variable.maze variable.post variable.remaze];
% temp = temp(~isnan(sum(temp, 2)), :);
% 
% [p_friedman, anovaTab, stats] = friedman(temp);
% 


% significance scores


p_value = nan(4);
z_value  = nan(4);
for iEpoch = 2:4
   var1 = variable.(epochNames{iEpoch}).data;
   
   for jEpoch = iEpoch+1:5
       var2 = variable.(epochNames{jEpoch}).data;

       [p_value(iEpoch, jEpoch),~, stats] = signrank(var1, var2);
       z_value(iEpoch, jEpoch) = stats.zval;
   end
end
       
% add the significance signs here ...  
    

legend([ax(2) ax(3) ax(4)], 'location', 'best')

  
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




function pvalue = customCDF1(variable) 

pooledVar    = variable.data;
pooledVar_ui = variable.ui;


med_data = nanmedian(pooledVar);
med_ui   = nanmedian(pooledVar_ui);

pvalue = numel(find(med_ui >= med_data))/size(med_ui, 2);


if nargin == 1
    cl = [0 0 0];
end

if pvalue < 0.001
    sigSign = '***';
elseif pvalue < 0.01
    sigSign = '**';
elseif pvalue < 0.05
    sigSign = '*';
else
    sigSign = '';
end


[cdf_pooled, x_pooled] = ecdf(pooledVar);

[x_pooled, idx] = unique(x_pooled);
cdf_pooled = cdf_pooled(idx);

plot(x_pooled, cdf_pooled, 'color', [cl 0.8], 'linewidth', 1);

medianCorr = nanmedian(pooledVar);

 
ylim([0 1])
xlim([-1 1])
yl = ylim;
xl = xlim;

line([xl(1) medianCorr], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [cl 0.8])
line([medianCorr medianCorr], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', [cl 0.8])

text(med_data, 0.1, sprintf('%.2f%s',med_data, sigSign), 'fontsize', 8, 'color', [cl 0.8])


end


function annotPval = calAnnotPval(pval)

if pval < 0.0001
    annotPval = 'p < 10^-4';
elseif pval >= 0.0001
    annotPval = sprintf('p = %.3f', pval);
end

end
     
    


