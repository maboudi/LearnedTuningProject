
clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 12 15:17 6]; 
nSessions    = numel(currSessions);
epochNames = {'pre';'post'};


for iEpoch = 1:numel(epochNames)
    currEpoch = epochNames{iEpoch};

    meanCorr_w_posterior_pooled.(currEpoch) = cell(nSessions, 1);
    jump_from_prev_pooled.(currEpoch)       = cell(nSessions, 1);

    highJump_dist_meanPosteriorCorr.(currEpoch) = cell(nSessions, 1);
    lowJump_dist_meanPosteriorCorr.(currEpoch)  = cell(nSessions, 1);

    wc_entireSess.(currEpoch)    = cell(nSessions, 1);
    wc_ts_entireSess.(currEpoch) = cell(nSessions, 1);

end


for iSess = 1:nSessions
    
sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name
basePath      = fullfile(parentDir, sessionName);
    

load(fullfile(parentDir, sessionName, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')

spikes   = spikes_pyr;
nUnits   = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

peakPosBins = nan(nUnits, 1);
for iUnit = 1:nUnits
    peakPosBins(iUnit) = spikes(iUnit).peakPosBin.uni;
end

[~, sortIdx] = sort(peakPosBins, 'ascend');


load(fullfile(parentDir, sessionName, 'BayesianDecoding', [sessionName '.PBEInfo_replayScores.mat']), 'PBEInfo_replayScores')

PBEInfo = PBEInfo_replayScores;
PBEInfo = PBEInfo([PBEInfo.startT] < fileInfo.behavior.time(3,1) + 4*60*60);


nPBEs   = numel(PBEInfo);


%% surrogate distributions for each active unit counts

spatialTunings = nan(nUnits, nPosBins);
peakFR  = nan(nUnits, 1);
for unit = 1:nUnits
    spatialTunings(unit, :) = spikes(unit).spatialTuning_smoothed.uni;
    peakFR(unit)            = max(spatialTunings(unit, :));
end
spatialTunings = spatialTunings + 1e-4;
units_w_PFs = find(peakFR > 1);

% 
% if ismember('RL', fieldnames(spikes(unit).spatialTuning_smoothed))
%     spatialTunings_RL = nan(nUnits, nPosBins);
%     spatialTunings_LR = nan(nUnits, nPosBins);
%     for unit = 1:nUnits
%         spatialTunings_RL(unit, :) = spikes(unit).spatialTuning_smoothed.RL;
%         spatialTunings_LR(unit, :) = spikes(unit).spatialTuning_smoothed.LR;        
%     end
%     spatialTunings_RL = spatialTunings_RL + 1e-4;
%     spatialTunings_LR = spatialTunings_LR + 1e-4;
% end

cnct_20msTimeBin   = cell2mat({PBEInfo.fr_20msbin});
% cnct_20msTimeBin_onlyPFs = cnct_20msTimeBin(units_w_PFs, :);

% one way to incorpotate overall firing rates of each unitsi n the sampling process 
% would be to create a long list by concatenating particiapting units across all time bins and then sample from that list 

[partipation_by_bin_pooled, ~] = find(cnct_20msTimeBin); % _onlyPFs
participation_idx = 1:numel(partipation_by_bin_pooled);

eachTimeBinFiringRates = cnct_20msTimeBin(cnct_20msTimeBin > 0); % _onlyPFs

maxActiveUnitCount     = max(sum(cnct_20msTimeBin > 0)); % _onlyPFs
unitsParticipationRate = sum(cnct_20msTimeBin > 0, 2)/size(cnct_20msTimeBin, 2); % _onlyPFs  this will be used as a weight vector in the sampling process


allUnitCounts = 2:maxActiveUnitCount;
nIter         = 10000;
binDur        = 0.02;
surr_meanCorr_w_posterior = cell(maxActiveUnitCount, 1);
surr_mean_meanCorr_w_posterior = nan(maxActiveUnitCount, 1);
surr_std_meanCorr_w_posterior  = nan(maxActiveUnitCount, 1);

for uc = 1:maxActiveUnitCount-1

    curr_uc = allUnitCounts(uc);

    % all combinations of curr_uc
    allCombs = nan(nIter, curr_uc);
    for iu = 1:curr_uc
        allCombs(:, iu) = participation_idx(randperm(numel(participation_idx), nIter)); 
    end


%     allCombs = nchoosek(participation_idx, curr_uc); % this should reflect the firing rates of the units during the PBEs
%     allCombs = allCombs(randi(size(allCombs, 1), nIter, 1), :);

    partic_units = partipation_by_bin_pooled(allCombs);

%     nUniqPartUnits = zeros(nIter, 1);
%     for jj = 1:nIter
%         nUniqPartUnits = numel(unique(partic_units(jj, :)));
%     end
%     partic_units(nUniqPartUnits < curr_uc, :) = [];


    corresp_FRs  = eachTimeBinFiringRates(allCombs);

    surrBinIdx = repmat((1:nIter)', [1, curr_uc]);
    surrBinIdx = surrBinIdx';
    surrBinIdx = surrBinIdx(:);
    
    partic_units = partic_units';
    partic_units = partic_units(:);

    corresp_FRs = corresp_FRs';
    corresp_FRs = corresp_FRs(:);

    surrogate_20msBins = [];
    surrogate_20msBins = full(sparse(partic_units, surrBinIdx, corresp_FRs, nUnits, nIter));
    
    surr_postPr = baysDecoder(surrogate_20msBins, spatialTunings, binDur);
    surr_postPr = surr_postPr./repmat(sum(surr_postPr, 1), [nPosBins, 1]);

    


    % calculate the overlap between the posterior and the PFs

    surr_meanCorr_w_posterior{curr_uc} = nan(nIter, 1);
    
    for surr_bin = 1:nIter
    
        firingUnits   = find(surrogate_20msBins(:, surr_bin));
        
        currPosterior = surr_postPr(:, surr_bin);

        firingUnits_w_PFs = firingUnits(ismember(firingUnits, units_w_PFs));
    
        if numel(firingUnits_w_PFs) >= curr_uc

            if any(firingUnits_w_PFs)
                allCorrs = corr(spatialTunings(firingUnits_w_PFs, :)', currPosterior); 
                surr_meanCorr_w_posterior{curr_uc}(surr_bin) = nanmean(allCorrs);

%                 allCorrs = corr(spatialTunings(firingUnits_w_PFs, :)', spatialTunings(firingUnits_w_PFs, :)');
%                 allCorrs = setdiff(allCorrs(:), diag(allCorrs));
% 
%                 surr_meanCorr_w_posterior{curr_uc}(surr_bin) = nanmean(allCorrs);
            end
    
        else
            continue
        end
    
    end

    surr_mean_meanCorr_w_posterior(curr_uc) = nanmean(surr_meanCorr_w_posterior{curr_uc});
    surr_std_meanCorr_w_posterior(curr_uc)  = nanstd(surr_meanCorr_w_posterior{curr_uc});

end



%% correlation between the place field of acive units within each time bin
% and the decoded posterior within the same time bin 

meanCorr_w_posterior = cell(nPBEs, 1);
wc                   = cell(nPBEs, 1);
wc_ts                = cell(nPBEs, 1);
jump_from_prev       = cell(nPBEs, 1);
jump_to_next         = cell(nPBEs, 1);
jump                 = cell(nPBEs, 1);

for pbe = 1:nPBEs

%     postPr = PBEInfo(pbe).posteriorProbMat;
    postPr = baysDecoder(PBEInfo(pbe).fr_20msbin, spatialTunings, binDur);
    postPr = postPr./repmat(sum(postPr, 1), [nPosBins, 1]);


    weightedCorr = abs(PBEInfo(pbe).weightedCorr);
    weightedCorr_timeswap = PBEInfo(pbe).wc_ts;
    
    nTimeBins = size(postPr, 2);

    silentBins = sum(postPr, 1) == 0;

    [~, peakPosBins] = max(postPr);
    peakPosBins(silentBins) = nan;
    

    jump_from_prev{pbe}        = nan(nTimeBins, 1);
    jump_from_prev{pbe}(2:end) = abs(peakPosBins(2:end) - peakPosBins(1:end-1))/nPosBins; 

    jump_to_next{pbe}          = nan(nTimeBins, 1);
    jump_to_next{pbe}(1:end-1) = abs(peakPosBins(1:end-1) - peakPosBins(2:end))/nPosBins; 

%     jump{pbe} = min(jump_from_prev{pbe}, jump_to_next{pbe});
    medianJumpDistance = nanmedian(jump_from_prev{pbe});
    jump{pbe}          = repelem(medianJumpDistance, nTimeBins)';

    wc{pbe}    = weightedCorr * ones(nTimeBins, 1);
    wc_ts{pbe} = weightedCorr_timeswap*ones(nTimeBins, 1);
    

    meanCorr_w_posterior{pbe} = nan(nTimeBins, 1);
    
    for bin = 1:nTimeBins
    
        firingUnits   = find(PBEInfo(pbe).fr_20msbin(:, bin));
        nFiringUnits  = numel(firingUnits);

        currPosterior = postPr(:, bin);

        firingUnits_w_PFs = firingUnits(ismember(firingUnits, units_w_PFs));

        if numel(firingUnits_w_PFs) >= 2

            if any(firingUnits_w_PFs)

                allCorrs = corr(spatialTunings(firingUnits_w_PFs, :)', currPosterior); 
                meanCorr_w_posterior{pbe}(bin) = (nanmean(allCorrs) - surr_mean_meanCorr_w_posterior(nFiringUnits))/surr_std_meanCorr_w_posterior(nFiringUnits);

%                 allCorrs = corr(spatialTunings(firingUnits_w_PFs, :)', spatialTunings(firingUnits_w_PFs, :)');
%                 allCorrs = setdiff(allCorrs(:), diag(allCorrs));

%                 meanCorr_w_posterior{pbe}(bin) = (nanmean(allCorrs) - surr_mean_meanCorr_w_posterior(nFiringUnits))/surr_std_meanCorr_w_posterior(nFiringUnits);

            end
    
        else
            continue
        end
    
    end

end


pooledJumpDistance_entire = cell2mat(jump); % _from_prev
lowThresh_jump  = prctile(pooledJumpDistance_entire, 50);
highThresh_jump = prctile(pooledJumpDistance_entire, 50);


epochs = {PBEInfo.epoch};

for epoch = 1:2
    currEpoch = epochNames{epoch};
    
    pbe_inclusion_idx = strcmp(epochs, currEpoch);
    
    meanCorr_w_posterior_pooled.(currEpoch){iSess} = cell2mat(meanCorr_w_posterior(pbe_inclusion_idx));
    jump_from_prev_pooled.(currEpoch){iSess}       = cell2mat(jump(pbe_inclusion_idx)); %_from_prev
    
%     idx = ~isnan(meanCorr_w_posterior_pooled) & ~isnan(jump_from_prev_pooled);
%     [a,b] = corr(meanCorr_w_posterior_pooled(idx), jump_from_prev_pooled(idx));

    highJump_dist_meanPosteriorCorr.(currEpoch){iSess} = meanCorr_w_posterior_pooled.(currEpoch){iSess}(jump_from_prev_pooled.(currEpoch){iSess} > highThresh_jump);
    lowJump_dist_meanPosteriorCorr.(currEpoch){iSess}  = meanCorr_w_posterior_pooled.(currEpoch){iSess}(jump_from_prev_pooled.(currEpoch){iSess} < lowThresh_jump);

    wc_entireSess.(currEpoch){iSess} = cell2mat(wc(pbe_inclusion_idx));
    wc_ts_entireSess.(currEpoch){iSess} = cell2mat(wc_ts(pbe_inclusion_idx));

end



end

%%
for iEpoch = 1:2

    currEpoch = epochNames{iEpoch};

    mean_posteriorCorr.(currEpoch)   = cell2mat(meanCorr_w_posterior_pooled.(currEpoch));
    jumpDistance.(currEpoch)         = cell2mat(jump_from_prev_pooled.(currEpoch));

    hj_meanPosteriorCorr.(currEpoch) = cell2mat(highJump_dist_meanPosteriorCorr.(currEpoch));
    lj_meanPosteriorCorr.(currEpoch) = cell2mat(lowJump_dist_meanPosteriorCorr.(currEpoch));

    wc_pooled.(currEpoch) = cell2mat(wc_entireSess.(currEpoch));
    wc_ts_pooled.(currEpoch) = cell2mat(wc_ts_entireSess.(currEpoch));

end


% The second way: pool across sessions and then take the median jump as the threshold 

allJumpDistances = [jumpDistance.pre; jumpDistance.post];
j_thresh_l = prctile(allJumpDistances, 50);
j_thresh_h = prctile(allJumpDistances, 50);


% allPosteriorCorr = [mean_posteriorCorr.pre; mean_posteriorCorr.post];

allWCs = [wc_pooled.pre; wc_pooled.pre];
wc_thresh = prctile(allWCs, 50);



plotwidth  = 50;
plotheight = 50;
fontsize   = 5;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

colors = {'#005CE9';'#D40035'};

hold on

for iEpoch = 1:2
    
    currEpoch = epochNames{iEpoch};
    
    currVariable = mean_posteriorCorr.(currEpoch);
    
    
    lowJump   = wc_ts_pooled.(currEpoch) > 75 & jumpDistance.(currEpoch) < j_thresh_l;
    highJump  = wc_ts_pooled.(currEpoch) < 25;

%     lowJump   = jumpDistance.(currEpoch) < j_thresh_l;
%     highJump  = jumpDistance.(currEpoch) >= j_thresh_h;
    
    
    lowJumpData.(currEpoch) = currVariable(lowJump); 
%     [xData, yData] = myViolin_oneSided(lowJumpData.(currEpoch), 0.03);
%     h1 = patch(iEpoch+20*xData, yData, 'k', 'FaceColor', colors{iEpoch}, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
%     
%     lq  = prctile(lowJumpData.(currEpoch), 25);
%     hq  = prctile(lowJumpData.(currEpoch), 75);
    med = prctile(lowJumpData.(currEpoch), 50);
% 
%     [med hq-lq]
%     
%     line([iEpoch iEpoch+0.25], [lq lq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
%     line([iEpoch iEpoch+0.25], [hq hq], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)
%     line([iEpoch iEpoch+0.25], [med med], 'linestyle', '-', 'color', 'w', 'linewidth', 0.5)

    [curr_cdf, curr_x] = ecdf(lowJumpData.(currEpoch));
    [curr_x, idx] = unique(curr_x);
    curr_cdf = curr_cdf(idx);
    
    plot(curr_x, curr_cdf, 'color', colors{iEpoch}, 'linewidth', 0.5, 'DisplayName', [currEpoch '-high'])

    ylim([0 1])
    xlim([-4 4])
    yl = ylim;
    xl = xlim;

    line([xl(1) med], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors{iEpoch})
    line([med med], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors{iEpoch})


    highJumpData.(currEpoch) = currVariable(highJump); 
%     [xData, yData] = myViolin_oneSided(highJumpData.(currEpoch), 0.03); 
%     h2 = patch(iEpoch-20*xData, yData, 'k', 'FaceColor', 'none','EdgeColor', colors{iEpoch}, 'FaceAlpha', 0.8);
%     
%     lq  = prctile(highJumpData.(currEpoch), 25);
%     hq  = prctile(highJumpData.(currEpoch), 75);
    med = prctile(highJumpData.(currEpoch), 50);
% 
%     [med hq-lq]
% 
%     line([iEpoch-0.25 iEpoch], [lq lq], 'linestyle', '-', 'color', colors{iEpoch}, 'linewidth', 0.5)
%     line([iEpoch-0.25 iEpoch], [hq hq], 'linestyle', '-', 'color', colors{iEpoch}, 'linewidth', 0.5)
%     line([iEpoch-0.25 iEpoch], [med med], 'linestyle', '-', 'color', colors{iEpoch}, 'linewidth', 0.5)
    [curr_cdf, curr_x] = ecdf(highJumpData.(currEpoch));
    [curr_x, idx] = unique(curr_x);
    curr_cdf = curr_cdf(idx);
    
    plot(curr_x, curr_cdf, 'color', colors{iEpoch}, 'linewidth', 0.5, 'LineStyle', ':', 'DisplayName', [currEpoch '-low'])

    ylim([0 1])
    xlim([-4 4])
    yl = ylim;
    xl = xlim;

    line([xl(1) med], [0.5 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors{iEpoch})
    line([med med], [yl(1) 0.5], 'linewidth', 0.5, 'linestyle', ':', 'color', colors{iEpoch})


    

    [pval.(currEpoch), ~,stats] = ranksum(lowJumpData.(currEpoch), highJumpData.(currEpoch));
    z_value.(currEpoch) = stats.zval;

end

[pval.l_pre, ~,stats] = ranksum(lowJumpData.pre, 0);
z_value.l_pre_h_post = stats.zval;

[pval.h_post, ~,stats] = ranksum(lowJumpData.post, 0);
z_value.h_post = stats.zval;


[pval.l_pre_h_post, ~,stats] = ranksum(lowJumpData.pre, highJumpData.post);
z_value.l_pre_h_post = stats.zval;

[pval.h_pre_h_post, ~,stats] = ranksum(highJumpData.pre, highJumpData.post);
z_value.h_pre_h_post = stats.zval;

[pval.l_pre_l_post, ~,stats] = ranksum(lowJumpData.pre, lowJumpData.post);
z_value.l_pre_l_post = stats.zval;

[pval.h_pre_l_post, ~,stats] = ranksum(highJumpData.pre, lowJumpData.post);
z_value.h_pre_l_post = stats.zval;

ax = gca;
ax.YGrid = 'on';

% xlim([0.5 2.5])



% legend([h1, h2], 'stable LTs', 'nonstable LTs', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')

% xlabel('epoch', 'fontsize', fontsize)
% ylabel('PF-LT Pearson correlation coeff. (z)', 'fontsize', fontsize)
% ylabel('PF-LT KL divergence (z)', 'fontsize', fontsize)


% ylim([-0.5 1])
axis square
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', [xl(1) 0 xl(2)], 'TickDir', 'out','TickLength',[0.01, 0.01] ) %'ytick', [-0.5 0 0.5 1]




%% functions

function [xData, yData] = myViolin_oneSided(data, binSize)

% data = data(data > prctile(data, 5) & data < prctile(data, 97.5));

% binEdges = (min(data)-eps):binSize:(max(data)+eps);

binEdges = -5:binSize:5;
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

