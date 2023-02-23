clear
clc
close all


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/';

rr = dir(parentDir);


currSessions = [1:5 7:9 10 13:15 6];
% currSessions = 7;
nSessions    = numel(currSessions);

nShuffles = 1000;

startT = cell(nSessions, 1);
endT   = cell(nSessions, 1);


spikes_pooled     = cell(nSessions, 1);
binCenters_pooled = cell(nSessions, 1);
binCenters_nonOverlap_pooled = cell(nSessions, 1);


assemblyTunings_time_pooled            = cell(nSessions, 1);
assemblyTunings_time_nonOverlap_pooled = cell(nSessions, 1);

assemblyTuningPFcorr_zscore_pooled = cell(nSessions, 1);
assemblyTuningPFcorr_pooled        = cell(nSessions, 1);


% assemblyTuningPFKLdiv_pooled        = cell(nSessions, 1);
% assemblyTuningPFKLdiv_zscore_pooled = cell(nSessions, 1);



epochsXCorr         = cell(nSessions, 1);
epochsXCorr_zscore  = cell(nSessions, 1);

epochsKLdiv        = cell(nSessions, 1);
epochsKLdiv_zscore = cell(nSessions, 1);


% epochsXCorr_pval    = cell(nSessions, 1);
% epochsXCorr_ui_pval = cell(nSessions, 1);


activeUnits = cell(nSessions, 1);

epochNames = {'pre'; 'run'; 'post'};



%%


for iSess = 1%:numel(currSessions)

sessionNumber = currSessions(iSess);
sessionName   = rr(sessionNumber+2).name

basePath      = fullfile(parentDir, sessionName);



% spikes and behavior data
load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes_pyr', 'fileInfo')
behavior = fileInfo.behavior;
behavior.time = behavior.time/3600;

startT{iSess}.pre  = behavior.time(1,1); endT{iSess}.pre  = behavior.time(2,1); 
startT{iSess}.run  = behavior.time(2,1); endT{iSess}.run  = behavior.time(2,2); 
startT{iSess}.post = behavior.time(2,2); endT{iSess}.post = behavior.time(3,1)+4; % in long recordings only the first 4 hours of post sleep is used in this figures 


spikes = spikes_pyr;
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);

nUnits = numel(spikes);



% bidirectional spatial tuning
spatialTunings_merge = zeros(nUnits, nPosBins);
for iUnit = 1:nUnits
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;    
end

spikes_pooled{iSess} = spikes;




% % % load learned tunings

% 15-minutes time windows

s = load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_time4.mat']), 'binCenters', 'assemblyTunings_time', 'activeUnits');

binCenters            = s.binCenters;
assemblyTunings_time  = s.assemblyTunings_time;
activeUnits{iSess}    = s.activeUnits;


for iEpoch = 1:numel(epochNames)
    binCenters.(epochNames{iEpoch}) = binCenters.(epochNames{iEpoch})/3600;
end

binCenters_pooled{iSess} = binCenters;

okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);



% % % alternative unit-ID shuffle for the 15 minutes time windows, we need
% % this new shuffle for checking the stability later
% 
% 
% for iEpoch = 1:3
%     epochName = epochNames{iEpoch};
%     
%     
%     currData = assemblyTunings_time.(epochName).data;
%     nTimeBins = size(currData, 3);
%     
%     assemblyTunings_time.(epochName).ui = nan([size(currData), nShuffles]);
%     for inst = 1:nShuffles
%         
%         for iTimeBin = 1:nTimeBins
% 
%             for iUnit = 1:numel(okUnits)
% 
%                 currUnit = okUnits(iUnit);
%                 otherUnits = setdiff(okUnits, currUnit);
%                 assemblyTunings_time.(epochName).ui(currUnit, :, iTimeBin, inst) = currData(otherUnits(randperm(numel(otherUnits), 1)), :, iTimeBin); % randperm(nTimeBins, 1)
%             end
%         end
%     end
% end



% remove the time windows with overlap
    
for iEpoch = 1:numel(epochNames)
    epochName = epochNames{iEpoch};
    
    binCenters_nonOverlap.(epochName) = binCenters.(epochName)(1:3:end);
    
    assemblyTunings_time_nonOverlap.(epochName).data = assemblyTunings_time.(epochName).data(:, :, 1:3:end);
    assemblyTunings_time_nonOverlap.(epochName).ui  = assemblyTunings_time.(epochName).ui(:, :, 1:3:end, :); 
end



% truncate the long recording sessions to the first 4 hours of post sleep

idx = find(binCenters.post < endT{iSess}.post);

binCenters.post = binCenters.post(idx);
assemblyTunings_time.post.data = assemblyTunings_time.post.data(:, :, idx);
assemblyTunings_time.post.ui   = assemblyTunings_time.post.ui(:, :, idx, :);


idx = find(binCenters_nonOverlap.post < endT{iSess}.post);
binCenters_nonOverlap.post = binCenters_nonOverlap.post(idx);
assemblyTunings_time_nonOverlap.post.data = assemblyTunings_time_nonOverlap.post.data(:, :, idx);
assemblyTunings_time_nonOverlap.post.ui   = assemblyTunings_time_nonOverlap.post.ui(:, :, idx, :);

binCenters_nonOverlap_pooled{iSess} = binCenters_nonOverlap;


assemblyTunings_time_pooled{iSess} = assemblyTunings_time;
assemblyTunings_time_nonOverlap_pooled{iSess} = assemblyTunings_time_nonOverlap;


load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']), 'assemblyTunings', 'assemblyTuningPFcorr', 'assemblyTuningPFcorr_zscore');




%% calculate an alternative shuffle PF correlations


for iEpoch = 1:3
    currEpoch = epochNames{iEpoch};
    
    currLT = assemblyTunings.(currEpoch).data;
         
%     corrMat = corr(currLT(okUnits, :)', spatialTunings_merge(okUnits, :)', 'type', 'pearson');

    corrMat = corr(currLT', spatialTunings_merge', 'type', 'pearson');
    KLdivergence = calKLDivergence(currLT', spatialTunings_merge');
    
    
    % % % data
    assemblyTuningPFcorr.(currEpoch).data  = diag(corrMat);
%     assemblyTuningPFKLdiv.(currEpoch).data = diag(KLdivergence);
    

    % % % shuffle: basically every non-disgoanal element of the correlation
    % matrix
    
    % PF correlation
    for ii = 1:nUnits
        corrMat(ii, ii) = nan;
        KLdivergence(ii, ii) = nan;
    end
    
    nonMatchingUnitCorrs = corrMat(:);
    nonMatchingUnitCorrs(isnan(nonMatchingUnitCorrs)) = []; 

    assemblyTuningPFcorr.(currEpoch).ui     = nonMatchingUnitCorrs(randperm(numel(nonMatchingUnitCorrs), min(numel(nonMatchingUnitCorrs), nShuffles)));
    assemblyTuningPFcorr_zscore.(currEpoch) = (assemblyTuningPFcorr.(currEpoch).data - nanmean(assemblyTuningPFcorr.(currEpoch).ui))./nanstd(assemblyTuningPFcorr.(currEpoch).ui);
     
    
    
%     % KL divergence
%     for ii = 1:nUnits
%         KLdivergence(ii, ii) = nan;
%     end
%     
%     nonMatchingUnitKLs = KLdivergence(:);
%     nonMatchingUnitKLs(isnan(nonMatchingUnitKLs)) = []; 
% 
%     assemblyTuningPFKLdiv.(currEpoch).ui     = nonMatchingUnitKLs(randperm(numel(nonMatchingUnitKLs), min(numel(nonMatchingUnitKLs), nShuffles))); 
%     assemblyTuningPFKLdiv_zscore.(currEpoch) = (assemblyTuningPFKLdiv.(currEpoch).data - nanmean(assemblyTuningPFKLdiv.(currEpoch).ui))./nanstd(assemblyTuningPFKLdiv.(currEpoch).ui);    
    
end      

assemblyTuningPFcorr_pooled{iSess}        = assemblyTuningPFcorr;
assemblyTuningPFcorr_zscore_pooled{iSess} = assemblyTuningPFcorr_zscore;


% assemblyTuningPFKLdiv_pooled{iSess}        = assemblyTuningPFKLdiv;
% assemblyTuningPFKLdiv_zscore_pooled{iSess} = assemblyTuningPFKLdiv_zscore;




% % calculate the PF correlations for data and ui for the 15-minutes moving
% time windows
% 
% dataTypes = {'data'; 'ui'};
% for iDataType = 1:2
%     currDataType = dataTypes{iDataType};
%     
%     for iEpoch = 1:3
%         currEpoch = epochNames{iEpoch};
%         
%         allLTs = assemblyTunings_time_nonOverlap.(currEpoch).(currDataType);
%         
%         nTimeBins = size(allLTs, 3);
%         nInstances = size(allLTs, 4); 
%         
%         assemblyTuningPFcorr_time.(currEpoch).(currDataType) = nan(nUnits, nTimeBins, nInstances);
%         
%         for ins = 1:nInstances
%         
%             for iTimeBin = 1:nTimeBins 
%                 currLT = allLTs(:, :, iTimeBin, ins);
%                 corrMat = corr(currLT', spatialTunings_merge', 'type', 'pearson');
%                 assemblyTuningPFcorr_time.(currEpoch).(currDataType)(:, iTimeBin, ins) = diag(corrMat);
%             end            
%         end
% 
%     end
%     
% end




%% for assemblyTunings_time calculate the correlation matrix of learned tunings

dataTypes = {'data'; 'ui'};
epochNames = {'pre'; 'run'; 'post'};

epochsXCorr{iSess}.data   = nan(3, 3, nUnits);
epochsXCorr{iSess}.ui     = nan(3, 3, nUnits, nShuffles);
epochsXCorr_zscore{iSess} = nan(3, 3, nUnits); 


% epochsKLdiv{iSess}.data   = nan(3,3, nUnits);
% epochsKLdiv{iSess}.ui     = nan(3,3, nUnits);
% epochsKLdiv_zscore{iSess} = nan(3,3, nUnits);


nInstances = struct('data', 1, 'ui', nShuffles);

        
for iEpoch = 1:3
    
    iEpoch
    
    currData1 = assemblyTunings_time_nonOverlap.(epochNames{iEpoch}).data;

    for jEpoch = iEpoch:3
        currData2 = assemblyTunings_time_nonOverlap.(epochNames{jEpoch}).data;

        for iUnit = 1:nUnits
            
%             currUnit   = okUnits(iUnit);
            currUnit   = iUnit;

%             otherUnits = setdiff(okUnits, currUnit);
            otherUnits = setdiff(1:nUnits, currUnit);
            
            
            % data
            tuning1 = squeeze(currData1(currUnit, :, :));
            tuning2 = squeeze(currData2(currUnit, :, :));
            
            corrMat      = corr(tuning1, tuning2, 'type', 'pearson');
            KLdivergence = calKLDivergence(tuning1, tuning2);
            
            if jEpoch == iEpoch           
                for ii = 1:size(corrMat, 1)
                    corrMat(ii, ii)= nan;
                    KLdivergence(ii, ii) = nan;
                end
            end
            
            epochsXCorr{iSess}.data(iEpoch, jEpoch, currUnit) = nanmedian(corrMat(:));
%             epochsKLdiv{iSess}.data(iEpoch, jEpoch, currUnit) = nanmedian(KLdivergence(:));
            
            
            
            % shuffle
            for inst = 1:nInstances.(currDataType)
                
                % generate shuffle tunings
                
                nTimeBins1 = size(currData1, 3);
                tuning1 = zeros(nPosBins, nTimeBins1);
                for iTimeBin = 1:size(currData1, 3)
                    tuning1(:, iTimeBin) = squeeze(currData1(otherUnits(randperm(numel(otherUnits), 1)), :, iTimeBin));
                end
                
                
                nTimeBins2 = size(currData2, 3);
                tuning2 = zeros(nPosBins, nTimeBins2);
                for iTimeBin = 1:size(currData2, 3) 
                    tuning2(:, iTimeBin) = squeeze(currData2(otherUnits(randperm(numel(otherUnits), 1)), :, iTimeBin));
                end


                corrMat = corr(tuning1, tuning2, 'type', 'pearson');
%                 KLdivergence = calKLDivergence(tuning1, tuning2);


                if jEpoch == iEpoch
                    for ii = 1:size(corrMat, 1)
                        corrMat(ii, ii)= nan;
%                         KLdivergence(ii, ii) = nan;
                    end
                end

                epochsXCorr{iSess}.ui(iEpoch, jEpoch, currUnit, inst) = nanmedian(corrMat(:));
%                 epochsKLdiv{iSess}.ui(iEpoch, jEpoch, currUnit, inst) = nanmedian(KLdivergence(:));

            end

        end
    end    
 
end

epochsXCorr_zscore{iSess} = ((epochsXCorr{iSess}.data) - mean(abs(epochsXCorr{iSess}.ui), 4)) ./ std(abs(epochsXCorr{iSess}.ui), [], 4);
% epochsKLdiv_zscore{iSess} = ((epochsKLdiv{iSess}.data) - mean(abs(epochsKLdiv{iSess}.ui), 4)) ./ std(abs(epochsKLdiv{iSess}.ui), [], 4);


end




%% plot the learning tunings across time for example units in an example session

iSess = 1;


spikes = spikes_pooled{iSess};

selectedUnits = [72 117 79 63 116 90]; % AChilles
% selectedUnits = [5 15 20 25 35 55]; % rat U


plotwidth  = 360;
plotheight = 235;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

nor = 3;
noc = 2;

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

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time_pooled{iSess}.(currEpoch).data(currUnit, :, :, 1));
        learnedTuning{iEpoch} = learnedTuning{iEpoch} ./ repmat(max(learnedTuning{iEpoch}, [], 1), [size(learnedTuning{iEpoch}, 1) 1]);

        imagesc(binCenters_nonOverlap_pooled{iSess}.(currEpoch), 1:nPosBins, learnedTuning{iEpoch}); colormap('jet')

    end
    
    plot(0.1*currSpatialTuning + startT{iSess}.run , 1:nPosBins, 'linewidth', 1.5, 'color', [0 0 0 0.6], 'DisplayName', 'track spatial tuning')


    yl = ylim(ax(iUnit));

    
    if iUnit == 4
        line(0.2*[0 2] + startT{iSess}.run, [yl(2)-10 yl(2)-10], 'linewidth', 1.5, 'color', [0 0 0 0.6])
        text(startT{iSess}.run, yl(2)+20, '2 Hz', 'fontsize', fontsize, 'color', 'k')
    end
    
    ylim([0 nPosBins]);
    yl = ylim;
    text(0, yl(2)+12, ['unit ' num2str(iUnit)], 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'left')
    
    xlim([0 binCenters_nonOverlap_pooled{iSess}.post(end)])
    
    set(ax(iUnit), 'YTick', [1 nPosBins], 'YTickLabel', [0 1], 'TickDir', 'out','TickLength',[0.01, 0.01])
    
    if ismember(iUnit, [1 4])
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           currBinCenters = binCenters_nonOverlap_pooled{iSess}.(currEpoch);
           text(median(currBinCenters), yl(2)+20, currEpoch, 'fontsize', fontsize, 'color', 'k', 'HorizontalAlignment', 'center')
       end
    end
        
    
    if iUnit == 2
        ylabel('position(normalized)', 'fontsize', fontsize)
    end
    
    if ismember(iUnit, [3 6])
        xlabel('time (hour)')
    end
    
    if iUnit >= 4
        set(ax(iUnit), 'YTickLabel', {})
    end
    
    if ~ismember(iUnit, [3 6])
       set(ax(iUnit), 'XTickLabel', {}) 
    end
    
    hold off
    
end



%% plot cross-correlation matrix of learned tunings for the entire session


plotwidth  = 375;
plotheight = 175;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


cnctBinCenters = [binCenters_pooled{iSess}.pre binCenters_pooled{iSess}.run binCenters_pooled{iSess}.post];
tickLabels = 0:2:cnctBinCenters(end);

for iUnit = 1:numel(selectedUnits)
    
    
    currUnit = selectedUnits(iUnit);
    
    
    learnedTuning = cell(3,1); 
    for iEpoch = 1:3

        currEpoch = epochNames{iEpoch};

        learnedTuning{iEpoch} = squeeze(assemblyTunings_time_pooled{iSess}.(currEpoch).data(currUnit, :, :, 1));

    end
    
    tmpConcat = cell2mat(learnedTuning');

    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');
    hold on
    
    
    
    for iEpoch = 1:3
        for jEpoch = 1:3
            
            corrMat = corr(learnedTuning{jEpoch}, learnedTuning{iEpoch});
            corrMat(isnan(corrMat)) = 0;
            imagesc(binCenters_pooled{iSess}.(epochNames{iEpoch}), binCenters_pooled{iSess}.(epochNames{jEpoch}), corrMat)
 
        end
    end
    caxis([-1 1])
    colormap('jet')
    

    ylim([0 cnctBinCenters(end)])
    xlim([0 cnctBinCenters(end)])
    
    yl = ylim;
    text(0, yl(2)+0.5, sprintf('unit %d', iUnit), 'fontsize', fontsize, 'HorizontalAlignment', 'left')
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    if iUnit == 1
       for iEpoch = 1:3
           currEpoch = epochNames{iEpoch};
           h = text(-4, median(binCenters_pooled{iSess}.(currEpoch)), currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           h = text(-2, median(cnctBinCenters), 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center');
           set(h, 'rotation', 90)
           
           
           text(median(binCenters_pooled{iSess}.(currEpoch)), -4 , currEpoch, 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 
           
          text(median(cnctBinCenters), -2 , 'time(hour)', 'fontsize', fontsize, 'HorizontalAlignment', 'center'); 

       end
    end
    
    
    xticks(tickLabels)
    yticks(tickLabels)
    axis square
    
    if iUnit == 6
        h = colorbar;
        h.Label.String = 'correlation';
    end
    set(gca, 'box', 'off', 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'TickDir', 'out','TickLength',[0.01, 0.01])
   
    
end





%% plot the distribution of mean cross-correlations for the example units



% for the example units (distribution of correlation for real data and ui)

plotwidth  = 330;
plotheight = 150;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nor = 1;
noc = 6;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 5;
gapr = 0;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


for iUnit = 1:numel(selectedUnits)
    
    currUnit = selectedUnits(iUnit);
    
    ax(iUnit) = axes('position',sub_pos{iUnit},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
    
    
    tt = 0;
    for iEpoch = [1 3] % pre and post within-epoch consistency of learned tunings over time
        tt = tt + 1; 
        % data
        data = epochsXCorr{iSess}.data(iEpoch, iEpoch, currUnit);
        
%         [xData, yData] = myViolin(data);
% 
%         patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')

        % ui
        ui = squeeze(epochsXCorr{iSess}.ui(iEpoch, iEpoch, currUnit, :));
        
        
        try 
        [xData, yData] = myViolin(ui);
        patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
        scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')

        catch
            continue
        end
        
    end
    
    tt = tt+1;
    
    % similarity of learned tunings between pre and post
    
    % data
    data = epochsXCorr{iSess}.data(1, 3, currUnit, 1);
        
%     [xData, yData] = myViolin(data);
% 
%     h1 = patch(tt-0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    h1 =scatter(tt-0.15, data, 5, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none'); 
    
    
    % ui

    ui = squeeze(epochsXCorr{iSess}.ui(1, 3, currUnit, :));
    
    try
    [xData, yData] = myViolin(ui);    
    h2 = patch(tt+0.15+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    scatter(tt+0.15, nanmedian(ui), 2, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none');
    catch
        continue
    end
    
    if iUnit == 1
        ylabel({'learned tuning'; 'stability'}, 'fontsize', fontsize)
    end
    
    if iUnit == 6
       legend([h1, h2], 'data', 'shuffle', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')
    end
    
    xlim([0.5 3.5])
%     ylim([-1 1])
%     grid on
    
    if iUnit >= 2
        set(gca, 'YTickLabel', {})
    end
    
    set(ax(iUnit), 'fontsize', fontsize, 'position', sub_pos{iUnit}, 'linewidth', 1, 'box', 'off', 'xtick', 1:3, 'xticklabel', {'pre'; 'post'; 'pre-post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
    xtickangle(ax(iUnit), 45)
    set(gca, 'ytick', -1:0.5:1)
    
    aa = gca;
    aa.YGrid = 'on';
    
end

linkaxes(ax, 'y')
yl = ylim;

ylim([yl(1)-0.1 yl(2)+0.1])




%% plot the distribution of cross-correlation for an example session (panel D)

iSess = 1;

plotwidth  = 70;
plotheight = 90;
fontsize   = 6;


okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

hold on

tt = 1;
for iEpoch = [1 3]
    
    
%     currMeds = nan(numel(okUnits), 1);
%     for iUnit = 1:numel(okUnits)
%         currData = squeeze(epochsXCorr{iSess}.data{iEpoch, iEpoch, okUnits(iUnit), 1});  
%         currMeds(iUnit) = nanmedian(currData);
%     end

    currMeds = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits, 1));

%     signIdx  = squeeze(epochsXCorr_ui_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4 & squeeze(epochsXCorr_pval{iSess}(iEpoch, iEpoch, okUnits)) < 1e-4;
%     scatter(tt+0.05*randn(numel(find(~signIdx)), 1), currMeds(~signIdx), 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
%     scatter(tt+0.05*randn(numel(find(signIdx)), 1), currMeds(signIdx), 1, 'MarkerFaceColor', 'g', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);

    scatter(tt+0.05*randn(numel(currMeds), 1), currMeds, 1, 'MarkerFaceColor', 'k', 'markerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
    
    [xData, yData] = myViolin(currMeds);
    patch(tt+2*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    
    tt = tt + 1;
end

xlim([0.5 2.5])

set(gca, 'box', 'off', 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
xtickangle(gca, 45)
ylabel({'median cross-correlation'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';



%% correlation between .data and _zscore

iEpoch = 1;
jEpoch = 2;

iSess = 6;
okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


figure; scatter(epochsXCorr{iSess}.data(iEpoch, jEpoch, okUnits, 1), epochsXCorr_zscore{iSess}(iEpoch, jEpoch, okUnits, 1), 5, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none')


grid on



%% plot the distribution of cross-corrrelation for all sessions

plotwidth  = 90;
plotheight = 90;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


hold on

nSessions = size(epochsXCorr_zscore, 1);

tt = 1;
for iEpoch = [1 3]
        
    currData = cell(nSessions, 1);
    for iSess = 1:nSessions
        curr_x = tt+(iSess-ceil(nSessions/2))*0.85/nSessions;
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);


        currData{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));

        med = nanmedian(currData{iSess});
        lq  = nanprctile(currData{iSess}, 25);
        hq  = nanprctile(currData{iSess}, 75);
        data_min = min(currData{iSess});
        data_max = max(currData{iSess});


        scatter(curr_x, med, 2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

        % whiskers
        line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 0 0 0.8])
        line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 0 0 0.8])

        % box
        patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    end
    
    pooledData = cell2mat(currData);

    med = nanmedian(pooledData);
    lq  = nanprctile(pooledData, 25);
    hq  = nanprctile(pooledData, 75);
    data_min = min(pooledData);
    data_max = max(pooledData);

    curr_x = tt+(iSess+1-ceil(nSessions/2))*0.85/nSessions;

    scatter(curr_x, med, 2, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8)

    % whiskers
    line([curr_x curr_x], [data_min lq], 'linewidth', 0.25, 'color', [0 1 0 0.8])
    line([curr_x curr_x], [data_max hq], 'linewidth', 0.25, 'color', [0 1 0 0.8])

    % box
    patch([curr_x-0.02 curr_x+0.02 curr_x+0.02 curr_x-0.02], [lq lq hq hq], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.8)
    
    tt = tt + 1;
end

xlim([0.5 2.5])

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', 1:2, 'xticklabel', {'pre';'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01])
xtickangle(gca, 45)
ylabel({'mean'; 'cross-correlation(z)'}, 'fontsize', fontsize)
ax = gca;
ax.YGrid = 'on';





%% correlation between LT consistency during an epoch and its correlation with PFs


plotwidth  = 50;
plotheight = 130;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);


nSessions = size(epochsXCorr_zscore, 1);
thresh_z = 4;

clear ax med LTstability_pooled PFcorrs_pooled

tt = 1;
for iEpoch = 1
    
    hold on
       
    
    LTstability  = cell(nSessions, 1); % LT consistency
    PFcorrs = cell(nSessions, 1);
    
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
%         PFcorrs{iSess}     = assemblyTuningPFcorr_zscore_pooled{iSess}.(epochNames{iEpoch})(okUnits);
        PFcorrs{iSess}     = squeeze(epochsXCorr_zscore{iSess}(1, 3, okUnits));

    end
    
    
    LTstability_pooled.(epochNames{iEpoch}) = cell2mat(LTstability);
    PFcorrs_pooled.(epochNames{iEpoch}) = cell2mat(PFcorrs);
    
    [xData, yData] = myViolin_oneSided(PFcorrs_pooled.(epochNames{iEpoch})(LTstability_pooled.(epochNames{iEpoch}) > thresh_z));
    h1 = patch(tt+5*xData, yData, [1 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    med.(epochNames{iEpoch}) = [nanmedian(PFcorrs_pooled.(epochNames{iEpoch})(LTstability_pooled.(epochNames{iEpoch}) < thresh_z)) nanmedian(PFcorrs_pooled.(epochNames{iEpoch})(LTstability_pooled.(epochNames{iEpoch}) > thresh_z))];

    
    [xData, yData] = myViolin_oneSided(PFcorrs_pooled.(epochNames{iEpoch})(LTstability_pooled.(epochNames{iEpoch}) < thresh_z)); 
    h2 = patch(tt-5*xData, yData, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    tt = tt+1;
    
end

xlim([0.5 1.5])


% legend([h1, h2], 'stable LTs', 'nonstable LTs', 'fontsize', fontsize, 'location', 'northout', 'box', 'off')

% xlabel('epoch', 'fontsize', fontsize)
ylabel('LT correlation bw pre and post(z)', 'fontsize', fontsize)

set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'xtick', [1 2], 'xticklabel', {'pre'; 'post'}, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 





%% correlation between LT consistency during an epoch and its correlation with PFs


plotwidth  = 200;
plotheight = 250;
fontsize   = 6;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



nor = 2;
noc = 1;

leftmargin = 30;  rightmargin = 50;    topmargin = 60;    bottommargin = 60;

gapc = 0;
gapr = 5;

sub_pos = subplot_pos(plotwidth, plotheight, leftmargin, rightmargin, bottommargin, topmargin, noc, nor, gapr, gapc);


nSessions = size(epochsXCorr_zscore, 1);


clear ax 

tt = 1;
for iEpoch = [1 3]
    
    currEpoch = epochNames{iEpoch};
    
    ax(tt) = axes('position',sub_pos{tt},'XGrid','off','XMinorGrid','off','FontSize',6,'Box','off','Layer','top');

    hold on
        
    LTstability     = cell(nSessions, 1); % LT consistency
    PFcorrs         = cell(nSessions, 1);
    PFcorrs_zscored = cell(nSessions, 1);
    
    for iSess = 1:nSessions
        
        okUnits = intersect(intersect(activeUnits{iSess}.pre, activeUnits{iSess}.post), activeUnits{iSess}.run);
        
        LTstability{iSess} = squeeze(epochsXCorr_zscore{iSess}(iEpoch, iEpoch, okUnits));
        PFcorrs{iSess} = assemblyTuningPFcorr_pooled{iSess}.(currEpoch)(okUnits);
        PFcorrs_zscored{iSess} = assemblyTuningPFcorr_zscore_pooled{iSess}.(currEpoch)(okUnits);
    end
    
    LTstability_pooled.(currEpoch)          = cell2mat(LTstability);
    PFcorrs_pooled.(currEpoch)         = cell2mat(PFcorrs);
    PFcorrs_zscored_pooled.(currEpoch) = cell2mat(PFcorrs_zscored);
    
    
    consisUnits = find(abs(LTstability_pooled.(currEpoch)) > 2);
    
    
%     allCorr = PFcorrs_pooled.(currEpoch)(consisUnits);
    signIdx = abs(PFcorrs_zscored_pooled.(currEpoch)(consisUnits)) > 1 & PFcorrs_zscored_pooled.(currEpoch)(consisUnits).*PFcorrs_pooled.(currEpoch)(consisUnits) > 0; 
    
    
    binEdges   = -1:0.05:1;
    binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;
    
    allCounts    = histc(PFcorrs_pooled.(currEpoch)(consisUnits), binEdges); allCounts(end) = [];
    nonSigCounts = histc(PFcorrs_pooled.(currEpoch)(consisUnits(~signIdx)), binEdges); nonSigCounts(end) = [];
    

    bar(binCenters, allCounts, 'FaceColor', 'k', 'FaceAlpha', 0.8, 'EdgeColor', 'k')
    bar(binCenters, nonSigCounts, 'FaceColor', 'w', 'EdgeColor', 'k')
    
    
    set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', fontsize, 'TickDir', 'out','TickLength',[0.01, 0.01] ) 
    
    tt = tt + 1;

end

% 
% xlim([0.5 2.5])
% 
% 
% xlabel('learned tuning consistency', 'fontsize', fontsize)
% ylabel('PF correlation', 'fontsize', fontsize)




    
%%
function output = nanprctile(data, percentile)

data = data(~isnan(data));
output = prctile(data, percentile);

end


function KLdiv = calKLDivergence(learnedTunings, spatialTunings)


% the individual learned tunnigs and spatial tunings should be in columns


% this function works well when when the learnedTunings matrix contains multipels instance of learned tunings (for examples if they belong to different time windows or different instances of unit identity shuffle)

nPosBins = size(spatialTunings, 1);

spatialTunings = spatialTunings ./ repmat(sum(spatialTunings, 1), [nPosBins 1]);
spatialTunings = spatialTunings + eps;

learnedTunings = learnedTunings ./ repmat(sum(learnedTunings, 1), [nPosBins 1]);
learnedTunings = learnedTunings + eps;


nSTs = size(spatialTunings, 2);
nLTs = size(learnedTunings, 2); % number of learned tunings


KLdiv = nan(nLTs, nSTs);

for ii = 1:nSTs
    
    currSpatialTuning = spatialTunings(:, ii);

    spatialTuningTerm = repmat(currSpatialTuning, [1 nLTs]);
    KLdiv(:, ii) = sum(spatialTuningTerm .* (log(spatialTuningTerm) - log(learnedTunings)), 1);

end

end


function [xData, yData] = myViolin(data)

binEdges = linspace(min(data)-0.1*range(data), max(data)+0.1*range(data), 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
count = count';

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

xData = [count -fliplr(count)];
yData = [binCenters fliplr(binCenters)];

end


function [xData, yData] = myViolin_oneSided(data)

binEdges = linspace(min(data)-0.1*range(data), max(data)+0.1*range(data), 50);

count = histc(data, binEdges); 
count(end) = [];
binCenters = binEdges(1:end-1) + diff(binEdges(1:2))/2;


count = count/sum(count);
count = count';

gw = gausswindow(2,6);
count = conv(count, gw, 'same');

xData = [count zeros(1, numel(count))];
yData = [binCenters fliplr(binCenters)];

end
