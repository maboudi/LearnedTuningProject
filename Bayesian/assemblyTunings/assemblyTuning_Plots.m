
parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/sessions_calculated_PBEinfo/';


rr = dir(parentDir);
sessionName = rr(sessionNumber+2).name;

basePath = fullfile(parentDir, sessionName);


load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'spikes')
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_allPBEs.mat']))
load(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTuning_vs_replayScores.mat']))


% for better plots sort the units based on their place field peak position
nUnits = numel(spikes);
nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
spatialTunings_merge = zeros(nUnits, nPosBins);

peakPosBins = zeros(nUnits, 1);
for iUnit = 1:nUnits
    peakPosBins(iUnit) = spikes(iUnit).peakPosBin.uni;
    spatialTunings_merge(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
end
    
[~, pf_sortIndices] = sort(peakPosBins, 'descend');


%% plot assembly tunings, distributions of their spatial information and correlation with corresponding place fields

figure;
periodNames = {'pre', 'run', 'post'};


% % actual place fields

subplot(5,3,4)

spatialTunings_merge = spatialTunings_merge./repmat(max(spatialTunings_merge, [], 2), [1, nPosBins]);
sortedPFs            = spatialTunings_merge(pf_sortIndices, :);

imagesc(sortedPFs); 

colormap(gca, 'jet'); 
xlabel('track position', 'fontsize', 12); 
ylabel('unit', 'fontsize', 12); 
title('Place Fields'); 
caxis([0 1])


subplot(5,3,5)
corr(sortedPFs', sortedPFs');
xlabel('unit place field', 'fontsize', 12)
ylabel('unit place field', 'fontsize', 12)
xlim([0 nActiveUnits])
ylim([0 nActiveUnits])
caxis([-1 1])
colormap(gca, 'jet');


% % assembly tunings in each behavioral epoch
for iperiod = 1:3
    
    currPeriod         = periodNames{iperiod};
    currAssemblyTuning = assemblyTunings.(currPeriod).data;
    currAssemblyTuning = currAssemblyTuning(pf_sortIndices, :);
    
    subplot(5,3, iperiod+6)
    
    imagesc(currAssemblyTuning ./ repmat(max(currAssemblyTuning, [], 2), [1, nPosBins]))
    colormap(gca, 'jet'); 
    xlabel('track position', 'fontsize', 12); 
    if iperiod == 1; ylabel('unit', 'fontsize', 12); end 
    title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
    caxis([0 1])
    
end


% % correlation matrix consisted of correlation values between assembly tunings and place fields

for iperiod = 1:3
    
    subplot(5, 3, iperiod+9)
    
    currPeriod = periodNames{iperiod};
    currCorrMat = assemblyTuningCorrMat.(currPeriod).data;
    currCorrMat = currCorrMat(pf_sortIndices, pf_sortIndices);
    
    imagesc(currCorrMat)
   
    xlabel('unit place field', 'fontsize', 12)
    if iperiod == 1; ylabel('unit assembly tuning', 'fontsize', 12); end
    xlim([0 nActiveUnits])
    ylim([0 nActiveUnits])
    caxis([-1 1])
    colormap(gca, 'jet');

end




% % correlation between assembly tuning and place field of each unit

subplot(5,3,13)

allCorrValues = [assemblyTuningPFcorr.PRE.(currDataType).data; assemblyTuningPFcorr.RUN.(currDataType).data; assemblyTuningPFcorr.POST.(currDataType).data];
bins = linspace(min(allCorrValues), max(allCorrValues)+0.01, 15);

colors = {'b'; 'none';'r'};
edgeColors = {'none'; 'k'; 'none'};
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    h = histc(assemblyTuningPFcorr.(currPeriod).(currDataType).data, bins); h(end) = [];
    a(iperiod ) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/nActiveUnits, 'Edgecolor', edgeColors{iperiod}, 'facecolor', colors{iperiod}, 'facealpha', 0.3);
    hold on
    
    rr = ylim;
    medianCorr = nanmedian(assemblyTuningPFcorr.(currPeriod).(currDataType).data);
    cl = colors{iperiod}; 
    if iperiod == 2; cl = 'k'; end
    line([medianCorr medianCorr], [0 rr(2)], 'color', cl, 'linewidth', 2)
    
end

legend(a, 'PRE', 'RUN', 'POST')
xlabel({'correlation b/w assembly tuning'; 'and place field of each unit'})
ylabel('ratio', 'fontsize', 12)
set(gca, 'box', 'off')




% % spatial information of the assembly tunings

subplot(5,3,14)


allSpatialInfo = [assemblyTuningSpatialInfo.PRE.(currDataType); assemblyTuningSpatialInfo.RUN.(currDataType); assemblyTuningSpatialInfo.POST.(currDataType)]; 
bins = linspace(min(allSpatialInfo), max(allSpatialInfo)+0.01, 50);


colors = {'b'; 'k';'r'};
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    h = histc(assemblyTuningSpatialInfo.(currPeriod).(currDataType).data, bins); h(end) = [];
    a(iperiod ) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/nActiveUnits, 'Edgecolor', edgeColors{iperiod}, 'facecolor', colors{iperiod}, 'facealpha', 0.3);

%     plot(bins(1:end-1)+diff(bins(1:2))/2, cumsum(h)/sum(h), colors{iperiod}, 'linewidth', 2)
    hold on
    
    rr = ylim;
    medianCorr = nanmedian(assemblyTuningSpatialInfo.(currPeriod).(currDataType).data);
    cl = colors{iperiod}; 
    if iperiod == 2; cl = 'k'; end
    line([medianCorr medianCorr], [0 rr(2)], 'color', cl, 'linewidth', 2)
    
end


legend(a, 'PRE', 'RUN', 'POST', 'place Fields') 
xlabel('spatial information', 'fontsize', 12)
ylabel('cumulative ratio', 'fontsize', 12)
set(gca, 'box', 'off')



subplot(17,3, 45)
scatter(assemblyTuningSpatialInfo.RUN.(currDataType).data, assemblyTuningPFcorr.RUN.(currDataType).data, 10, 'MarkerFaceColor', 'k', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')
xlim([0 prctile(assemblyTuningSpatialInfo.RUN.(currDataType).data, 90)])
legend('RUN')


subplot(17,3, [48 51])
scatter(assemblyTuningSpatialInfo.PRE.(currDataType).data, assemblyTuningPFcorr.PRE.(currDataType).data, 10, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')
hold on
scatter(assemblyTuningSpatialInfo.POST.(currDataType).data, assemblyTuningPFcorr.POST.(currDataType).data, 10, 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')

xlim([0 prctile([assemblyTuningSpatialInfo.PRE.(currDataType).data; assemblyTuningSpatialInfo.POST.(currDataType).data], 90)])
xlabel('spatial information', 'fontsize', 12)
ylabel('correlation with PFs')
legend('PRE', 'POST')







%%   Assembly tunings versus BD replay scores     
        
% plot figures

figure;
set(gcf, 'position', [3000 50 410 900])
colors = {'b'; 'k'; 'r'};
edgeColors = {'k', 'g', 'k'};
panelNumbers = [3 4 5; 7 8 9; 11 12 13];


% spatial information

for iperiod = 1:3

    currPeriod = periodNames{iperiod};

    subplot(14,1, panelNumbers(iperiod, :))

    h = violinplot(asTuningSpatialInfo_sub.(currPeriod).data);

    for ii = 1:4

        h(ii).ViolinColor = colors{iperiod};
        h(ii).EdgeColor   = edgeColors{iperiod};

    end

    xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})

    ylabel('spatial information (bits)')

    if iperiod == 3; xlabel('Bayesian replay scores'); end
    xtickangle(45)

    if ismember(iperiod, [1 3])
        ylim([0 0.65])
    end


%     hold on

%     for ii = 1:4

%         dataMedian    = median(asTuningSpatialInfo_sub.(currPeriod).data(:, ii));
%         shuffleMedian = median(asTuningSpatialInfo_sub.(currPeriod).data(:, ii));
        
%         line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iperiod}, 'linestyle', ':', 'linewidth', 2)

%         pval = ranksum(asTuningSpatialInfo_sub.(currPeriod)(activeUnits, ii), spatialInfo_ui.(currPeriod)(activeUnits, ii), 'tail', 'right');
%         text(ii-0.5, dataMedian+0.2, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
%         text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriod)(ii)), 'color', 'b', 'fontsize', 8)
%     end
    
%     hold off

    tt = ylim;
    text(0.5, tt(2), currPeriod, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

    dim = [0.1 0.8 0.5 .1];
    annotation('textbox',dim,'String', {sessionName; sprintf('%s - %s', currShuffleMethod, currReplayMethod)}, ...
          'FitBoxToText','on', 'Interpreter', 'none');

end


% correlation of assembly tunings with place fields


figure;
set(gcf, 'position', [3000 50 410 900])
colors = {'b'; 'k'; 'r'};
edgeColors = {'k', 'g', 'k'};
panelNumbers = [3 4 5; 7 8 9; 11 12 13];

for iperiod = 1:3

    currPeriod = periodNames{iperiod};

    subplot(14,1, panelNumbers(iperiod, :))

    h = violinplot(asTuningPFcorr_sub.(currPeriod));

    for ii = 1:4

        h(ii).ViolinColor = colors{iperiod};
        h(ii).EdgeColor   = edgeColors{iperiod};

    end

    xticklabels({'0-25%'; '25-50%'; '50-75%'; '75-100%'})

    ylabel('correlation with PF')

    if iperiod == 3; xlabel('Bayesian replay scores'); end
    xtickangle(45)


    ylim([-1 1])


    hold on
    for ii = 1:4

%                 validUnits = ~isnan(asTuningPFcorr_sub.(currPeriod)(:, ii));
        dataMedian    = median(asTuningPFcorr_sub.(currPeriod)(activeUnits, ii));
        shuffleMedian = median(corrWithPFs_ui.(currPeriod)(activeUnits, ii));
        line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iperiod}, 'linestyle', ':', 'linewidth', 2)

        pval = ranksum(asTuningPFcorr_sub.(currPeriod)(activeUnits, ii), corrWithPFs_ui.(currPeriod)(activeUnits, ii), 'tail', 'right');
        text(ii-0.5, dataMedian+0.2, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
        text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriod)(ii)), 'color', 'b', 'fontsize', 8)
    end

    hold off

    tt = ylim;
    text(0.5, tt(2), currPeriod, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

    dim = [0.1 0.8 0.5 .1];
    annotation('textbox',dim,'String', {sessionName; sprintf('%s - %s', currShuffleMethod, currReplayMethod)}, ...
          'FitBoxToText','on', 'Interpreter', 'none');

end




%% spatial tunings calculated including only PBEs with high BD scores 


shuffleMethods       = {'wPBEtimeswap'; 'unitIDshuffle'};
replayScoringMethods = {'weightedCorr'; 'replayScore'};
currPeriodName       = 'POST';

replayScoreThresh    = 75;

numPBEsEachMethod = zeros(2, 2);
for ii = 1:2
    currShuffleMethod = shuffleMethods{ii};
    for jj = 1:2
        
        currRepScorMethod = replayScoringMethods{jj};
        
        currBDscores = BDseqscore.(currPeriodName).data.(currShuffleMethod).(currRepScorMethod).prctilescore;
        
        numPBEsEachMethod(ii, jj) = numel(find(currBDscores > replayScoreThresh));
        
        
    end
end


PBEsetSize = min(numPBEsEachMethod(:));

normAssemblyTunings_subset = [];

for ii = 1:2
    currShuffleMethod = shuffleMethods{ii};
    for jj = 1:2
        
        currRepScorMethod = replayScoringMethods{jj};
        
        currBDscores = BDseqscore.(currPeriodName).data.(currShuffleMethod).(currRepScorMethod).prctilescore;
        allPBEs = find(currBDscores > replayScoreThresh);
        
        currBDscores = allPBEs(randperm(numel(allPBEs), PBEsetSize));
        
        assemblyTunings.(currShuffleMethod).(currRepScorMethod)     = calculateAssemblyTuning(binnedPBEs.(currPeriodName), spatialTunings_A, spatialTunings_B, currBDscores);
        normAssemblyTunings_subset.(currShuffleMethod).(currRepScorMethod) = assemblyTunings.(currShuffleMethod).(currRepScorMethod) ./ repmat(max(assemblyTunings.(currShuffleMethod).(currRepScorMethod), [], 2), [1, nPosBins]);

    end
end



figure;
set(gcf, 'position', [3200 300 570 770])

for ii = 1:numel(shuffleMethods)
    
    currShuffleMethod = shuffleMethods{ii};
    for jj = 1:numel(replayScoringMethods)     
        currRepScorMethod = replayScoringMethods{jj};
        
        currPanel = ii*2 + jj;
        subplot(3, 2, currPanel)
        
        BDscores = BDseqscore.(currPeriodName).data.(currShuffleMethod).(currRepScorMethod).prctilescore;
        
        plotseqscoredist(BDscores, [], [], [], 0, sprintf('%s\n  (n=%d)', currPeriodName, length(binnedPBEs.(currPeriodName))))
        ylabel('ratio of PBEs', 'fontsize', 10)
        
        title([currRepScorMethod '--' currShuffleMethod])
        
    end
end

dim = [0.1 0.8 0.5 .1];
annotation('textbox', dim,'String', {sessionName; currPeriodName}, ...
      'FitBoxToText','on', 'Interpreter', 'none', 'fontsize', 16);



figure; 
set(gcf, 'position', [3200 300 570 770])

for ii = 1:2
    
    currShuffleMethod = shuffleMethods{ii};
    for jj = 1:2
        
        currPanel = ii*2 + jj;
        subplot(3, 2, currPanel)
        
        currRepScorMethod = replayScoringMethods{jj};
        
        currMat = normAssemblyTunings_subset.(currShuffleMethod).(currRepScorMethod)(sortInd(1:nActiveUnits), :);
        imagesc(currMat)
        xlabel('position bin')
        ylabel('PF-sorted units')
        
        title([currRepScorMethod '--' currShuffleMethod])
        if currPanel == 6; colorbar; end
    end
end
colormap('jet')

dim = [0.1 0.8 0.5 .1];
annotation('textbox', dim,'String', {sessionName; currPeriodName}, ...
      'FitBoxToText','on', 'Interpreter', 'none', 'fontsize', 16);
  
  
  
  




%% Supplementary analysis to check whether the number of coactive units with a given unit is a significant factor in the unit's spatial tuning 
 

% For each unit calculate whether there is a significant difference between
% any two recording blocks

nUnits = numel(nCoactiveUnits_wEachUnitRUN);

baseStruct1 = struct('mean', zeros(nUnits, 1));
baseStruct2 = struct('meanDiff', zeros(nUnits, 1), 'pval', zeros(nUnits, 1));
nCoactiveUnitsComaparsions = struct('PRE', baseStruct1, 'RUN', baseStruct1, 'POST', baseStruct1, 'POSTvsPRE', baseStruct2, 'RUNvsPRE', baseStruct2, 'RUNvsPOST', baseStruct2);


for unit = 1:nUnits
    
    PREvalues  = nCoactiveUnits_wEachUnitPRE{unit};
    RUNvalues  = nCoactiveUnits_wEachUnitRUN{unit};
    POSTvalues = nCoactiveUnits_wEachUnitPOST{unit};
    
    PREvalues  = PREvalues(PREvalues > 0); 
    RUNvalues  = RUNvalues(RUNvalues > 0);
    POSTvalues = POSTvalues(POSTvalues > 0);
    
    
    if ~isempty(PREvalues)
    
        nCoactiveUnitsComaparsions.PRE.mean(unit) = mean(PREvalues);
        nCoactiveUnitsComaparsions.RUN.mean(unit) = mean(RUNvalues);
        nCoactiveUnitsComaparsions.POST.mean(unit) = mean(POSTvalues);
        

        nCoactiveUnitsComaparsions.POSTvsPRE.meanDiff(unit) = nCoactiveUnitsComaparsions.POST.mean(unit) - nCoactiveUnitsComaparsions.PRE.mean(unit);
        nCoactiveUnitsComaparsions.POSTvsPRE.pval(unit)       = ranksum(POSTvalues, PREvalues);

        nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(unit)  = nCoactiveUnitsComaparsions.RUN.mean(unit) - nCoactiveUnitsComaparsions.PRE.mean(unit);
        nCoactiveUnitsComaparsions.RUNvsPRE.pval(unit)        = ranksum(RUNvalues, PREvalues);

        nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff(unit) = nCoactiveUnitsComaparsions.RUN.mean(unit) - nCoactiveUnitsComaparsions.POST.mean(unit);
        nCoactiveUnitsComaparsions.RUNvsPOST.pval(unit)       = ranksum(RUNvalues, POSTvalues);

    end  
    
end


%%


figure;
fontsize = 8;

% panel A

% first, compare number of co-active units within a time bin during each
% period

nCoactiveUnits_all = [nCoactiveUnitsPRE nCoactiveUnitsRUN nCoactiveUnitsPOST];
nCoactiveUnits_all = nCoactiveUnits_all(nCoactiveUnits_all > 0); % we exclude the silent bins because they are not included when we calculate spike-weighted average of posteriors

bins = unique(nCoactiveUnits_all);

PRECounts = hist(nCoactiveUnitsPRE, bins);
RUNCounts = hist(nCoactiveUnitsRUN, bins);
POSTCounts = hist(nCoactiveUnitsPOST, bins);



subplot(6,3, [4 5])

bar(bins, PRECounts/numel(nCoactiveUnitsPRE), 'Edgecolor', 'none', 'facecolor', 'b', 'facealpha', 0.3)
hold on
bar(bins, RUNCounts/numel(nCoactiveUnitsRUN), 'Edgecolor', 'k', 'facecolor', 'none', 'facealpha', 0.3)
hold on
bar(bins, POSTCounts/numel(nCoactiveUnitsPOST), 'Edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.3)

legend('PRE', 'RUN', 'POST')
xlabel('number of coactive units in a PBE time bin', 'fontsize', fontsize)
ylabel('ratio of time bins', 'fontsize', fontsize)
set(gca, 'box', 'off')

POST_PRE_pvalue = ranksum(nCoactiveUnitsPOST, nCoactiveUnitsPRE, 'tail', 'right');
RUN_PRE_pvalue  = ranksum(nCoactiveUnitsRUN, nCoactiveUnitsPRE, 'tail', 'right');
RUN_POST_pvalue = ranksum(nCoactiveUnitsRUN, nCoactiveUnitsPOST, 'tail', 'right');

hold on 
text(mean(bins), 0.1, sprintf('POST vs PRE pval = %.3f', POST_PRE_pvalue), 'fontsize', fontsize)
hold on
text(mean(bins), 0.075, sprintf('RUN vs PRE pval = %.3f', RUN_PRE_pvalue), 'fontsize', fontsize)
hold on
text(mean(bins), 0.05, sprintf('RUN vs POST pval = %.3f', RUN_POST_pvalue), 'fontsize', fontsize)



% panel B
% % second, see how average number of coactive unit with a unit is
% different when comparing any two periods

allValues = [nCoactiveUnitsComaparsions.POSTvsPRE.meanDiff; nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff; nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff];

bins       = linspace(min(allValues), 1.01*max(allValues), 30);
binCenters = bins(1:end-1) + diff(bins(1:2))/2;

% POST vs PRE
POSTvsPRE_sig = nCoactiveUnitsComaparsions.POSTvsPRE.meanDiff(nCoactiveUnitsComaparsions.POSTvsPRE.pval > 0 & nCoactiveUnitsComaparsions.POSTvsPRE.pval <= 0.01);
POSTvsPRE_ns  = nCoactiveUnitsComaparsions.POSTvsPRE.meanDiff(nCoactiveUnitsComaparsions.POSTvsPRE.pval > 0 & nCoactiveUnitsComaparsions.POSTvsPRE.pval > 0.01);

POSTvsPRE_sig_counts = histc(POSTvsPRE_sig, bins); POSTvsPRE_sig_counts(end) = [];
POSTvsPRE_ns_counts  = histc(POSTvsPRE_ns, bins); POSTvsPRE_ns_counts(end) = [];


% RUN vs PRE
RUNvsPRE_sig = nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(nCoactiveUnitsComaparsions.RUNvsPRE.pval > 0 & nCoactiveUnitsComaparsions.RUNvsPRE.pval <= 0.01);
RUNvsPRE_ns  = nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(nCoactiveUnitsComaparsions.RUNvsPRE.pval > 0 & nCoactiveUnitsComaparsions.RUNvsPRE.pval > 0.01);

RUNvsPRE_sig_counts = histc(RUNvsPRE_sig, bins); RUNvsPRE_sig_counts(end) = [];
RUNvsPRE_ns_counts  = histc(RUNvsPRE_ns, bins); RUNvsPRE_ns_counts(end) = [];



% RUN vs POST
RUNvsPOST_sig = nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff(nCoactiveUnitsComaparsions.RUNvsPOST.pval > 0 & nCoactiveUnitsComaparsions.RUNvsPOST.pval <= 0.01);
RUNvsPOST_ns  = nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff(nCoactiveUnitsComaparsions.RUNvsPOST.pval > 0 & nCoactiveUnitsComaparsions.RUNvsPOST.pval > 0.01);

RUNvsPOST_sig_counts = histc(RUNvsPOST_sig, bins); RUNvsPOST_sig_counts(end) = [];
RUNvsPOST_ns_counts  = histc(RUNvsPOST_ns, bins); RUNvsPOST_ns_counts(end) = [];



subplot(6,3,7)
h1 = bar(binCenters, [POSTvsPRE_ns_counts POSTvsPRE_sig_counts], 1, 'stacked', 'FaceAlpha', 0.6);
h1(1).FaceColor = [1 1 1];
h1(2).FaceColor = [.1 .1 .1];
title('POST-PRE', 'fontsize', fontsize)
ylabel('unit counts', 'fontsize', fontsize)
% xlabel({'\Deltamean number of co-active units' ; 'with each unit'}, 'fontsize', fontsize)



subplot(6,3,8)
h2 = bar(binCenters, [RUNvsPRE_ns_counts RUNvsPRE_sig_counts], 1, 'stacked', 'FaceAlpha', 0.6);
h2(1).FaceColor = [1 1 1];
h2(2).FaceColor = [.1 .1 .1];
title('RUN-PRE', 'fontsize', fontsize)
% ylabel('unit counts', 'fontsize', fontsize)
xlabel('\Deltamean number of co-active units with each unit', 'fontsize', fontsize)

ylim([0 20])
legend('ns','pval< 0.01')



subplot(6,3,9)
h3 = bar(binCenters, [RUNvsPOST_ns_counts RUNvsPOST_sig_counts], 1, 'stacked', 'FaceAlpha', 0.6);
h3(1).FaceColor = [1 1 1];
h3(2).FaceColor = [.1 .1 .1];
title('RUN-POST', 'fontsize', fontsize)
% ylabel('unit counts', 'fontsize', fontsize)
% xlabel({'\Deltamean number of co-active units' ; 'with each unit'}, 'fontsize', fontsize)

ylim([0 20])


% panel C
% % the difference with RUN in nCoactive units in less during POST compared
% with PRE

subplot(6,3,6)

includedUnits = find(nCoactiveUnitsComaparsions.PRE.mean > 0);

plot(nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(includedUnits), nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff(includedUnits), '.k', 'markersize', 5)
hold on
plot(-4:2, -4:2, 'color', [0.5 0.5 0.5])
hold on
a = polyfit(nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(includedUnits), nCoactiveUnitsComaparsions.RUNvsPOST.meanDiff(includedUnits), 1);
plot(nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(includedUnits), nCoactiveUnitsComaparsions.RUNvsPRE.meanDiff(includedUnits)*a(1) + a(2), '-r')

xlabel({'\Delta mean no. coactive units'; 'RUN-PRE'}, 'fontsize', fontsize)
ylabel({'\Delta mean no. coactive units'; 'RUN-POST'}, 'fontsize', fontsize)

legend({'units'; 'diagonal'; 'best fit'})


% for the following panels

sigPOSTminusPREunits = (nCoactiveUnitsComaparsions.POSTvsPRE.pval > 0 & nCoactiveUnitsComaparsions.POSTvsPRE.pval <= 0.01);

% panel D

subplot(6,3,10)
plotWithfitLine(nCoactiveUnitsComaparsions.PRE.mean, assemblyTuningSpatialInfo.PRE, includedUnits, sigPOSTminusPREunits, 1, [], 'spatial info (bits)', 'PRE', fontsize)

legend('sig POST-PRE no. coactive units', 'ns POST-PRE no. coactive units', 'best fit', 'median all units', 'median only ns units')


subplot(6,3,11)
plotWithfitLine(nCoactiveUnitsComaparsions.RUN.mean, spatialInfoRUN, includedUnits, sigPOSTminusPREunits, 1, 'no. coactive units with each unit', [], 'RUN', fontsize)

subplot(6,3,12)
plotWithfitLine(nCoactiveUnitsComaparsions.POST.mean, spatialInfoPOST, includedUnits, sigPOSTminusPREunits, 1, [], [], 'POST', fontsize)



% panel E

ax1 = subplot(6,3,13);
plotWithfitLine(nCoactiveUnitsComaparsions.PRE.mean, assemblyTuningPFcorr.PRE.(dataType), includedUnits, sigPOSTminusPREunits, 1, [], 'correlation with PF', 'PRE', fontsize)

ax2 = subplot(6,3,14);
plotWithfitLine(nCoactiveUnitsComaparsions.RUN.mean, assemblyTuningPFcorr.RUN.(dataType), includedUnits, sigPOSTminusPREunits, 1, 'no. coactive units with each unit', [], 'RUN', fontsize)

ax3 = subplot(6,3,15);
plotWithfitLine(nCoactiveUnitsComaparsions.POST.mean, assemblyTuningPFcorr.POST.(dataType), includedUnits, sigPOSTminusPREunits, 1, [], [], 'POST', fontsize)

linkaxes([ax1, ax2, ax3], 'y')

% panel F

ax1 = subplot(6,3,16);
plotWithfitLine(assemblyTuningSpatialInfo.PRE, assemblyTuningPFcorr.PRE.(dataType), includedUnits, sigPOSTminusPREunits, 0, [], 'correlation with PF', 'PRE', fontsize)

ax2 = subplot(6,3,17);
plotWithfitLine(spatialInfoRUN, assemblyTuningPFcorr.RUN.(dataType), includedUnits, sigPOSTminusPREunits, 0,'spatial info (bits)', [], 'RUN', fontsize)

ax3 = subplot(6,3,18);
plotWithfitLine(spatialInfoPOST, assemblyTuningPFcorr.POST.(dataType), includedUnits, sigPOSTminusPREunits, 0, [], [], 'POST', fontsize)

linkaxes([ax1, ax2, ax3], 'y')




%% functions

function plotWithfitLine(xVar, yVar, inlcudedUnits, sigDiffUnits, plotyMedians, xl, yl, ttl, fontsize)

% units with significant POST vs PRE difference in number of coactive units 
scatter(xVar(sigDiffUnits), yVar(sigDiffUnits), 6, 'MarkerFaceColor', 'k')

hold on

% rest of the units
scatter(xVar(setdiff(inlcudedUnits, sigDiffUnits)), yVar(setdiff(inlcudedUnits, sigDiffUnits)), 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');

hold on

a = polyfit(xVar(inlcudedUnits), yVar(inlcudedUnits), 1);
[rho, pval]       = corr(xVar(inlcudedUnits), yVar(inlcudedUnits));
[rho_ns, pval_ns] = corr(xVar(setdiff(inlcudedUnits, sigDiffUnits)), yVar(setdiff(inlcudedUnits, sigDiffUnits)));


plot(xVar(inlcudedUnits), xVar(inlcudedUnits)*a(1) + a(2), '-r')
hold on

if pval > 0.001
    text(mean(xVar(inlcudedUnits)), mean(yVar(inlcudedUnits)), sprintf('corr = %.2f(%.3f)', rho, pval), 'fontsize', fontsize)
else
    text(mean(xVar(inlcudedUnits)), mean(yVar(inlcudedUnits)), sprintf('corr = %.2f(<%.3f)', rho, 0.001), 'fontsize', fontsize)
end

if pval_ns > 0.001
    text(mean(xVar(inlcudedUnits)), mean(yVar(inlcudedUnits))+0.5, sprintf('corr-ns = %.2f(%.3f)', rho_ns, pval_ns), 'fontsize', fontsize, 'color', 'b')
else
    text(mean(xVar(inlcudedUnits)), mean(yVar(inlcudedUnits))+0.5, sprintf('corr-ns = %.2f(<%.3f)', rho_ns, 0.001), 'fontsize', fontsize, 'color', 'b')
end

if plotyMedians
   
    % all units
    
    ymedian = mean(yVar(inlcudedUnits));
    hold on
    xrange = xlim;
    
    line([xrange(1) xrange(2)], [ymedian ymedian], 'color', 'k')
    
    % all units excluding the units with significant difference in number
    % of coactive units
    
    ymedian = mean(yVar(setdiff(inlcudedUnits, sigDiffUnits)));
    hold on
    xrange = xlim;
    
    line([xrange(1) xrange(2)], [ymedian ymedian], 'color', 'b')   
    
end



if ~isempty(yl); ylabel(yl, 'fontsize', fontsize); end
if ~isempty(ttl); title(ttl, 'fontsize', fontsize); end
if ~isempty(xl); xlabel(xl, 'fontsize', fontsize); end

end