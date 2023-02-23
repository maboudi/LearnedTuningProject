


%% Replay evaluation using Bayesian decoding
% 
% % spatial tunings unit ID shuffle DONE FOR EACH PBE SEPARATELY
% % 
% % nShuffles = 200;
% % shuffleInclusionFRthresh = 1; % shuffling of place fields is done only between the units with firing rates above a certain threshold
% % 
% 
% spatialTunings_RL(~spatialTunings_RL) = 1e-4;
% spatialTunings_LR(~spatialTunings_LR) = 1e-4;
% 
% 
% BDinitialize = struct('PRE', [], 'RUN', [], 'POST', []);
% BDinitialize = structfun(@(x)(struct('data', [], 'pts', [], 'p', [], 'ts', [])), BDinitialize,'UniformOutput', false);
% 
% posteriorProbMatrix = BDinitialize;
% 
% posteriorProbMatrix_nn = BDinitialize;
% 
% weightedCorr        = BDinitialize;
% weightedDistCorr    = BDinitialize;
% BDseqscore          = BDinitialize;
% 
% 
% %%% PBEs
% 
% binDur = 0.02;
% 
% % PRE
% 
% % [posteriorProbMatrix.PRE.data, weightedCorr.PRE.data, weightedDistCorr.PRE.data, BDseqscore.PRE.data, weightedCorr_null.PRE.data, weightedDistCorr_null.PRE.data] = BDreplayDetect_weightedDist(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% [posteriorProbMatrix.PRE.data, posteriorProbMatrix_nn.PRE.data, weightedCorr.PRE.data, weightedDistCorr.PRE.data, BDseqscore.PRE.data, weightedCorr_null.PRE.data, weightedDistCorr_null.PRE.data] = BDreplayDetect_weightedDist(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % posteriorProbMatrix_nnPRE_p =  BDreplayDetect_weightedDist(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% % RUN
% [posteriorProbMatrix.RUN.data, posteriorProbMatrix_nn.RUN.data, weightedCorr.RUN.data, weightedDistCorr.RUN.data, BDseqscore.RUN.data, weightedCorr_null.RUN.data, weightedDistCorr_null.RUN.data] = BDreplayDetect_weightedDist(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % posteriorProbMatrix_nnRUN_p =  BDreplayDetect_weightedDist(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% 
% % POST
% [posteriorProbMatrix.POST.data, posteriorProbMatrix_nn.POST.data, weightedCorr.POST.data, weightedDistCorr.POST.data, BDseqscore.POST.data, weightedCorr_null.POST.data, weightedDistCorr_null.POST.data] = BDreplayDetect_weightedDist(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % posteriorProbMatrix_nnPOST_p =  BDreplayDetect_weightedDist(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);



% BDinitialize = struct('PRE', [], 'RUN', [], 'POST', []);
% BDinitialize = structfun(@(x)(struct('data', [], 'pts', [], 'p', [], 'ts', [])), BDinitialize,'UniformOutput', false);
% 
% posteriorProbMatrix = BDinitialize;
% replayScore         = BDinitialize;
% weightedCorr        = BDinitialize;
% jumpDist            = BDinitialize;
% BDseqscore          = BDinitialize;
% sigMatrix           = BDinitialize;
% onlyLineElements    = BDinitialize;
% coveredLen          = BDinitialize;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% two place field templates %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% PBEs
% 
% binDur = 0.02;
% 
% % PRE
% % PREsubset = randperm(length(PREbinnedPBEs.data), min(1000, length(PREbinnedPBEs.data)));
% 
% [posteriorProbMatrix.PRE.data, replayScore.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
%     onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, replayScore_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data] = BDreplayDetect_vFeb20(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% % [posteriorProbMatrix.PRE.p, replayScore.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p, ...
% %     onlyLineElements.PRE.p, begPosition.PRE.p, endPosition.PRE.p, replayScore_null.PRE.p, weightedCorr_null.PRE.p, coveredLen.PRE.p, coveredLen_null.PRE.p] = BDreplayDetect_vFeb20(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % 
% 
% 
% 
% % RUN
% % RUNsubset = randperm(length(RUNbinnedPBEs.data), min(1000, length(RUNbinnedPBEs.data)));
% 
% [posteriorProbMatrix.RUN.data, replayScore.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
%     onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, replayScore_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data] = BDreplayDetect_vFeb20(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% % [posteriorProbMatrix.RUN.p, replayScore.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p, ...
% %     onlyLineElements.RUN.p, begPosition.RUN.p, endPosition.RUN.p, replayScore_null.RUN.p, weightedCorr_null.RUN.p, coveredLen.RUN.p, coveredLen_null.RUN.p] = BDreplayDetect_vFeb20(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% 
% 
% % POST
% % POSTsubset = randperm(length(POSTbinnedPBEs.data), min(1000, length(POSTbinnedPBEs.data)));
% 
% [posteriorProbMatrix.POST.data, replayScore.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
%     onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, replayScore_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data] = BDreplayDetect_vFeb20(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% % [posteriorProbMatrix.POST.p, replayScore.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p, ...
% %     onlyLineElements.POST.p, begPosition.POST.p, endPosition.POST.p, replayScore_null.POST.p, weightedCorr_null.POST.p, coveredLen.POST.p, coveredLen_null.POST.p] = BDreplayDetect_vFeb20(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% one place field template %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%% PBEs
% 
% binDur = 0.02;
% 
% % PRE
% % PREsubset = randperm(length(PREbinnedPBEs.data), min(1000, length(PREbinnedPBEs.data)));
% 
% [posteriorProbMatrix.PRE.data, replayScore.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
%     onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, replayScore_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data] = BDreplayDetect_1D_vFeb20(PREbinnedPBEs.data, spatialTunings, [], binDur);
% 
% % [posteriorProbMatrix.PRE.p, replayScore.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p, ...
% %     onlyLineElements.PRE.p, begPosition.PRE.p, endPosition.PRE.p, replayScore_null.PRE.p, weightedCorr_null.PRE.p, coveredLen.PRE.p, coveredLen_null.PRE.p] = BDreplayDetect_vFeb20(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % 
% 
% 
% 
% % RUN
% % RUNsubset = randperm(length(RUNbinnedPBEs.data), min(1000, length(RUNbinnedPBEs.data)));
% 
% [posteriorProbMatrix.RUN.data, replayScore.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
%     onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, replayScore_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data] = BDreplayDetect_1D_vFeb20(RUNbinnedPBEs.data, spatialTunings, [], binDur);
% 
% % [posteriorProbMatrix.RUN.p, replayScore.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p, ...
% %     onlyLineElements.RUN.p, begPosition.RUN.p, endPosition.RUN.p, replayScore_null.RUN.p, weightedCorr_null.RUN.p, coveredLen.RUN.p, coveredLen_null.RUN.p] = BDreplayDetect_vFeb20(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% 
% 
% % POST
% % POSTsubset = randperm(length(POSTbinnedPBEs.data), min(1000, length(POSTbinnedPBEs.data)));
% 
% [posteriorProbMatrix.POST.data, replayScore.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
%     onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, replayScore_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data] = BDreplayDetect_1D_vFeb20(POSTbinnedPBEs.data, spatialTunings, [], binDur);
% 
% % [posteriorProbMatrix.POST.p, replayScore.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p, ...
% %     onlyLineElements.POST.p, begPosition.POST.p, endPosition.POST.p, replayScore_null.POST.p, weightedCorr_null.POST.p, coveredLen.POST.p, coveredLen_null.POST.p] = BDreplayDetect_vFeb20(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% % 



%% How firings of each neuron is tuned to the position decoded by other neurons within a time bin


nUnits = 52;
nPosBins = 159;

spatialTunings_merge = (spatialTunings_LR + spatialTunings_RL)/2;
% spatialTunings_merge = spatialTunings;

spatialTunings_merge = spatialTunings_merge./repmat(max(spatialTunings_merge, [], 2), [1, nPosBins]);

[~, peakPos] = max(spatialTunings_merge, [], 2);
[~, sortInd] = sort(peakPos, 'descend');


nActiveUnits = 48;
sortInd = sortInd(1:nActiveUnits);

% spatialTunings_RL = spatialTunings;
% spatialTunings_LR = [];

%%

clear nPBEs

% % PRE

% currBDscore = BDseqscore.PRE.data.unitIDshuffle.replayScore.prctilescore;
% selectPBEs = find(currBDscore > 75 && coveredLen.PRE.data > 0.05);

selectPBEs = find(secondaryPBEs_PRE(:, 6) == 1 & secondaryPBEs_PRE(:, 5) == 1); % Incluing only PBEs within NREM

[assemblyTuningsPRE, assemblyTuningsPRE_ui, spikeCountProbsPRE, pxPRE, nCoactiveUnits_wEachUnitPRE, nCoactiveUnitsPRE, nPBEs.PRE] = calculateAssemblyTuning(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, selectPBEs);

normAssemblyTunings.PRE.data = assemblyTuningsPRE.integrated ./ repmat(max(assemblyTuningsPRE.integrated, [], 2), [1, nPosBins]);
normAssemblyTunings.PRE.ui   = assemblyTuningsPRE_ui.integrated ./ repmat(max(assemblyTuningsPRE_ui.integrated, [], 2), [1, nPosBins]);



% % RUN

% currBDscore = BDseqscore.RUN.data.unitIDshuffle.replayScore.prctilescore;
% selectPBEs = find(currBDscore > 75 && coveredLen.RUN.data > 0.05);

selectPBEs = find(secondaryPBEs_RUN(:, 5) == 1); % PBEs overlapping with the ripples

[assemblyTuningsRUN, assemblyTuningsRUN_ui, spikeCountProbsRUN, pxRUN, nCoactiveUnits_wEachUnitRUN, nCoactiveUnitsRUN, nPBEs.RUN] = calculateAssemblyTuning(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, selectPBEs);

normAssemblyTunings.RUN.data = assemblyTuningsRUN.integrated ./ repmat(max(assemblyTuningsRUN.integrated, [], 2), [1, nPosBins]);
normAssemblyTunings.RUN.ui   = assemblyTuningsRUN_ui.integrated ./ repmat(max(assemblyTuningsRUN_ui.integrated, [], 2), [1, nPosBins]);




% % POST

% currBDscore = BDseqscore.POST.data.unitIDshuffle.replayScore.prctilescore;
% selectPBEs = find(currBDscore > 75 && coveredLen.POST.data > 0.05);

selectPBEs = find(secondaryPBEs_POST(:, 6) == 1 & secondaryPBEs_POST(:, 5) == 1); % Incluing only PBEs within NREM

[assemblyTuningsPOST, assemblyTuningsPOST_ui, spikeCountProbsPOST, pxPOST, nCoactiveUnits_wEachUnitPOST, nCoactiveUnitsPOST, nPBEs.POST] = calculateAssemblyTuning(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, selectPBEs);

normAssemblyTunings.POST.data = assemblyTuningsPOST.integrated ./ repmat(max(assemblyTuningsPOST.integrated, [], 2), [1, nPosBins]);
normAssemblyTunings.POST.ui   = assemblyTuningsPOST_ui.integrated ./ repmat(max(assemblyTuningsPOST_ui.integrated, [], 2), [1, nPosBins]);




%% How assembly tunigns varies as a function of BD replay score


periodNames = {'PRE', 'RUN', 'POST'};
normassemblyTunings = struct('PRE', [], 'RUN', [], 'POST', []); 


spatialInfo = struct('PRE', [], 'RUN', [], 'POST', []); 
corrWithPFs = struct('PRE', [], 'RUN', [], 'POST', []); 

spatialInfo_ui = struct('PRE', [], 'RUN', [], 'POST', []); 
corrWithPFs_ui = struct('PRE', [], 'RUN', [], 'POST', []); 


datasets = struct('PRE', [], 'RUN', [], 'POST', []);
datasets.PRE  = PREbinnedPBEs.data;
datasets.RUN  = RUNbinnedPBEs.data;
datasets.POST = POSTbinnedPBEs.data;

numPBEsinEachQuartile = struct('PRE', [], 'RUN', [], 'POST', []);


shuffleMethods       = {'wPBEtimeswap'; 'unitIDshuffle'};
replayScoringMethods = {'weightedCorr'; 'replayScore'};


secondaryPBEs2.PRE = secondaryPBEs_PRE;
secondaryPBEs2.RUN = secondaryPBEs_RUN;
secondaryPBEs2.POST = secondaryPBEs_POST;

for imethod = 1:2
    for ishuffle = 1:2


        currShuffleMethod = shuffleMethods{ishuffle};
        currReplayMethod  = replayScoringMethods{imethod};

%         for iperiod = 1:3
% 
%             currPeriodName = periodNames{iperiod};
% 
%             currBDscores = BDseqscore.(currPeriodName).data.(currShuffleMethod).(currReplayMethod).prctilescore;
% 
%             numPBEsinEachQuartile.(currPeriodName) = zeros(4,1);
% 
%             for ii = 1:4
%                 numPBEsinEachQuartile.(currPeriodName)(ii) = numel(find(currBDscores >= (ii-1)*25 & currBDscores <= ii*25));
%             end
% 
%         end
% 
% 
%         minPBEnumbers_PREandPOST = min([numPBEsinEachQuartile.PRE; numPBEsinEachQuartile.POST]);
%         minPBEnumbers_RUN        = min([numPBEsinEachQuartile.RUN]);

        for iperiod = 1:3

            currPeriodName = periodNames{iperiod};
            
            nPBEs_subset.(currPeriodName) = zeros(4,1);
            
            spatialInfo.(currPeriodName)    = zeros(nUnits, 4);
            corrWithPFs.(currPeriodName)    = zeros(nUnits, 4);

            spatialInfo_ui.(currPeriodName) = zeros(nUnits, 4);
            corrWithPFs_ui.(currPeriodName) = zeros(nUnits, 4);

            currBDscores = BDseqscore.(currPeriodName).data.(currShuffleMethod).(currReplayMethod).prctilescore;
            for ii = 1:4
                
                if iperiod == 2
                    currPBEs = find(currBDscores >= (ii-1)*25 & currBDscores <= ii*25 & secondaryPBEs2.(currPeriodName)(:, 5) == 1);
                else
                    currPBEs = find(currBDscores >= (ii-1)*25 & currBDscores <= ii*25 & secondaryPBEs2.(currPeriodName)(:, 6) == 1);
                end
                
                nPBEs_subset.(currPeriodName)(ii) = numel(currPBEs);
%                 [~, sortInd] = sort(currBDscores(currPBEs), 'descend');
% 
%                 if iperiod == 2
%                    currPBEs_equal = currPBEs(sortInd(1: minPBEnumbers_RUN));
%                 else
%                    currPBEs_equal = currPBEs(sortInd(1: minPBEnumbers_PREandPOST));
%                 end

%                 numel(currPBEs_equal)
                [assemblyTunings, assemblyTunings_ui] = calculateAssemblyTuning(datasets.(currPeriodName),  spatialTunings_RL, spatialTunings_LR, currPBEs);

                normAssemblyTunings_subset    = assemblyTunings.integrated ./ repmat(max(assemblyTunings.integrated, [], 2), [1, nPosBins]);
                normAssemblyTunings_subset_ui = assemblyTunings_ui.integrated ./ repmat(max(assemblyTunings_ui.integrated, [], 2), [1, nPosBins]);

                occupancy = mean(cell2mat(posteriorProbMatrix.(currPeriodName).data'), 2);
                occupancy = occupancy / sum(occupancy);


                % actual assembly tunings
                overallMeanTuning = sum(repmat(occupancy', [nUnits 1]).* normAssemblyTunings_subset, 2) ./ sum(occupancy);
                spatialInfo.(currPeriodName)(:, ii) = sum(repmat(occupancy', [nUnits 1]).* (normAssemblyTunings_subset ./ repmat(overallMeanTuning, [1 nPosBins])) .* log2(normAssemblyTunings_subset ./ repmat(overallMeanTuning, [1 nPosBins])), 2);

                allCorrs = corr(permute(normAssemblyTunings_subset, [2 1]), spatialTunings_merge');
                corrWithPFs.(currPeriodName)(:, ii) = diag(allCorrs);
                
                

                % assembly tunings from unit-ID shuffled PBEs
                overallMeanTuning = sum(repmat(occupancy', [nUnits 1]).* normAssemblyTunings_subset_ui, 2) ./ sum(occupancy);
                spatialInfo_ui.(currPeriodName)(:, ii) = sum(repmat(occupancy', [nUnits 1]).* (normAssemblyTunings_subset_ui ./ repmat(overallMeanTuning, [1 nPosBins])) .* log2(normAssemblyTunings_subset_ui ./ repmat(overallMeanTuning, [1 nPosBins])), 2);

                allCorrs = corr(permute(normAssemblyTunings_subset_ui, [2 1]), spatialTunings_merge');
                corrWithPFs_ui.(currPeriodName)(:, ii) = diag(allCorrs);

            end

        end



        % spatial information

        figure;
        set(gcf, 'position', [3000 50 410 900])
        colors = {'b'; 'k'; 'r'};
        edgeColors = {'k', 'g', 'k'};
        panelNumbers = [3 4 5; 7 8 9; 11 12 13];


        for iperiod = 1:3

            currPeriodName = periodNames{iperiod};

            subplot(14,1, panelNumbers(iperiod, :))

            h = violinplot(spatialInfo.(currPeriodName));

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



            hold on

            for ii = 1:4

                validUnits = ~isnan(spatialInfo.(currPeriodName)(:, ii));

                dataMedian    = median(spatialInfo.(currPeriodName)(validUnits, ii));
                shuffleMedian = median(spatialInfo_ui.(currPeriodName)(validUnits, ii));
                line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iperiod}, 'linestyle', ':', 'linewidth', 2)

                pval = ranksum(spatialInfo.(currPeriodName)(validUnits, ii), spatialInfo_ui.(currPeriodName)(validUnits, ii), 'tail', 'right');
                text(ii-0.5, dataMedian+0.2, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
                text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriodName)(ii)), 'color', 'b', 'fontsize', 8)
            end
            hold off

            tt = ylim;
            text(0.5, tt(2), currPeriodName, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

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

            currPeriodName = periodNames{iperiod};

            subplot(14,1, panelNumbers(iperiod, :))

            h = violinplot(corrWithPFs.(currPeriodName));

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

                validUnits = ~isnan(corrWithPFs.(currPeriodName)(:, ii));
                dataMedian    = median(corrWithPFs.(currPeriodName)(validUnits, ii));
                shuffleMedian = median(corrWithPFs_ui.(currPeriodName)(validUnits, ii));
                line([ii-0.3 ii+0.3], [shuffleMedian shuffleMedian], 'color', colors{iperiod}, 'linestyle', ':', 'linewidth', 2)

                pval = ranksum(corrWithPFs.(currPeriodName)(validUnits, ii), corrWithPFs_ui.(currPeriodName)(validUnits, ii), 'tail', 'right');
                text(ii-0.5, dataMedian+0.2, sprintf('p = %.3f', pval), 'color', [0.1 0.1 0.1], 'fontsize', 8)
                text(ii-0.5, dataMedian+0.4, sprintf('n = %d', nPBEs_subset.(currPeriodName)(ii)), 'color', 'b', 'fontsize', 8)
            end

            hold off

            tt = ylim;
            text(0.5, tt(2), currPeriodName, 'fontsize', 12, 'fontweight', 'bold', 'horizontalAlignment', 'center')

            dim = [0.1 0.8 0.5 .1];
            annotation('textbox',dim,'String', {sessionName; sprintf('%s - %s', currShuffleMethod, currReplayMethod)}, ...
                  'FitBoxToText','on', 'Interpreter', 'none');

        end
    end
end



%% vicarious place fields calculated icluding only PBEs with high BD scores 


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
        
        assemblyTunings.(currShuffleMethod).(currRepScorMethod)     = calculateAssemblyTuning(datasets.(currPeriodName), spatialTunings_RL, spatialTunings_LR, currBDscores);
        normAssemblyTunings_subset.(currShuffleMethod).(currRepScorMethod) = assemblyTunings.(currShuffleMethod).(currRepScorMethod).integrated ./ repmat(max(assemblyTunings.(currShuffleMethod).(currRepScorMethod).integrated, [], 2), [1, nPosBins]);

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
        
        plotseqscoredist(BDscores, [], [], [], 0, sprintf('%s\n  (n=%d)', currPeriodName, length(datasets.(currPeriodName))))
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

%% difference of the spatial tunings from 5% confidence interval

figure; 
cm = redblue;

% PRE
subplot(2,3,4)
currMat = spatialTuningPRE_ci.integrated;
currMat = currMat(sortInd, :)./repmat(max(currMat(sortInd, :), [], 2), [1, nPosBins]);
imagesc(currMat);
% caxis([0 0.24])

colorbar
xlabel('position bin')
ylabel('units')

% RUN
subplot(2,3,5)
currMat = spatialTuningRUN_ci.integrated;
currMat = currMat(sortInd, :)./repmat(max(currMat(sortInd, :), [], 2), [1, nPosBins]);
imagesc(currMat);
% caxis([0 0.24])

colorbar
xlabel('position bin')


% POST
subplot(2,3,6)
currMat = spatialTuningPOST_ci.integrated;
currMat = currMat(sortInd, :)./repmat(max(currMat(sortInd, :), [], 2), [1, nPosBins]);
imagesc(currMat);
% caxis([0 0.24])

colorbar
xlabel('position bin')

colormap('jet')




%% spatial information and correlation with place fields of assesmbly tunings


dataType = 'data'; % data or unit-ID shuffle PBEs

periodNames = {'PRE'; 'RUN'; 'POST'};

for iperiod = 1:3
    currPeriod = periodNames{iperiod};
    
    
    % % correlation of the assembly tunings with the place fields
    assemblyTuningCorrMat.(currPeriod).(dataType) = corr(normAssemblyTunings.(currPeriod).(dataType)(sortInd(1:nActiveUnits), :)', spatialTunings_merge(sortInd(1:nActiveUnits), :)');
    assemblyTuningPFcorr.(currPeriod).(dataType)  = diag(assemblyTuningCorrMat.(currPeriod).(dataType));
    
    
    % % spatial information of the assembly tunings
    
    virtualOccupancy  = mean(cell2mat(posteriorProbMatrix.(currPeriod).data'), 2);
    
    meanFR = sum(repmat(virtualOccupancy', [nActiveUnits 1]).* normAssemblyTunings.(currPeriod).(dataType)(sortInd(1:nActiveUnits), :), 2) ./ sum(virtualOccupancy);
    assemblyTuningSpatialInfo.(currPeriod).(dataType) = sum(repmat(virtualOccupancy', [nActiveUnits 1]).* (normAssemblyTunings.(currPeriod).(dataType)(sortInd(1:nActiveUnits), :) ./ repmat(meanFR, [1 nPosBins])) .* log2(normAssemblyTunings.(currPeriod).(dataType)(sortInd(1:nActiveUnits), :) ./ repmat(meanFR, [1 nPosBins])), 2);

end


%% plot assembly tunings and their quality features


figure;

periodNames = {'PRE'; 'RUN'; 'POST'};


% % actual place fields
subplot(5,3,4)
imagesc(spatialTunings_merge(sortInd, :)./repmat(max(spatialTunings_merge(sortInd, :), [], 2), [1, nPosBins])); colormap(gca, 'jet'); 
xlabel('track position', 'fontsize', 12); 
ylabel('unit', 'fontsize', 12); 
title('Place Fields'); 

caxis([0 1])

% % assembly tunings in each behavioral epoch
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    subplot(5,3, iperiod+6)
    
    imagesc(normAssemblyTunings.(currPeriod).(dataType)(sortInd, :)./repmat(max(normAssemblyTunings.(currPeriod).(dataType)(sortInd, :), [], 2), [1, nPosBins]))
    colormap(gca, 'jet'); 
    xlabel('track position', 'fontsize', 12); 
    if iperiod == 1; ylabel('unit', 'fontsize', 12); end 
    title([currPeriod sprintf('(n=%d)', nPBEs.(currPeriod))]);
    caxis([0 1])
    
end


% % correlation matrix consisted of correlation values between assemblt tunings and place fields

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    subplot(5, 3, iperiod+9)
    
    imagesc(assemblyTuningCorrMat.(currPeriod).(dataType))
   
    xlabel('unit place field', 'fontsize', 12)
    if iperiod == 1; ylabel('unit assembly tuning', 'fontsize', 12); end
    xlim([0 nActiveUnits])
    ylim([0 nActiveUnits])
    caxis([-1 1])
    colormap(gca, 'jet');

end


% % correlation between assembly tuning and place field of each unit


subplot(5,3,13)

allCorrValues = [assemblyTuningPFcorr.PRE.(dataType); assemblyTuningPFcorr.RUN.(dataType); assemblyTuningPFcorr.POST.(dataType)];
bins = linspace(min(allCorrValues), max(allCorrValues)+0.01, 15);

colors = {'b'; 'none';'r'};
edgeColors = {'none'; 'k'; 'none'};
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    h = histc(assemblyTuningPFcorr.(currPeriod).(dataType), bins); h(end) = [];
    a(iperiod ) = bar(bins(1:end-1)+diff(bins(1:2))/2, h/nActiveUnits, 'Edgecolor', edgeColors{iperiod}, 'facecolor', colors{iperiod}, 'facealpha', 0.3);
    hold on
    
    rr = ylim;
    medianCorr = nanmedian(assemblyTuningPFcorr.(currPeriod).(dataType));
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


allSpatialInfo = [assemblyTuningSpatialInfo.PRE.(dataType); assemblyTuningSpatialInfo.RUN.(dataType); assemblyTuningSpatialInfo.POST.(dataType)]; 

% allSpatialInfo = [assemblyTuningSpatialInfo.PRE.(dataType); assemblyTuningSpatialInfo.RUN.(dataType); assemblyTuningSpatialInfo.POST.(dataType); spatialInfo_LR; spatialInfo_RL]; 
% allSpatialInfo = [assemblyTuningSpatialInfo.PRE; spatialInfoRUN; spatialInfoPOST; spatialInfo]; % 

bins = linspace(min(allSpatialInfo), max(allSpatialInfo)+0.01, 50);

% PFs  = hist([spatialInfo_RL; spatialInfo_LR], bins); 

colors = {'b'; 'k';'r'};
for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    h = histc(assemblyTuningSpatialInfo.(currPeriod).(dataType), bins); h(end) = [];
    
    plot(bins(1:end-1)+diff(bins(1:2))/2, cumsum(h)/sum(h), colors{iperiod}, 'linewidth', 2)
    hold on
end

% plot(bins+diff(bins(1:2))/2, cumsum(PFs)/numel(spatialInfo_RL)/2, ':k', 'linewidth', 2)

legend('PRE', 'RUN', 'POST', 'place Fields') 
xlabel('spatial information', 'fontsize', 12)
ylabel('cumulative ratio', 'fontsize', 12)
set(gca, 'box', 'off')



subplot(17,3, 45)
scatter(assemblyTuningSpatialInfo.RUN.(dataType), assemblyTuningPFcorr.RUN.(dataType), 10, 'MarkerFaceColor', 'k', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')
xlim([0 prctile(assemblyTuningSpatialInfo.RUN.(dataType), 90)])
legend('RUN')

subplot(17,3, [48 51])
scatter(assemblyTuningSpatialInfo.PRE.(dataType), assemblyTuningPFcorr.PRE.(dataType), 10, 'MarkerFaceColor', 'b', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')
hold on
scatter(assemblyTuningSpatialInfo.POST.(dataType), assemblyTuningPFcorr.POST.(dataType), 10, 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',.5, 'MarkerEdgeColor','none')

xlim([0 prctile([assemblyTuningSpatialInfo.PRE.(dataType); assemblyTuningSpatialInfo.POST.(dataType)], 90)])
xlabel('spatial information', 'fontsize', 12)
ylabel('correlation with PFs')
legend('PRE', 'POST')




%%
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

%% Plot the figures




%% (1) Example PBEs sorted based on their sequence scores 

close all
binDur = 0.02;


datasets.PRE = PREbinnedPBEs;
datasets.RUN = RUNbinnedPBEs;
datasets.POST = POSTbinnedPBEs;

periods = {'PRE';'RUN'; 'POST'};

PBELen = struct('PRE', [], 'RUN', [], 'POST', []);

for pp = 1:3
    
    currentPBEs = datasets.(periods{pp});
    nPBEs = size(currentPBEs, 1);
    
    PBELen.(periods{pp}) = zeros(nPBEs, 1);
    nFiringUnits.(periods{pp}) = zeros(nPBEs, 1);
    
    for pbe = 1:nPBEs
        
        PBELen.(periods{pp})(pbe) = size(currentPBEs{pbe, 2}, 2);
        
        nFiringUnits.(periods{pp})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));
        
    end
    
end


secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.PREEpoch(1, 2), :);
secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 2) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 1), :);
secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.POSTEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 2), :);


acceptedEvts.PRE = find(secondaryPBEs_1.PRE(:, 8) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5);
acceptedEvts.RUN = find(secondaryPBEs_1.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
acceptedEvts.POST = find(secondaryPBEs_1.POST(:, 8) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);



shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
BDscoreMethods = {'weightedCorr'; 'replayScore'};

% 
% datasets.PRE = PREbinnedPBEs.data;
% datasets.RUN = RUNbinnedPBEs.data;
% datasets.POST = POSTbinnedPBEs.data;



for pp = 1: numel(periods)
    for ii = 1: numel(shuffleMethods)
        for jj = 1: numel(BDscoreMethods)

            period = periods{pp};
            shuffleMethod = shuffleMethods{ii};
            BDscoreMethod = BDscoreMethods{jj};


%             directory    = fullfile(mainDir, 'Bayesian_combDir_weightedDistCorr1', period);
%             mkdir(directory)

            
            originalIndices = 1:numel(acceptedEvts.(period));
            

            score2sortPBEs = BDseqscore.(period).data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.(period));
            score2sortPBEs = max(score2sortPBEs, [], 2);
            
            
            
            % Indices of PBEs with highest sequence scores
%             
%             aboveThreshPBEs = find(score2sortPBEs > 90);
%             
%             if numel(aboveThreshPBEs) < 60
%                aboveThreshPBEs = find(score2sortPBEs > 70);
%             end
%             
%             PBEidx2plot = aboveThreshPBEs(randperm(numel(aboveThreshPBEs), min(60, numel(aboveThreshPBEs))));
        
            [~, idx]    = sort(score2sortPBEs, 'descend'); 
%             PBEidx2plot = idx(1:60);
            PBEidx2plot = acceptedEvts.(period)(idx(1:min(60, numel(idx))));

            % % plot individual events
            
            sessionName = 'Achilles_10252013';
            directory = ['/home/kouroshmaboudi/Documents/HMM_project/GreatLakesNov17-2020/Grosmark_reclustered/' sessionName '/BayesianDecoding'];
            
            textOnTop = sprintf('Session: %s\nPeriod: %s\nPBEs were sorted based on their scores using %s and %s', sessionName, period, shuffleMethod, BDscoreMethod);
%             textOnTop = sprintf('Session: %s\nPeriod: %s\nrandom PBEs with scores above 90 percentile were plotted. The scores were calculated using %s and %s', sessionName, period, shuffleMethod, BDscoreMethod);

            filename = [period '_' shuffleMethod '_' BDscoreMethod];
            
            tempPlotBD3(datasets.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, [], [], [], runTemplate_LR, runTemplate_RL, PBEidx2plot, binDur, directory, filename, textOnTop) % directory
%             tempPlotBD3_1D(datasets.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, [], [], runTemplate, PBEidx2plot, binDur, directory, filename, textOnTop) % directory

        end
    end
end
%%
directory    = fullfile(mainDir, 'Bayesian_combDir_weightedDistCorr');
save(fullfile(directory, 'replayDetectionResults'), 'posteriorProbMatrix', 'weightedCorr', 'weightedDistCorr', 'BDseqscore', 'weightedCorr_null', 'weightedDistCorr_null')


%% within-PBE time swap positive and place field unit ID shuffle negative PBEs

% PRE
% 
% 
% PBEidx2plot = find(BDseqscore.POST.data.weightedCorr.prctilescore(:, 1) >= thresh & BDseqscore2.POST.data.weightedCorr.prctilescore >= thresh);
% 
% textOnTop = sprintf('Session: %s\nPeriod: %s\nPBEs classified as sequences using both within-PBE time swap and place field unit ID shuffle', sessionName, 'PRE');
% tempPlotBD3(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore2.POST.data, weightedCorr.POST.data, begPosition.POST.data, endPosition.POST.data, runTemplate_LR, runTemplate_RL, PBEidx2plot(randperm(numel(PBEidx2plot), 60)), binDur, 'POST_wc_tsPpfP', textOnTop) % directory
% 



%% (2) distribution of sequence scores


% directory = fullfile(mainDir, 'Bayesian_combDir_weightedDistCorr', 'scoreDistributions');

if ~exist(directory, 'dir'); mkdir(directory); end

shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
BDscoreMethods = {'weightedCorr'; 'replayScore'};


for ii = 1:numel(shuffleMethods)
    for jj = 1:numel(BDscoreMethods)

        shuffleMethod = shuffleMethods{ii};
        BDscoreMethod = BDscoreMethods{jj};
        
        
        figure;
        set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


        % % PRE
        subplot(3,3,[4 7])
        plotseqscoredist(BDseqscore.PRE.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.PRE, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'PRE', length(acceptedEvts.PRE)))
        ylabel('ratio of PBEs', 'fontsize', 10)

        % % RUN
        subplot(3,3,[5 8])
        plotseqscoredist(BDseqscore.RUN.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.RUN, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'RUN', length(acceptedEvts.RUN)))

        % % POST
        subplot(3,3,[6 9])
        plotseqscoredist(BDseqscore.POST.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.POST, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'POST', length(acceptedEvts.POST))) 

        
        % legend
        legendSub = subplot(3,3,3);

        hold on

%         p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
        p_pts = plot(1, nan, 'color', 'r', 'linestyle', ':', 'linewidth', 2);
        p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

        h_legend = legend([p, p_pts],{'data', 'diagonal'}, 'Location', 'South');


        set(h_legend, 'fontsize', 8)
        legend boxoff 

        set(legendSub, 'Visible', 'off');

        filename = fullfile(directory, [BDscoreMethod '-' shuffleMethod]);
        
        print(gcf, filename, '-dpdf')
        
%         savepdf(gcf, filename, '-dsvg')
%         saveas(gcf, filename,  'epsc')

    end
end

close all

%% scatter plots for comparing scores between different shuffling methods
% (for a given BD replay measure) and different BD replay measures (for a given shuffling method)

fileinfo.name = sessionName;
plotScatterHistogram_RSvsWC(BDseqscore, acceptedEvts, fileinfo, directory)

plotScatterHistogram_TSvsPF(BDseqscore, acceptedEvts, fileinfo, directory)


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