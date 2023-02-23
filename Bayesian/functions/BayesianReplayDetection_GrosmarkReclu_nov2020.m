

%%% linearized positions

% Validate BD
% Decode the postions during animal's running and comapre it with the
% actual position of the animal. The less the error the more reliable are
% the templates to be used for replay detection during the PBEs


binDur      = 0.25; % duration of active RUN time bins
posBinSize  = 2; % in cm


directory = fullfile(storageDir, 'placeFields');


thetaPeriods = [];

% LR

laps_LR_idx = 2:2:length(laps);
laps_LR_idx = find(ismember(laps(:,3), laps_LR_idx));

laps_LR    = laps(laps_LR_idx, :);

individualAndPopulationSpatialInfo2(spikes_pyr, laps_LR, spatialTunings_LR, behavior, thetaPeriods, turningPeriods, speed, runSpeedThresh, posBinSize, binDur, 'LR', fileInfo, directory);




% RL

laps_RL_idx = 1:2:length(laps);
laps_RL_idx = find(ismember(laps(:,3), laps_RL_idx));

laps_RL    = laps(laps_RL_idx, :);


individualAndPopulationSpatialInfo2(spikeStruct, laps_RL, spatialTunings_RL, behavior, thetaPeriods, nan_or_turningPeriods, speed, runSpeedThresh, posBinSize, binDur, 'RL', fileinfo, directory);



%% Replay evaluation using Bayesian decoding

% sz = getenv('SLURM_CPUS_PER_TASK');
% p  = parpool('local', str2num(sz));

subfolder = fullfile(mainDir, 'BayesianDecoding');
mkdir(subfolder)


BDinitialize = struct('PRE', [], 'RUN', [], 'POST', []);
BDinitialize = structfun(@(x)(struct('data', [], 'pts', [], 'p', [], 'ts', [])), BDinitialize,'UniformOutput', false);

posteriorProbMatrix = BDinitialize;
RadonTF             = BDinitialize;
weightedCorr        = BDinitialize;
jumpDist            = BDinitialize;
BDseqscore          = BDinitialize;
sigMatrix           = BDinitialize;
onlyLineElements    = BDinitialize;
coveredLen          = BDinitialize;


% PREbinnedPBEs  = PREbinnedPBEs.data;
% RUNbinnedPBEs  = RUNbinnedPBEs.data;
% POSTbinnedPBEs = POSTbinnedPBEs.data;



if exist('spatialTunings', 'var') %ismember(sessionNumber, [2 5 8]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% one place field template %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% PBEs

binDur = 0.02;
% 
% % PRE
% 
% [posteriorProbMatrix.PRE.data, RadonTF.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
%     onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, RadonTF_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data, jumpDist_null.PRE.data] = BDreplayDetect_1D_vMarch21(PREbinnedPBEs, spatialTunings, [], binDur);
% 
% % [posteriorProbMatrix.PRE.p, RadonTF.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p, ...
% %     onlyLineElements.PRE.p, begPosition.PRE.p, endPosition.PRE.p, RadonTF_null.PRE.p, weightedCorr_null.PRE.p, coveredLen.PRE.p, coveredLen_null.PRE.p, jumpDist_null.PRE.p] = BDreplayDetect_1D_vMarch21(PREbinnedPBEs.p, spatialTunings, [], binDur);
% 
% save(fullfile(subfolder, ['Bayesian_1_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
%     'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
%     'jumpDist_null', 'spatialTunings', 'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate', 'spatialInfo', 'conslapsRatio', 'diffWithAvg', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')



% RUN

% [posteriorProbMatrix.RUN.data, RadonTF.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
%     onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, RadonTF_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data, jumpDist_null.RUN.data] = BDreplayDetect_1D_vMarch21(RUNbinnedPBEs, spatialTunings, [], binDur);

% [posteriorProbMatrix.RUN.p, RadonTF.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p, ...
%     onlyLineElements.RUN.p, begPosition.RUN.p, endPosition.RUN.p, RadonTF_null.RUN.p, weightedCorr_null.RUN.p, coveredLen.RUN.p, coveredLen_null.RUN.p, jumpDist_null.RUN.p] = BDreplayDetect_1D_vMarch21(RUNbinnedPBEs.p, spatialTunings, [], binDur);


% save(fullfile(subfolder, ['Bayesian_2_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
%     'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
%     'jumpDist_null', 'spatialTunings', 'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate', 'spatialInfo', 'conslapsRatio', 'diffWithAvg', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')


% POST

[posteriorProbMatrix.POST.data, RadonTF.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
    onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, RadonTF_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data, jumpDist_null.POST.data] = BDreplayDetect_1D_vMarch21(POSTbinnedPBEs(1:3000, :), spatialTunings, [], binDur);

% [posteriorProbMatrix.POST.p, RadonTF.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p, ...
%     onlyLineElements.POST.p, begPosition.POST.p, endPosition.POST.p, RadonTF_null.POST.p, weightedCorr_null.POST.p, coveredLen.POST.p, coveredLen_null.POST.p, jumpDist_null.POST.p] = BDreplayDetect_1D_vMarch21(POSTbinnedPBEs.p, spatialTunings, [], binDur);



save(fullfile(subfolder, ['Bayesian_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
    'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
    'jumpDist_null', 'spatialTunings', 'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate', 'spatialInfo', 'conslapsRatio', 'diffWithAvg', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')


else


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% two place field templates %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PBEs

binDur = 0.02;
% 
% % PRE
% 
% [posteriorProbMatrix.PRE.data, RadonTF.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
%     onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, RadonTF_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data, jumpDist_null.PRE.data] = BDreplayDetect_vMarch21(PREbinnedPBEs, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% % [posteriorProbMatrix.PRE.p, RadonTF.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p, ...
% %     onlyLineElements.PRE.p, begPosition.PRE.p, endPosition.PRE.p, RadonTF_null.PRE.p, weightedCorr_null.PRE.p, coveredLen.PRE.p, coveredLen_null.PRE.p, jumpDist_null.PRE.p] = BDreplayDetect_vMarch21(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
% 
% save(fullfile(subfolder, ['Bayesian_1_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
%     'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
%     'jumpDist_null', 'spatialTunings_RL', 'spatialTunings_LR',  'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate_LR', 'runTemplate_RL', 'conslapsRatio_LR', 'diffWithAvg_LR', 'conslapsRatio_RL', 'diffWithAvg_RL', ...
%     'spatialInfo_RL', 'spatialInfo_LR', 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')



% RUN

% [posteriorProbMatrix.RUN.data, RadonTF.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
%     onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, RadonTF_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data, jumpDist_null.RUN.data] = BDreplayDetect_vMarch21(RUNbinnedPBEs, spatialTunings_RL, spatialTunings_LR, [], binDur);

% [posteriorProbMatrix.RUN.p, RadonTF.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p, ...
%     onlyLineElements.RUN.p, begPosition.RUN.p, endPosition.RUN.p, RadonTF_null.RUN.p, weightedCorr_null.RUN.p, coveredLen.RUN.p, coveredLen_null.RUN.p, jumpDist_null.RUN.p] = BDreplayDetect_vMarch21(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);

% save(fullfile(subfolder, ['Bayesian_2_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
%     'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
%     'jumpDist_null', 'spatialTunings_RL', 'spatialTunings_LR',  'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate_LR', 'runTemplate_RL', 'conslapsRatio_LR', 'diffWithAvg_LR', 'conslapsRatio_RL', 'diffWithAvg_RL', ...
%     'spatialInfo_RL', 'spatialInfo_LR', 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')


% POST
 
[posteriorProbMatrix.POST.data, RadonTF.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
    onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, RadonTF_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data, jumpDist_null.POST.data] = BDreplayDetect_vMarch21(POSTbinnedPBEs(1:3000, :), spatialTunings_RL, spatialTunings_LR, [], binDur);
    
% [posteriorProbMatrix.POST.p, RadonTF.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p, ...
%     onlyLineElements.POST.p, begPosition.POST.p, endPosition.POST.p, RadonTF_null.POST.p, weightedCorr_null.POST.p, coveredLen.POST.p, coveredLen_null.POST.p, jumpDist_null.POST.p] = BDreplayDetect_vMarch21(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);



save(fullfile(subfolder, ['Bayesian_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'secondaryPBEs', 'behavior', 'posteriorProbMatrix', 'RadonTF', 'weightedCorr', ...
    'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'RadonTF_null', 'weightedCorr_null', 'coveredLen_null', ...
    'jumpDist_null', 'spatialTunings_RL', 'spatialTunings_LR',  'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate_LR', 'runTemplate_RL', 'conslapsRatio_LR', 'diffWithAvg_LR', 'conslapsRatio_RL', 'diffWithAvg_RL', ...
    'spatialInfo_RL', 'spatialInfo_LR', 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir', 'fileinfo', 'shanks', 'clus', 'secondaryPBEs_RUN', 'secondaryPBEs_PRE', 'secondaryPBEs_POST')


end




% 
% %% Plot the figures
% 
% 
% 
% %% (1) Example PBEs sorted based on their sequence scores 
% 
% close all
% binDur = 0.02;
% 
% 
% datasets.PRE = PREbinnedPBEs_p;
% datasets.RUN = RUNbinnedPBEs_p;
% datasets.POST = POSTbinnedPBEs_p;
% 
% periods = {'PRE';'RUN'; 'POST'};
% 
% PBELen = struct('PRE', [], 'RUN', [], 'POST', []);
% 
% for pp = 1:3
%     
%     currentPBEs = datasets.(periods{pp});
%     nPBEs = size(currentPBEs, 1);
%     
%     PBELen.(periods{pp}) = zeros(nPBEs, 1);
%     nFiringUnits.(periods{pp}) = zeros(nPBEs, 1);
%     
%     for pbe = 1:nPBEs
%         
%         PBELen.(periods{pp})(pbe) = size(currentPBEs{pbe, 2}, 2);
%         
%         nFiringUnits.(periods{pp})(pbe) = numel(find(sum(currentPBEs{pbe, 2}, 2)));
%         
%     end
%     
% end
% 
% 
% secondaryPBEs_1.PRE = secondaryPBEs(secondaryPBEs(:, 1) > behavior.PREEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.PREEpoch(1, 2), :);
% secondaryPBEs_1.RUN = secondaryPBEs(secondaryPBEs(:, 1) > behavior.MazeEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.MazeEpoch(1, 2), :);
% secondaryPBEs_1.POST = secondaryPBEs(secondaryPBEs(:, 1) > behavior.POSTEpoch(1, 1) & secondaryPBEs(:, 2) < behavior.POSTEpoch(1, 2), :);
% 
% 
% acceptedEvts.PRE = find(secondaryPBEs_1.PRE(:, 8) == 1 & PBELen.PRE >= 4 & nFiringUnits.PRE >= 5);
% acceptedEvts.RUN = find(secondaryPBEs_1.RUN(:, 5) == 1 & PBELen.RUN >= 4 & nFiringUnits.RUN >= 5);
% acceptedEvts.POST = find(secondaryPBEs_1.POST(:, 8) == 1 & PBELen.POST >= 4 & nFiringUnits.POST >= 5);
% 
% 
% 
% shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
% BDscoreMethods = {'weightedCorr'; 'RadonTF'};
% 
% 
% 
% for pp = 1%: numel(periods)
%     for ii = 1: numel(shuffleMethods)
%         for jj = 1: numel(BDscoreMethods)
% 
%             period = periods{pp};
%             shuffleMethod = shuffleMethods{ii};
%             BDscoreMethod = BDscoreMethods{jj};
% 
% 
% %             directory    = fullfile(mainDir, 'Bayesian_combDir_weightedDistCorr1', period);
% %             mkdir(directory)
% 
%             
%             originalIndices = 1:numel(acceptedEvts.(period));
%             
% 
%             score2sortPBEs = BDseqscore.(period).p.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.(period));
%             score2sortPBEs = max(score2sortPBEs, [], 2);
%             
%             
%             
%             % Indices of PBEs with highest sequence scores
% %             
% %             aboveThreshPBEs = find(score2sortPBEs > 90);
% %             
% %             if numel(aboveThreshPBEs) < 60
% %                aboveThreshPBEs = find(score2sortPBEs > 70);
% %             end
% %             
% %             PBEidx2plot = aboveThreshPBEs(randperm(numel(aboveThreshPBEs), min(60, numel(aboveThreshPBEs))));
%         
%             [~, idx]    = sort(score2sortPBEs, 'descend'); 
% %             PBEidx2plot = idx(1:60);
%             PBEidx2plot = acceptedEvts.(period)(idx(1:min(60, numel(idx))));
% 
%             % % plot individual events
%                         
%             textOnTop = sprintf('Session: %s\nPeriod: %s\nPBEs were sorted based on their scores using %s and %s', fileinfo.name, period, shuffleMethod, BDscoreMethod);
% %             textOnTop = sprintf('Session: %s\nPeriod: %s\nrandom PBEs with scores above 90 percentile were plotted. The scores were calculated using %s and %s', sessionName, period, shuffleMethod, BDscoreMethod);
% 
%             filename = [period '_' shuffleMethod '_' BDscoreMethod];
%             
%             sessionNumber = 1;
%             subfolder = ['/home/kouroshmaboudi/Documents/HMM_project/GreatLakes_firstRun_Nov2020/Grosmark_originalClusters/' fileinfo.name '/BayesianDecoding'];
%             
%             
%             if ismember(sessionNumber, [2 5 8])
%                 tempPlotBD3_1D(datasets.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, [], [], runTemplate, PBEidx2plot, binDur, subfolder, filename, textOnTop) % directory
%             else
%                 tempPlotBD3(datasets.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, [], [], [], runTemplate_LR, runTemplate_RL, PBEidx2plot, binDur, subfolder, filename, textOnTop) % directory
%             end
% %             
% 
%         end
%     end
% end
% 
% close all
% 
% 
% %% (2) distribution of sequence scores
% 
% 
% shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
% BDscoreMethods = {'weightedCorr'; 'RadonTF'};
% 
% 
% for ii = 1:numel(shuffleMethods)
%     for jj = 1:numel(BDscoreMethods)
% 
%         shuffleMethod = shuffleMethods{ii};
%         BDscoreMethod = BDscoreMethods{jj};
%         
%         
%         figure;
%         set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])
% 
% 
%         % % PRE
%         subplot(3,3,[4 7])
%         plotseqscoredist(BDseqscore.PRE.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.PRE, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'PRE', length(acceptedEvts.PRE)))
%         ylabel('ratio of PBEs', 'fontsize', 10)
% 
%         % % RUN
%         subplot(3,3,[5 8])
%         plotseqscoredist(BDseqscore.RUN.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.RUN, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'RUN', length(acceptedEvts.RUN)))
% 
%         % % POST
%         subplot(3,3,[6 9])
%         plotseqscoredist(BDseqscore.POST.data.(shuffleMethod).(BDscoreMethod).prctilescore(acceptedEvts.POST, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'POST', length(acceptedEvts.POST))) 
% 
%         
%         % legend
%         legendSub = subplot(3,3,3);
% 
%         hold on
% 
% %         p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
%         p_pts = plot(1, nan, 'color', 'r', 'linestyle', ':', 'linewidth', 2);
%         p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);
% 
%         h_legend = legend([p, p_pts],{'data', 'diagonal'}, 'Location', 'South');
% 
% 
%         set(h_legend, 'fontsize', 8)
%         legend boxoff 
% 
%         set(legendSub, 'Visible', 'off');
% 
%         filename = fullfile(subfolder, [BDscoreMethod '-' shuffleMethod]);
%         
%         print(gcf, filename, '-dpdf')
%         
% 
% 
%     end
% end
% 
% close all
% 
% 
% %% scatter plots for comparing scores between different shuffling methods
% % (for a given BD replay measure) and different BD replay measures (for a given shuffling method)
% 
% % fileinfo.name = sessionName;
% plotScatterHistogram_RSvsWC(BDseqscore, acceptedEvts, fileinfo, subfolder)
% 
% plotScatterHistogram_TSvsPF(BDseqscore, acceptedEvts, fileinfo, subfolder)
% 
% close all
% 
