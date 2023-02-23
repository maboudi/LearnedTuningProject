

%% validation of encoding model


%%% linearized positions

% Validate BD
% Decode the postions during animal's running and comapre it with the
% actual position of the animal. The less the error the more reliable are
% the templates to be used for replay detection during the PBEs


binDur      = 0.5; % duration of active RUN time bins
posBinSize  = 2; % in cm


directory = fullfile(mainDir, 'PlaceFields');


% LR

laps_LR_idx = 2:2:length(laps);
laps_LR_idx = find(ismember(laps(:,3), laps_LR_idx));

laps_LR    = laps(laps_LR_idx, :);

individualAndPopulationSpatialInfo2(spikeStruct, laps_LR, spatialTunings_LR, behavior, thetaPeriods, turningPeriods, speed, runSpeedThresh, posBinSize, binDur, 'LR', fileinfo, directory);


% RL

laps_RL_idx = 1:2:length(laps);
laps_RL_idx = find(ismember(laps(:,3), laps_RL_idx));

laps_RL    = laps(laps_RL_idx, :);

individualAndPopulationSpatialInfo2(spikeStruct, laps_RL, spatialTunings_RL, behavior, thetaPeriods, turningPeriods, speed, runSpeedThresh, posBinSize, binDur, 'RL', fileinfo, directory);




%% Distribution of decoded (Bayesian decoding) track locations pooled across the PBEs 

binDur     = 0.02; % in sec
posBinSize = 2; % in cm


% For each time bin consdiering the peak position from the most likely
% direction. The other way is to calculate the decoded positions separately
% for each direction: take a look at function PBEsLocationCoding 

% PRE PBEs
peakDecodedPositionPRE  = PBEsLocationCodingV2(PREbinnedPBEs.data, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);

% RUN PBEs
peakDecodedPositionRUN  = PBEsLocationCodingV2(RUNbinnedPBEs.data, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);

% POST PBEs
peakDecodedPositionPOST = PBEsLocationCodingV2(POSTbinnedPBEs.data, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);


%%% plot

PosBinWidth = 10;

plotDecodedPositionDuringPBEs(peakDecodedPositionPRE, peakDecodedPositionRUN, peakDecodedPositionPOST, PosBinWidth, fileinfo, directory)



%% Range of preferred position of units firing in a time bin 
% 
% varianceOfPlaceFieldsInPBEbins
% 
% 
% %% 
% 
% PFOverlaps_examplePBEs
% 
% 
% %% PF overlap
% 
% PFmultiComparedtoChance_v2


%% Replay evaluation using Bayesian decoding

% % spatial tunings unit ID shuffle DONE FOR EACH PBE SEPARATELY
% 
% nShuffles = 200;
% shuffleInclusionFRthresh = 1; % shuffling of place fields is done only between the units with firing rates above a certain threshold
% 
% 
% spatialTunings_RL(~spatialTunings_RL) = 1e-4;
% spatialTunings_LR(~spatialTunings_LR) = 1e-4;
% 
% % [spatialTunings_RLshuffles, spatialTunings_LRshuffles] = tuningsUnitIDshuffle(spatialTunings_RL, spatialTunings_LR, shuffleInclusionFRthresh, nShuffles);
% 
% [spatialTunings_RLshuffles, spatialTunings_LRshuffles] = tuningsUnitIDshuffle_normalize(spatialTunings_RL, spatialTunings_LR, shuffleInclusionFRthresh, nShuffles, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'
% 

BDinitialize = struct('PRE', [], 'RUN', [], 'POST', []);
BDinitialize = structfun(@(x)(struct('data', [], 'pts', [], 'p', [], 'ts', [])), BDinitialize,'UniformOutput', false);

posteriorProbMatrix = BDinitialize;
replayScore         = BDinitialize;
weightedCorr        = BDinitialize;
jumpDist            = BDinitialize;
BDseqscore          = BDinitialize;
sigMatrix           = BDinitialize;
onlyLineElements    = BDinitialize;
coveredLen          = BDinitialize;


%%% PBEs

binDur = 0.02;

% PRE
% PREsubset = randperm(length(PREbinnedPBEs.data), min(1000, length(PREbinnedPBEs.data)));

[posteriorProbMatrix.PRE.data, replayScore.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
    onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, replayScore_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data] = BDreplayDetect_vFeb20(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

[posteriorProbMatrix.PRE.pts, replayScore.PRE.pts, weightedCorr.PRE.pts, jumpDist.PRE.pts, BDseqscore.PRE.pts, sigMatrix.PRE.pts, ...
    onlyLineElements.PRE.pts, begPosition.PRE.pts, endPosition.PRE.pts, replayScore_null.PRE.pts, weightedCorr_null.PRE.pts, coveredLen.PRE.pts, coveredLen_null.PRE.pts] = BDreplayDetect_vFeb20(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);






% RUN
% RUNsubset = randperm(length(RUNbinnedPBEs.data), min(1000, length(RUNbinnedPBEs.data)));

[posteriorProbMatrix.RUN.data, replayScore.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
    onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, replayScore_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data] = BDreplayDetect_vFeb20(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

[posteriorProbMatrix.RUN.pts, replayScore.RUN.pts, weightedCorr.RUN.pts, jumpDist.RUN.pts, BDseqscore.RUN.pts, sigMatrix.RUN.pts, ...
    onlyLineElements.RUN.pts, begPosition.RUN.pts, endPosition.RUN.pts, replayScore_null.RUN.pts, weightedCorr_null.RUN.pts, coveredLen.RUN.pts, coveredLen_null.RUN.pts] = BDreplayDetect_vFeb20(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);



% POST
% POSTsubset = randperm(length(POSTbinnedPBEs.data), min(1000, length(POSTbinnedPBEs.data)));

[posteriorProbMatrix.POST.data, replayScore.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
    onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, replayScore_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data] = BDreplayDetect_vFeb20(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

[posteriorProbMatrix.POST.pts, replayScore.POST.pts, weightedCorr.POST.pts, jumpDist.POST.pts, BDseqscore.POST.pts, sigMatrix.POST.pts, ...
    onlyLineElements.POST.pts, begPosition.POST.pts, endPosition.POST.pts, replayScore_null.POST.pts, weightedCorr_null.POST.pts, coveredLen.POST.pts, coveredLen_null.POST.pts] = BDreplayDetect_vFeb20(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);



% saving the variables
mkdir(fullfile(mainDir, 'Bayesian_combDir2', 'scoreDistributions2'))
save(fullfile(mainDir, 'Bayesian_combDir2', 'scoreDistributions2', 'BDoutputs'),  'posteriorProbMatrix', 'replayScore', 'weightedCorr', 'jumpDist', 'BDseqscore', 'begPosition', 'endPosition', 'sigMatrix', 'replayScore_null', 'weightedCorr_null', 'onlyLineElements', 'coveredLen', 'coveredLen_null') 



%% Plot the figures


close all
binDur = 0.02;

%% (1) Example PBEs sorted based on their sequence scores 


periods = {'PRE'}; % ;'RUN';'POST'
% periods = {'RUN'};
shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
BDscoreMethods = {'replayScore'; 'weightedCorr'};


binnedPBEs2.PRE = PREbinnedPBEs.data;
binnedPBEs2.RUN = RUNbinnedPBEs.data;
binnedPBEs2.POST = POSTbinnedPBEs.data;


for pp = 1%: numel(periods)
    for ii = 1%: numel(shuffleMethods)
        for jj = 1%: numel(BDscoreMethods)

            period = periods{pp};
            shuffleMethod = shuffleMethods{ii};
            BDscoreMethod = BDscoreMethods{jj};


            directory    = fullfile(mainDir, 'Bayesian_combDir_withSeqMatchIdx2', period);
            mkdir(directory)


            score2sortPBEs = BDseqscore.(period).data.(shuffleMethod).(BDscoreMethod).prctilescore;
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
            PBEidx2plot = idx(1:60);


            % % plot individual events

            textOnTop = sprintf('Session: %s\nPeriod: %s\nPBEs were sorted based on their scores using %s and %s', sessionName, period, shuffleMethod, BDscoreMethod);
%             textOnTop = sprintf('Session: %s\nPeriod: %s\nrandom PBEs with scores above 90 percentile were plotted. The scores were calculated using %s and %s', sessionName, period, shuffleMethod, BDscoreMethod);

            filename = [period '_' shuffleMethod '_' BDscoreMethod];
            
%             tempPlotBD3(binnedPBEs2.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, begPosition.(period).data, endPosition.(period).data, seqMatchingIdx.(period).data, runTemplate_LR, runTemplate_RL, PBEidx2plot, binDur, directory, filename, textOnTop) % directory
            tempPlotBD3(binnedPBEs2.(period), posteriorProbMatrix.(period).data, BDseqscore.(period).data, weightedCorr.(period).data, begPosition.(period).data, endPosition.(period).data, [], runTemplate_LR, runTemplate_RL, PBEidx2plot, binDur, directory, filename, textOnTop) % directory

        end
    end
end



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


directory = fullfile(mainDir, 'Bayesian_combDir2', 'scoreDistributions');

if ~exist(directory, 'dir'); mkdir(directory); end

shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};
BDscoreMethods = {'replayScore'; 'weightedCorr'};

for ii = 1:numel(shuffleMethods)
    for jj = 1:numel(BDscoreMethods)

        shuffleMethod = shuffleMethods{ii};
        BDscoreMethod = BDscoreMethods{jj};
        
        
        figure;
        set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


        % % PRE
        subplot(3,3,[4 7])
        plotseqscoredist(BDseqscore.PRE.data.(shuffleMethod).(BDscoreMethod).prctilescore(:, 1), [], [], [], 0, sprintf('%s\n  (n=%d)', 'PRE', length(PREbinnedPBEs.data)))
        ylabel('ratio of PBEs', 'fontsize', 10)

        % % RUN
        subplot(3,3,[5 8])
        plotseqscoredist(BDseqscore.RUN.data.(shuffleMethod).(BDscoreMethod).prctilescore(:, 1), [], [], BDseqscore.RUN.pts.(shuffleMethod).(BDscoreMethod).prctilescore(:, 1), 0, sprintf('%s\n  (n=%d)', 'RUN', length(RUNbinnedPBEs.data)))

        % % POST
        subplot(3,3,[6 9])
        plotseqscoredist(BDseqscore.POST.data.(shuffleMethod).(BDscoreMethod).prctilescore(:, 1), [], [], BDseqscore.POST.pts.(shuffleMethod).(BDscoreMethod).prctilescore(:, 1), 0, sprintf('%s\n  (n=%d)', 'POST', length(POSTbinnedPBEs.data))) 

        
        % legend
        legendSub = subplot(3,3,3);

        hold on

        p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
        p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

        h_legend = legend([p, p_pts],{'data', 'poisson'}, 'Location', 'South');


        set(h_legend, 'fontsize', 8)
        legend boxoff 

        set(legendSub, 'Visible', 'off');

        filename = fullfile(directory, [BDscoreMethod '-' shuffleMethod]);
        
        print(gcf, filename, '-dpdf')
        
%         savepdf(gcf, filename, '-dsvg')
%         saveas(gcf, filename,  'epsc')

    end
end




%% scatter plots for comparing scores between different shuffling methods
% (for a given BD replay measure) and different BD replay measures (for a given shuffling method)


plotScatterHistogram_RSvsWC(BDseqscore, fileinfo, directory)

plotScatterHistogram_TSvsPF(BDseqscore, fileinfo, directory)



%% (3) 2D significance matrices



% plot all significance matrices corresponding to different jump distance
% criteria

shuffleMethods = {'wPBEtimeswap'; 'unitIDshuffle'};

sigMatrices = fieldnames(sigMatrix.RUN.data.(shuffleMethods{1})); % the different significance matrix 
jumpDistMatNames = {'max jump distance';'normalized max jump distance';'median jump distance'};

for ii = 1:numel(shuffleMethods)
    
    shuffleMethod = shuffleMethods{ii};
    
    for jj = 1:3
        
        currSigMatrix = sigMatrices{jj};
        
       
        figure;
        set(gcf, 'Units', 'centimeters', 'position', [0 0 20 20])


        % PRE
        subplot(3,2,1);  plotsigMatrix(sigMatrix.PRE.data.(shuffleMethod).(currSigMatrix), {fileinfo.name;shuffleMethod;'';'Data'}, 0); xlabel(jumpDistMatNames{jj}); ylabel({'PRE', '','weighted correlation'});
        subplot(3,2,2);  plotsigMatrix(sigMatrix.PRE.p.(shuffleMethod).(currSigMatrix), 'Poisson', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation'); 


        % RUN
        subplot(3,2,3);  plotsigMatrix(sigMatrix.RUN.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel({'RUN', '','weighted correlation'});
        subplot(3,2,4);  plotsigMatrix(sigMatrix.RUN.p.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation');


        % POST
        subplot(3,2,5);  plotsigMatrix(sigMatrix.POST.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel({'POST', '','weighted correlation'});
        subplot(3,2,6); plotsigMatrix(sigMatrix.POST.p.(shuffleMethod).(currSigMatrix), '', 0); xlabel(jumpDistMatNames{jj}); ylabel('weighted correlation'); 


    
%         filename = fullfile(directory, ['Bayesian__' currSigMatrix '__' shuffleMethod]);
            filename = ['Bayesian__' currSigMatrix '__' shuffleMethod];

        savepdf(gcf, filename, '-dpng')
%         saveas(gcf, filename,  'epsc')
    end
    
    jj = 4;
    
    currSigMatrix = sigMatrices{jj};
    figure;
    set(gcf, 'Units', 'centimeters', 'position', [0 0 20 20])


    % PRE
    subplot(3,2,1);  plotsigMatrix_rs(sigMatrix.PRE.data.(shuffleMethod).(currSigMatrix), {fileinfo.name; shuffleMethod;'';'Data'}, 0); xlabel('covered track length'); ylabel({'PRE', '','replay score'});
    subplot(3,2,2);  plotsigMatrix_rs(sigMatrix.PRE.p.(shuffleMethod).(currSigMatrix), 'poisson', 0); xlabel('covered track length'); ylabel('replay score'); 


    % RUN
    subplot(3,2,3);  plotsigMatrix_rs(sigMatrix.RUN.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel({'RUN', '','replay score'});
    subplot(3,2,4);  plotsigMatrix_rs(sigMatrix.RUN.p.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel('replay score');


    % POST
    subplot(3,2,5);  plotsigMatrix_rs(sigMatrix.POST.data.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel({'POST', '','replay score'});
    subplot(3,2,6); plotsigMatrix_rs(sigMatrix.POST.p.(shuffleMethod).(currSigMatrix), '', 0); xlabel('covered track length'); ylabel('replay score'); 



%     filename = fullfile(directory, ['Bayesian__' currSigMatrix '__' shuffleMethod]);
    
    filename = ['Bayesian__' currSigMatrix '__' shuffleMethod];
    savepdf(gcf, filename, '-dpng')
%     saveas(gcf, filename,  'epsc')


end



% 
% %% (4) time course of BD sequence score (calculated based on weighted correlation) during sleep (PRE and POST) and RUN 
% 
% % directory    = fullfile(mainDir, 'Bayesian', 'scoreDistributions');
% 
% 
% figure;
% set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])
% 
% 
% % PRE
% 
% subplot(1,3,1)
% 
% load(fullfile(mainDir, 'PBEs', 'PRE', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(BDseqscore.PRE.data.prctilescore, secondaryPBEs, behavior.time(1, :), 'PRE', fileinfo)
% ylabel('sequence score', 'fontsize', 10)
% 
% % RUN
% 
% subplot(1,3,2)
% 
% load(fullfile(mainDir, 'PBEs', 'RUN', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(BDseqscore.RUN.data.prctilescore, secondaryPBEs, behavior.time(2, :), 'RUN', fileinfo)
% 
% 
% % POST
% 
% subplot(1,3,3)
% 
% load(fullfile(mainDir, 'PBEs', 'POST', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(BDseqscore.POST.data.prctilescore, secondaryPBEs, behavior.time(3, :), 'POST', fileinfo)
% 
% 
% filename = fullfile(directory, 'seqScoreTimeCourse');
% 
% savepdf(gcf, filename, '-dsvg')
% saveas(gcf, filename,  'epsc')
% 
% 
% 
% %% Some other distributions similar to Grosmark's Science 2015
% 
% %% (5) Distributions of weighted correlation 
% % comparison between the periods and also of each period with its surrogate shuffles
% 
% plotBargraph3(weightedCorr.PRE.data, weightedCorr_null.PRE.data, weightedCorr.RUN.data, weightedCorr_null.RUN.data, weightedCorr.POST.data, weightedCorr_null.POST.data) 
% ylabel('absolute correlation')
% set(gca, 'fontsize', 10)
% 
% 
% %% (6) Distribution of sequence score
% 
% plotBargraph2(weightedCorr.PRE.data, weightedCorr_null.PRE.data, weightedCorr_null.POST.data, weightedCorr_null.POST.data)



