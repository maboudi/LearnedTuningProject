

%% validation of encoding model


%%% linearized positions


% Validate BD
% Decode the postions during animal's running and comapre it with the
% actual position of the animal. The less the error the more reliable are
% the templates to be used for replay detection during the PBEs

binDur      = 0.5; % duration of active RUN time bins
posBinSize  = 2; % in cm


directory = fullfile(mainDir, 'PlaceFields');


% % LR
% 
% laps_LR_idx = 2:2:length(laps);
% laps_LR_idx = find(ismember(laps(:,3), laps_LR_idx));
% 
% laps_LR    = laps(laps_LR_idx, :);
% 
% individualAndPopulationSpatialInfo2(spikeStruct, laps_LR, spatialTunings_LR, PF_peak_LR, sparsity_LR, behavior, speed, runSpeedThresh, posBinSize, binDur, 'LR', fileinfo, directory);
% 
% 
% close all
% 
% % RL
% 
% laps_RL_idx = 1:2:length(laps);
% laps_RL_idx = find(ismember(laps(:,3), laps_RL_idx));
% 
% laps_RL    = laps(laps_RL_idx, :);
% 
% individualAndPopulationSpatialInfo2(spikeStruct, laps_RL, spatialTunings_RL, PF_peak_RL, sparsity_RL, behavior, speed, runSpeedThresh, posBinSize, binDur, 'RL', fileinfo, directory);


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


%% Replay evaluation using Bayesian decoding


BDinitialize = struct('PRE', [], 'RUN', [], 'POST', []);
BDinitialize = structfun(@(x)(struct('data', [], 'pts', [], 'p', [], 'ts', [])), BDinitialize,'UniformOutput', false);

posteriorProbMatrix = BDinitialize;
weightedCorr        = BDinitialize;
jumpDist            = BDinitialize;
BDseqscore          = BDinitialize;
sigMatrix           = BDinitialize;



%%% PBEs

binDur = 0.02;

% PRE
[posteriorProbMatrix.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ~,weightedCorr_null.PRE.data] = BDreplayDetect_v2(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

% ____ surrogates
[posteriorProbMatrix.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p]             = BDreplayDetect_v2(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.PRE.ts, weightedCorr.PRE.ts, jumpDist.PRE.ts, BDseqscore.PRE.ts, sigMatrix.PRE.ts]        = BDreplayDetect_v2(PREbinnedPBEs.ts, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.PRE.pts, weightedCorr.PRE.pts, jumpDist.PRE.pts, BDseqscore.PRE.pts, sigMatrix.PRE.pts]   = BDreplayDetect_v2(PREbinnedPBEs.pts, spatialTunings_RL, spatialTunings_LR, [], binDur);



% RUN
[posteriorProbMatrix.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ~,weightedCorr_null.RUN.data] = BDreplayDetect_v2(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

% ____ surrogates
[posteriorProbMatrix.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p]             = BDreplayDetect_v2(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.RUN.ts, weightedCorr.RUN.ts, jumpDist.RUN.ts, BDseqscore.RUN.ts, sigMatrix.RUN.ts]        = BDreplayDetect_v2(RUNbinnedPBEs.ts, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.RUN.pts, weightedCorr.RUN.pts, jumpDist.RUN.pts, BDseqscore.RUN.pts, sigMatrix.RUN.pts]   = BDreplayDetect_v2(RUNbinnedPBEs.pts, spatialTunings_RL, spatialTunings_LR, [], binDur);



% POST
[posteriorProbMatrix.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ~, weightedCorr_null.POST.data] = BDreplayDetect_v2(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

% ____ surrogates
[posteriorProbMatrix.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p]              = BDreplayDetect_v2(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.POST.ts, weightedCorr.POST.ts, jumpDist.POST.ts, BDseqscore.POST.ts, sigMatrix.POST.ts]         = BDreplayDetect_v2(POSTbinnedPBEs.ts, spatialTunings_RL, spatialTunings_LR, [], binDur);
[posteriorProbMatrix.POST.pts, weightedCorr.POST.pts, jumpDist.POST.pts, BDseqscore.POST.pts, sigMatrix.POST.pts]    = BDreplayDetect_v2(POSTbinnedPBEs.pts, spatialTunings_RL, spatialTunings_LR, [], binDur);




% saving the variables

save(fullfile(mainDir, 'Bayesian', 'scoreDistributions', 'BDoutputs'),  'posteriorProbMatrix', 'weightedCorr', 'jumpDist', 'BDseqscore', 'sigMatrix', 'weightedCorr_null') 



%% Plot the figures


close all
binDur = 0.02;

% %% (1) Example PBEs sorted based on their sequence scores 
% 
% 
% %%% PRE
% 
% directory    = fullfile(mainDir, 'Bayesian', 'PRE', 'data');
% 
% seqScores    = BDseqscore.PRE.data.zscore;
% maxseqScores = max(seqScores, [], 2);
% 
% % Indices of PBEs with highest sequence scores
% 
% [~, idx]    = sort(maxseqScores, 'descend'); 
% PBEidx2plot = idx([1:50 60:10:length(idx)]);
% 
% 
% % plot individual events
% tempPlotBD(PREbinnedPBEs.data, posteriorProbMatrix.PRE.data, seqScores, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)
% 
% 
% 
% %%% RUN
% 
% directory    = fullfile(mainDir, 'Bayesian', 'RUN', 'data');
% 
% seqScores    = BDseqscore.RUN.data.zscore;
% maxseqScores = max(seqScores, [], 2);
% 
% % Indices of PBEs with highest sequence scores
% 
% [~, idx]    = sort(maxseqScores, 'descend'); 
% PBEidx2plot = idx([1:50 60:10:length(idx)]);
% 
% 
% % plot individual events
% tempPlotBD(RUNbinnedPBEs.data, posteriorProbMatrix.RUN.data, seqScores, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)
% 
% 
% 
% %%% POST
% 
% directory    = fullfile(mainDir, 'Bayesian', 'POST', 'data');
% 
% seqScores    = BDseqscore.POST.data.zscore;
% maxseqScores = max(seqScores, [], 2);
% 
% % Indices of PBEs with highest sequence scores
% 
% [~, idx]    = sort(maxseqScores, 'descend'); 
% PBEidx2plot = idx([1:50 60:10:length(idx)]);
% 
% 
% % plot individual events
% tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, seqScores, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)



%% (2) distribution of sequence scores
% comparing data with shuffles

directory    = fullfile(mainDir, 'Bayesian', 'scoreDistributions');
mkdir(directory)

figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


% PRE
subplot(3,3,[4 7])
plotseqscoredist(BDseqscore.PRE.data.prctilescore(:), BDseqscore.PRE.p.prctilescore(:), BDseqscore.PRE.ts.prctilescore(:), BDseqscore.PRE.pts.prctilescore(:), 0, sprintf('%s\n  (n=%d)', 'PRE', length(BDseqscore.PRE.data.prctilescore)))
ylabel('ratio of PBEs', 'fontsize', 10)

% RUN
subplot(3,3,[5 8])
plotseqscoredist(BDseqscore.RUN.data.prctilescore(:), BDseqscore.RUN.p.prctilescore(:), BDseqscore.RUN.ts.prctilescore(:), BDseqscore.RUN.pts.prctilescore(:), 0, sprintf('%s\n  (n=%d)', 'RUN', length(BDseqscore.RUN.data.prctilescore)))

% POST
subplot(3,3,[6 9])
plotseqscoredist(BDseqscore.POST.data.prctilescore(:), BDseqscore.POST.p.prctilescore(:), BDseqscore.POST.ts.prctilescore(:), BDseqscore.POST.pts.prctilescore(:), 0, sprintf('%s\n  (n=%d)', 'POST', length(BDseqscore.POST.data.prctilescore)))

% legend
legendSub = subplot(3,3,3);

hold on
p_ts  = plot(1, nan, 'color', [.7 .7 .7], 'linestyle', ':', 'linewidth', 3);
p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
p_p   = plot(1, nan, 'color', [.7 .7 .7], 'linestyle', '-', 'linewidth', 3);
p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

h_legend = legend([p, p_pts, p_ts, p_p],{'data', 'pooled time swap', 'time swap', 'poisson'}, 'Location', 'South');
set(h_legend, 'fontsize', 8)
legend boxoff 

set(legendSub, 'Visible', 'off');

filename = fullfile(directory, 'Bayesian__sequence_score_distribution');

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename,  'epsc')



%% (3) 2D significance matrices

% directory    = fullfile(mainDir, 'Bayesian', 'scoreDistributions');

% plot all significance matrices corresponding to different jump distance
% criteria

fn = fieldnames(sigMatrix.PRE.data);

for ii = 1: length(fn)


    figure;
    set(gcf, 'Units', 'centimeters', 'position', [0 0 30 20])


    % PRE
    subplot(3,4,1);  plotsigMatrix(sigMatrix.PRE.data.(fn{ii}), 'Data', 0)
    subplot(3,4,2);  plotsigMatrix(sigMatrix.PRE.pts.(fn{ii}), 'pooled time swap', 0)
    subplot(3,4,3);  plotsigMatrix(sigMatrix.PRE.ts.(fn{ii}), 'time swap', 0)
    subplot(3,4,4);  plotsigMatrix(sigMatrix.PRE.p.(fn{ii}), 'Poisson', 0)

    % RUN
    subplot(3,4,5);  plotsigMatrix(sigMatrix.RUN.data.(fn{ii}), '', 0)
    subplot(3,4,6);  plotsigMatrix(sigMatrix.RUN.pts.(fn{ii}), '', 0)
    subplot(3,4,7);  plotsigMatrix(sigMatrix.RUN.ts.(fn{ii}), '', 0)
    subplot(3,4,8);  plotsigMatrix(sigMatrix.RUN.p.(fn{ii}), '', 0)

    % POST
    subplot(3,4,9);  plotsigMatrix(sigMatrix.POST.data.(fn{ii}), '', 0)
    subplot(3,4,10); plotsigMatrix(sigMatrix.POST.pts.(fn{ii}), '', 0)
    subplot(3,4,11); plotsigMatrix(sigMatrix.POST.ts.(fn{ii}), '', 0)
    subplot(3,4,12); plotsigMatrix(sigMatrix.POST.p.(fn{ii}), '', 1)



    filename = fullfile(directory, ['Bayesian__' fn{ii}]);

    savepdf(gcf, filename, '-dsvg')
    saveas(gcf, filename,  'epsc')

end


%% (4) time course of BD sequence score (calculated based on weighted correlation) during sleep (PRE and POST) and RUN 

% directory    = fullfile(mainDir, 'Bayesian', 'scoreDistributions');


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


% PRE

subplot(1,3,1)

load(fullfile(mainDir, 'PBEs', 'PRE', [sessionName '-PBEs.mat']), 'secondaryPBEs')
seqScoreTimeCourse(BDseqscore.PRE.data.prctilescore, secondaryPBEs, behavior.time(1, :), 'PRE', fileinfo)
ylabel('sequence score', 'fontsize', 10)

% RUN

subplot(1,3,2)

load(fullfile(mainDir, 'PBEs', 'RUN', [sessionName '-PBEs.mat']), 'secondaryPBEs')
seqScoreTimeCourse(BDseqscore.RUN.data.prctilescore, secondaryPBEs, behavior.time(2, :), 'RUN', fileinfo)


% POST

subplot(1,3,3)

load(fullfile(mainDir, 'PBEs', 'POST', [sessionName '-PBEs.mat']), 'secondaryPBEs')
seqScoreTimeCourse(BDseqscore.POST.data.prctilescore, secondaryPBEs, behavior.time(3, :), 'POST', fileinfo)


filename = fullfile(directory, 'seqScoreTimeCourse');

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename,  'epsc')



%% Some other distributions similar to Grosmark's Science 2015

%% (5) Distributions of weighted correlation 
% comparison between the periods and also of each period with its surrogate shuffles

plotBargraph3(weightedCorr.PRE.data, weightedCorr_null.PRE.data, weightedCorr.RUN.data, weightedCorr_null.RUN.data, weightedCorr.POST.data, weightedCorr_null.POST.data) 
ylabel('absolute correlation')
set(gca, 'fontsize', 10)


%% (6) Distribution of sequence score

plotBargraph2(weightedCorr.PRE.data, weightedCorr_null.PRE.data, weightedCorr_null.POST.data, weightedCorr_null.POST.data)



