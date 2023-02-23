

%% Replay evaluation using Bayesian decoding


subfolder = fullfile(mainDir, 'BayesianDecoding');
mkdir(subfolder)


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% two place field templates %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PBEs

binDur = 0.02;

% PRE
% PREsubset = randperm(length(PREbinnedPBEs.data), min(1000, length(PREbinnedPBEs.data)));

[posteriorProbMatrix.PRE.data, replayScore.PRE.data, weightedCorr.PRE.data, jumpDist.PRE.data, BDseqscore.PRE.data, sigMatrix.PRE.data, ...
    onlyLineElements.PRE.data, begPosition.PRE.data, endPosition.PRE.data, replayScore_null.PRE.data, weightedCorr_null.PRE.data, coveredLen.PRE.data, coveredLen_null.PRE.data] = BDreplayDetect_vFeb20(PREbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

[posteriorProbMatrix.PRE.p, replayScore.PRE.p, weightedCorr.PRE.p, jumpDist.PRE.p, BDseqscore.PRE.p, sigMatrix.PRE.p, ...
    onlyLineElements.PRE.p, begPosition.PRE.p, endPosition.PRE.p, replayScore_null.PRE.p, weightedCorr_null.PRE.p, coveredLen.PRE.p, coveredLen_null.PRE.p] = BDreplayDetect_vFeb20(PREbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);




% RUN
% RUNsubset = randperm(length(RUNbinnedPBEs.data), min(1000, length(RUNbinnedPBEs.data)));

[posteriorProbMatrix.RUN.data, replayScore.RUN.data, weightedCorr.RUN.data, jumpDist.RUN.data, BDseqscore.RUN.data, sigMatrix.RUN.data, ...
    onlyLineElements.RUN.data, begPosition.RUN.data, endPosition.RUN.data, replayScore_null.RUN.data, weightedCorr_null.RUN.data, coveredLen.RUN.data, coveredLen_null.RUN.data] = BDreplayDetect_vFeb20(RUNbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);

[posteriorProbMatrix.RUN.p, replayScore.RUN.p, weightedCorr.RUN.p, jumpDist.RUN.p, BDseqscore.RUN.p, sigMatrix.RUN.p, ...
    onlyLineElements.RUN.p, begPosition.RUN.p, endPosition.RUN.p, replayScore_null.RUN.p, weightedCorr_null.RUN.p, coveredLen.RUN.p, coveredLen_null.RUN.p] = BDreplayDetect_vFeb20(RUNbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);



% POST
% POSTsubset = randperm(length(POSTbinnedPBEs.data), min(1000, length(POSTbinnedPBEs.data)));
 
[posteriorProbMatrix.POST.data, replayScore.POST.data, weightedCorr.POST.data, jumpDist.POST.data, BDseqscore.POST.data, sigMatrix.POST.data, ...
    onlyLineElements.POST.data, begPosition.POST.data, endPosition.POST.data, replayScore_null.POST.data, weightedCorr_null.POST.data, coveredLen.POST.data, coveredLen_null.POST.data] = BDreplayDetect_vFeb20(POSTbinnedPBEs.data, spatialTunings_RL, spatialTunings_LR, [], binDur);
    
[posteriorProbMatrix.POST.p, replayScore.POST.p, weightedCorr.POST.p, jumpDist.POST.p, BDseqscore.POST.p, sigMatrix.POST.p, ...
    onlyLineElements.POST.p, begPosition.POST.p, endPosition.POST.p, replayScore_null.POST.p, weightedCorr_null.POST.p, coveredLen.POST.p, coveredLen_null.POST.p] = BDreplayDetect_vFeb20(POSTbinnedPBEs.p, spatialTunings_RL, spatialTunings_LR, [], binDur);


PREbinnedPBEs  = PREbinnedPBEs.data;
RUNbinnedPBEs  = RUNbinnedPBEs.data;
POSTbinnedPBEs = POSTbinnedPBEs.data;

save(fullfile(subfolder, ['Bayesian_' datestr(now, 'dd-mmm-yyyy') '.mat']), 'posteriorProbMatrix', 'replayScore', 'weightedCorr', 'jumpDist', 'BDseqscore', 'sigMatrix', 'onlyLineElements', 'coveredLen', 'begPosition', 'endPosition', 'replayScore_null', 'weightedCorr_null', 'coveredLen_null', 'spatialTunings_RL', 'spatialTunings_LR', 'PREbinnedPBEs', 'RUNbinnedPBEs', 'POSTbinnedPBEs', 'runTemplate_LR', 'runTemplate_RL')

