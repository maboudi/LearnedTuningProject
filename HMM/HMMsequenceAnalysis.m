
directory   = fullfile(mainDir, 'HMM');
mkdir(directory)


noActiveUnits = size(PREbinnedPBEs.data{1,2}, 1);

nPBEsPRE  = length(PREbinnedPBEs.data);
nPBEsRUN  = length(RUNbinnedPBEs.data);
nPBEsPOST = length(POSTbinnedPBEs.data);

nEvtsaRUN = length(aRUNbinnedBouts.data); % it should be periods of high theta instead of the whole lap


%% Example Hidden markov Models trained separately on PRE, RUN, and POST periods


% RUN (active RUN/theta periods)

nStates = 40;

[transmat0, lambda0, prior0]  = initialModel(nStates, noActiveUnits);

[~, prior, transmatActiveRUN, lambdaActiveRUN]     = phmm_em(aRUNbinnedBouts.data(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 200);
[transmatActiveRUN, lambdaActiveRUN]               = sortStates(transmatActiveRUN, lambdaActiveRUN, prior);



% PRE
nStates   = 40; % arbitrary number of states: could be optimised 

[transmat0, lambda0, prior0]       = initialModel(nStates, noActiveUnits);

[~, prior, transmatPRE, lambdaPRE] = phmm_em(PREbinnedPBEs.data(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 200);
[transmatPRE, lambdaPRE]           = sortStates(transmatPRE, lambdaPRE, prior);


% RUN (PBEs)
nStates = 40;

[transmat0, lambda0, prior0]       = initialModel(nStates, noActiveUnits);

[~, prior, transmatRUN, lambdaRUN] = phmm_em(RUNbinnedPBEs.data(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 200);
[transmatRUN, lambdaRUN]           = sortStates(transmatRUN, lambdaRUN, prior);


% POST
nStates = 40;

[transmat0, lambda0, prior0]         = initialModel(nStates, noActiveUnits);

[~, prior, transmatPOST, lambdaPOST] = phmm_em(POSTbinnedPBEs.data(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 200);
[transmatPOST, lambdaPOST]           = sortStates(transmatPOST, lambdaPOST, prior);


%%% plot example transition and observation matrices from each period

plotExampleHMM(transmatPRE, lambdaPRE, transmatRUN, lambdaRUN, transmatPOST, lambdaPOST, transmatActiveRUN, lambdaActiveRUN, directory)



%% HMM sequence scores


HMMSeqScores = struct('aRUN_PRE', [], 'aRUN_RUN', [], 'aRUN_POST', [], 'aRUN_CV', [], ...
    'PRE_aRUN', [], 'PRE_RUN', [],'PRE_POST', [], 'PRE_CV', [], ...
    'RUN_PRE', [], 'RUN_aRUN', [], 'RUN_POST', [], 'RUN_CV', [], ...
    'POST_PRE', [], 'POST_aRUN', [], 'POST_RUN', [], 'POST_CV', []);


HMMSeqScores = structfun(@(x)(struct('data', [], 'pts', [])), HMMSeqScores,'UniformOutput', false);




%% Cross-validation (training and testing models on same periods) to characterize consistency of sequential structures within each period


nFolds    = 5; %% 5-fold cross-validation
nStates   = 40;

nShuffles      = 500;
shuffleMethods = {'tmat-shuffle', 'WithinPBE-timeSwap'};


for sm = 1:2 % shuffle method = 1: transition matrix shuffle, 2: within-PBE time swap


% PRE 
[HMMSeqScores.PRE_CV.data(:, sm), HMMSeqScores.PRE_CV.pts(:, sm), PRE_statesProbDists]    = HMMCrossValidAnalysis(PREbinnedPBEs.data,  nStates, nShuffles, nFolds, directory, 'PRE', sm);

% RUN
[HMMSeqScores.RUN_CV.data(:, sm), HMMSeqScores.RUN_CV.pts(:, sm), RUN_statesProbDists]    = HMMCrossValidAnalysis(RUNbinnedPBEs.data, nStates, nShuffles, nFolds, directory, 'RUN', sm);

% POST
[HMMSeqScores.POST_CV.data(:, sm), HMMSeqScores.POST_CV.pts(:, sm), POST_statesProbDists] = HMMCrossValidAnalysis(POSTbinnedPBEs.data, nStates, nShuffles, nFolds, directory, 'POST', sm);

% active RUN
[HMMSeqScores.aRUN_CV.data(:, sm), HMMSeqScores.aRUN_CV.pts(:, sm), aRUN_statesProbDists] = HMMCrossValidAnalysis(aRUNbinnedBouts.data, nStates, nShuffles, nFolds, directory, 'aRUN', sm);


% note that StatesProbDists are only for visualizing if the states decoded
% in consequitive time bins follow a sequential structure. We have 5 folds and so 5 different models with states not been guarantieed to be the same  



%% plot the distributions (CDFs) of HMM cross-valiation scores

figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 21 9])



% PRE
subplot(3,4,[5 9]); plotseqscoredist_HMM(HMMSeqScores.PRE_CV.data(:, sm), HMMSeqScores.PRE_CV.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'PRE', nPBEsPRE))
ylabel('ratio of PBEs', 'fontsize', 10)

% RUN
subplot(3,4,[6 10]); plotseqscoredist_HMM(HMMSeqScores.RUN_CV.data(:, sm), HMMSeqScores.RUN_CV.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'RUN', nPBEsRUN))

% POST
subplot(3,4,[7 11]); plotseqscoredist_HMM(HMMSeqScores.POST_CV.data(:, sm), HMMSeqScores.POST_CV.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'POST', nPBEsPOST))

% active RUN
subplot(3,4,[8 12]); plotseqscoredist_HMM(HMMSeqScores.aRUN_CV.data(:, sm), HMMSeqScores.aRUN_CV.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'aRUN', nEvtsaRUN))



% legend
legendSub = subplot(3,4,4);


hold on

p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

h_legend = legend([p, p_pts],{'data', 'pooled time swap'}, 'Location', 'South');
set(h_legend, 'fontsize', 8)
legend boxoff 

set(legendSub, 'Visible', 'off');


% title
subplot(3,4,[1 2 3])
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; sprintf('shuffling method = %s', shuffleMethods{sm})}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))

filename = fullfile(directory, 'scoreDistributions', ['cross-validation_' shuffleMethods{sm}] );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')


end


%% Train and test HMMs on different periods



nShuffle = 500;

for sm = 1:2 % 1: tmat shuffle, 2: within-event time swap


%% train on activeRUN


%%% train activeRUN___ test PRE
HMMSeqScores.aRUN_PRE.data(:, sm)  = modelCongruence(PREbinnedPBEs.data, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'PRE', 0, sm);
HMMSeqScores.aRUN_PRE.pts(:, sm)   = modelCongruence(PREbinnedPBEs.pts, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'PRE', 0, sm);


%%% train activeRUN___ test RUN(PBEs)
HMMSeqScores.aRUN_RUN.data(:, sm)  = modelCongruence(RUNbinnedPBEs.data, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'RUN', 0, sm);
HMMSeqScores.aRUN_RUN.pts(:, sm)   = modelCongruence(RUNbinnedPBEs.pts, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'RUN', 0, sm);


%%% train activeRUN___ test POST
HMMSeqScores.aRUN_POST.data(:, sm) = modelCongruence(POSTbinnedPBEs.data, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'POST', 0, sm);
HMMSeqScores.aRUN_POST.pts(:, sm)  = modelCongruence(POSTbinnedPBEs.pts, transmatActiveRUN, lambdaActiveRUN, nShuffle, directory, 'aRUN', 'POST', 0, sm);



%% train on PRE

%%% train PRE ___ test active RUN
HMMSeqScores.PRE_aRUN.data(:, sm)  = modelCongruence(aRUNbinnedBouts.data, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'aRUN', 0, sm);
HMMSeqScores.PRE_aRUN.pts(:, sm)   = modelCongruence(aRUNbinnedBouts.pts, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'aRUN', 0, sm);


%%% train PRE ___ test RUN
HMMSeqScores.PRE_RUN.data(:, sm)   = modelCongruence(RUNbinnedPBEs.data, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'RUN', 0, sm);
HMMSeqScores.PRE_RUN.pts(:, sm)    = modelCongruence(RUNbinnedPBEs.pts, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'RUN', 0, sm);


%%% train PRE ___ test POST
HMMSeqScores.PRE_POST.data(:, sm)  = modelCongruence(POSTbinnedPBEs.data, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'POST', 0, sm);
HMMSeqScores.PRE_POST.pts(:, sm)   = modelCongruence(POSTbinnedPBEs.pts, transmatPRE, lambdaPRE, nShuffle, directory, 'PRE', 'POST', 0, sm);




%% train on RUN

%%% train RUN___ test PRE
HMMSeqScores.RUN_PRE.data(:, sm)   = modelCongruence(PREbinnedPBEs.data, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'PRE', 0, sm);
HMMSeqScores.RUN_PRE.pts(:, sm)    = modelCongruence(PREbinnedPBEs.pts, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'PRE', 0, sm);

%%% train RUN___ test aRUN
HMMSeqScores.RUN_aRUN.data(:, sm)  = modelCongruence(aRUNbinnedBouts.data, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'aRUN', 0, sm);
HMMSeqScores.RUN_aRUN.pts(:, sm)   = modelCongruence(aRUNbinnedBouts.pts, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'aRUN', 0, sm);

%%% train RUN___ test POST
HMMSeqScores.RUN_POST.data(:, sm)  = modelCongruence(POSTbinnedPBEs.data, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'POST', 0, sm);
HMMSeqScores.RUN_POST.pts(:, sm)   = modelCongruence(POSTbinnedPBEs.pts, transmatRUN, lambdaRUN, nShuffle, directory, 'RUN', 'POST', 0, sm);



%% train on POST

%%% train POST ___ test PRE
HMMSeqScores.POST_PRE.data(:, sm)  = modelCongruence(PREbinnedPBEs.data, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'PRE', 0, sm);
HMMSeqScores.POST_PRE.pts(:, sm)   = modelCongruence(PREbinnedPBEs.pts, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'PRE', 0, sm);


%%% train POST ___ test active RUN
HMMSeqScores.POST_aRUN.data(:, sm) = modelCongruence(aRUNbinnedBouts.data, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'aRUN', 0, sm);
HMMSeqScores.POST_aRUN.pts(:, sm)  = modelCongruence(aRUNbinnedBouts.pts, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'aRUN', 0, sm);


%%% train POST ___ test RUN
HMMSeqScores.POST_RUN.data(:, sm)  = modelCongruence(RUNbinnedPBEs.data, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'RUN', 0, sm);
HMMSeqScores.POST_RUN.pts(:, sm)   = modelCongruence(RUNbinnedPBEs.pts, transmatPOST, lambdaPOST, nShuffle, directory, 'POST', 'RUN', 0, sm);






%% plot the distributions (CDFs) of HMM sequence scores


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 30])

subplot(5,3,4); plotseqscoredist_HMM(HMMSeqScores.aRUN_PRE.data(:, sm), HMMSeqScores.aRUN_PRE.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow PRE', nPBEsPRE)); ylabel('ratio of PBEs', 'fontsize', 10)
subplot(5,3,5); plotseqscoredist_HMM(HMMSeqScores.aRUN_RUN.data(:, sm), HMMSeqScores.aRUN_RUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow RUN', nPBEsRUN))
subplot(5,3,6); plotseqscoredist_HMM(HMMSeqScores.aRUN_POST.data(:, sm), HMMSeqScores.aRUN_POST.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow POST', nPBEsPOST))

subplot(5,3,7); plotseqscoredist_HMM(HMMSeqScores.PRE_aRUN.data(:, sm), HMMSeqScores.PRE_aRUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'PRE \rightarrow aRUN', nEvtsaRUN)); ylabel('ratio of PBEs', 'fontsize', 10)
subplot(5,3,8); plotseqscoredist_HMM(HMMSeqScores.PRE_RUN.data(:, sm), HMMSeqScores.PRE_RUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'PRE \rightarrow RUN', nPBEsRUN))
subplot(5,3,9); plotseqscoredist_HMM(HMMSeqScores.PRE_POST.data(:, sm), HMMSeqScores.PRE_POST.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'PRE \rightarrow POST', nPBEsPOST))

subplot(5,3,10); plotseqscoredist_HMM(HMMSeqScores.RUN_PRE.data(:, sm), HMMSeqScores.RUN_PRE.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'RUN \rightarrow PRE', nPBEsPRE)); ylabel('ratio of PBEs', 'fontsize', 10)
subplot(5,3,11); plotseqscoredist_HMM(HMMSeqScores.RUN_aRUN.data(:, sm), HMMSeqScores.RUN_aRUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'RUN \rightarrow aRUN', nEvtsaRUN))
subplot(5,3,12); plotseqscoredist_HMM(HMMSeqScores.RUN_POST.data(:, sm), HMMSeqScores.RUN_POST.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'RUN \rightarrow POST', nPBEsPOST))

subplot(5,3,13); plotseqscoredist_HMM(HMMSeqScores.POST_PRE.data(:, sm), HMMSeqScores.POST_PRE.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'POST \rightarrow PRE', nPBEsPRE)); ylabel('ratio of PBEs', 'fontsize', 10)
subplot(5,3,14); plotseqscoredist_HMM(HMMSeqScores.POST_aRUN.data(:, sm), HMMSeqScores.POST_aRUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'POST \rightarrow aRUN', nEvtsaRUN))
subplot(5,3,15); plotseqscoredist_HMM(HMMSeqScores.POST_RUN.data(:, sm), HMMSeqScores.POST_RUN.pts(:, sm), 0, sprintf('%s\n  (n=%d)', 'POST \rightarrow RUN', nPBEsRUN))

legendSub = subplot(5,3,3);


hold on

p_pts = plot(1, nan, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
p     = plot(1, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 3);

h_legend = legend([p, p_pts],{'data', 'pooled time swap'}, 'Location', 'South');
set(h_legend, 'fontsize', 8)
legend boxoff 

set(legendSub, 'Visible', 'off');


% title

subplot(5,3,[1 2])
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; sprintf('shuffling method = %s', shuffleMethods{sm})}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))
filename = fullfile(directory, 'scoreDistributions', ['diffPeriods' shuffleMethods{sm}] );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')


end

%% analyzing the relationship between the scores calculated using different shuffling methods



% cross-validation results


figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 21 11])



% PRE
subplot(2,4,5); HMMSeqScores.PRE_CV.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.PRE_CV.data, sprintf('%s\n  (n=%d)', 'PRE', nPBEsPRE));
ylabel('t-mat shuffle scores', 'fontsize', 10)

% RUN
subplot(2,4,6); HMMSeqScores.RUN_CV.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.RUN_CV.data, sprintf('%s\n  (n=%d)', 'RUN', nPBEsRUN));

% POST
subplot(2,4,7); HMMSeqScores.POST_CV.data(:, 3) = HMMScores2DPlot(HMMSeqScores.POST_CV.data, sprintf('%s\n  (n=%d)', 'POST', nPBEsPOST));

% active RUN
subplot(2,4,8); HMMSeqScores.aRUN_CV.data(:, 3) = HMMScores2DPlot(HMMSeqScores.aRUN_CV.data, sprintf('%s\n  (n=%d)', 'aRUN', nEvtsaRUN));



% title
subplot(2,4,1:4)
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; 'Correlation between the scores from different shuffling methods'}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))
filename = fullfile(directory, 'scoreDistributions', 'ScoresCorrelationBwMethods_cv' );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')

%%

figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 30])


subplot(5,3,4); HMMSeqScores.aRUN_PRE.data(:, 3)   = HMMScores2DPlot(HMMSeqScores.aRUN_PRE.data, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow PRE', nPBEsPRE));ylabel('t-mat shuffle scores', 'fontsize', 10)
subplot(5,3,5); HMMSeqScores.aRUN_RUN.data(:, 3)   = HMMScores2DPlot(HMMSeqScores.aRUN_RUN.data, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow RUN', nPBEsRUN));
subplot(5,3,6); HMMSeqScores.aRUN_POST.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.aRUN_POST.data, sprintf('%s\n  (n=%d)', 'aRUN \rightarrow POST', nPBEsPOST));

subplot(5,3,7); HMMSeqScores.PRE_aRUN.data(:, 3)   = HMMScores2DPlot(HMMSeqScores.PRE_aRUN.data, sprintf('%s\n  (n=%d)', 'PRE \rightarrow aRUN', nEvtsaRUN));ylabel('t-mat shuffle scores', 'fontsize', 10)
subplot(5,3,8); HMMSeqScores.PRE_RUN.data(:, 3)    = HMMScores2DPlot(HMMSeqScores.PRE_RUN.data, sprintf('%s\n  (n=%d)', 'PRE \rightarrow RUN', nPBEsRUN));
subplot(5,3,9); HMMSeqScores.PRE_POST.data(:, 3)   = HMMScores2DPlot(HMMSeqScores.PRE_POST.data, sprintf('%s\n  (n=%d)', 'PRE \rightarrow POST', nPBEsPOST));

subplot(5,3,10); HMMSeqScores.RUN_PRE.data(:, 3)   = HMMScores2DPlot(HMMSeqScores.RUN_PRE.data, sprintf('%s\n  (n=%d)', 'RUN \rightarrow PRE', nPBEsPRE));ylabel('t-mat shuffle scores', 'fontsize', 10)
subplot(5,3,11); HMMSeqScores.RUN_aRUN.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.RUN_aRUN.data, sprintf('%s\n  (n=%d)', 'RUN \rightarrow aRUN', nEvtsaRUN));
subplot(5,3,12); HMMSeqScores.RUN_POST.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.RUN_POST.data, sprintf('%s\n  (n=%d)', 'RUN \rightarrow POST', nPBEsPOST));

subplot(5,3,13); HMMSeqScores.POST_PRE.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.POST_PRE.data, sprintf('%s\n  (n=%d)', 'POST \rightarrow PRE', nPBEsPRE));ylabel('t-mat shuffle scores', 'fontsize', 10)
subplot(5,3,14); HMMSeqScores.POST_aRUN.data(:, 3) = HMMScores2DPlot(HMMSeqScores.POST_aRUN.data, sprintf('%s\n  (n=%d)', 'POST \rightarrow aRUN', nEvtsaRUN));
subplot(5,3,15); HMMSeqScores.POST_RUN.data(:, 3)  = HMMScores2DPlot(HMMSeqScores.POST_RUN.data, sprintf('%s\n  (n=%d)', 'POST \rightarrow RUN', nPBEsRUN));



% title
subplot(5,3,[1 2])
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; 'Correlation between the scores from different shuffling methods'}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))
filename = fullfile(directory, 'scoreDistributions', 'ScoresCorrelationBwMethods' );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')


% save the variables

filename = fullfile(directory, 'scoreDistributions', 'HMMSeqScores.mat');
save(filename, 'HMMSeqScores')


%%%%


% % time course of consistency of HMM sequence scores in each sleep/run chunk with the remaining chunks in the same period
% 
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
% seqScoreTimeCourse(HMMprctilePRE, secondaryPBEs, behavior.time(1, :), 'PRE', fileinfo)
% ylabel('sequence score', 'fontsize', 10)
% 
% % RUN
% 
% subplot(1,3,2)
% 
% load(fullfile(mainDir, 'PBEs', 'RUN', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(HMMprctileRUN, secondaryPBEs, behavior.time(2, :), 'RUN', fileinfo)
% 
% 
% % POST
% 
% subplot(1,3,3)
% 
% load(fullfile(mainDir, 'PBEs', 'POST', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(HMMprctilePOST, secondaryPBEs, behavior.time(3, :), 'POST', fileinfo)
% 
% 
% filename = fullfile(directory, 'scoreDistributions', 'seqScoreTimeCourse_crossValidation');
% 
% savepdf(gcf, filename, '-dsvg')
% saveas(gcf, filename,  'epsc')










% %% time course of consistency of HMM sequence scores in each sleep/run chunk with the remaining chunks in the same period
% 
% 
% % directory    = fullfile(mainDir, 'HMM', 'populationResults');
% % mkdir(directory)
% 
% figure;
% set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])
% 
% 
% 
% % PRE
% 
% subplot(1,3,1)
% 
% load(fullfile(mainDir, 'PBEs', 'PRE', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(HMMprctile_RUN_PRE, secondaryPBEs, behavior.time(1, :), 'RUN \rightarrow PRE', fileinfo)
% ylabel('sequence score', 'fontsize', 10)
% 
% 
% 
% % POST
% 
% subplot(1,3,2)
% 
% load(fullfile(mainDir, 'PBEs', 'POST', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(HMMprctile_RUN_POST, secondaryPBEs, behavior.time(3, :), 'RUN \rightarrow POST', fileinfo)
% 
% 
% 
% subplot(1,3,3)
% load(fullfile(mainDir, 'PBEs', 'PRE', [sessionName '-PBEs.mat']), 'secondaryPBEs')
% seqScoreTimeCourse(HMMprctile_POST_PRE, secondaryPBEs, behavior.time(1, :), 'POST \rightarrow PRE', fileinfo)
% 
% 
% 
% filename = fullfile(directory, 'scoreDistributions', 'seqScoreTimeCourse');
% 
% savepdf(gcf, filename, '-dsvg')
% saveas(gcf, filename,  'epsc')
% 
% 
% 
% %% Sort the states and rasters 
% 
% 
% [HMMprctilePOST2, statesProbDistsPOST2, unitSortIdx]    = HMMCrossValidAnalysisV2(POSTbinnedPBEs.data(1:1000, :), nStates, nShuffles, nFolds, directory, 'POST');
% 
% 
% % plot PBEs sorted based on their HMM prctile scores
% 
% 
% [~, idx]    = sort(HMMprctilePOST2, 'descend'); 
% PBEidx2plot = idx([1:50 60:10:length(idx)]);
% 
% 
% tempPlotHMM(POSTbinnedPBEs.data(1:1000, :), statesProbDistsPOST2, HMMprctilePOST2, activeUnits_sorted, PBEidx2plot, binDur, fullfile(directory, 'POST2'))
% 



%% functions

function [HMMprctile, HMMprctile_pts, statesProbDists] = HMMCrossValidAnalysis(eventsBinnedfiring, nStates, nShuffles, nFolds, directory, pbePeriod, sm)


subDirectory = fullfile(directory, pbePeriod);
mkdir(subDirectory)


% randomizing the temporal order of the PBEs
rndgen                     = randperm(size(eventsBinnedfiring, 1)); 
eventsBinnedfiring_rnd     = eventsBinnedfiring(rndgen, :); 



% train cross-validation models build on 80% of PBEs
[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiring_rnd, nStates, nFolds); 

% calculating congruence of remaining PBEs with the cross-validation models
[HMMprctile, statesProbDists]      = HMMcongruence_crossvalid(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'actual', subDirectory, sm);


% back to the original order
HMMprctile        = HMMprctile(rndgen); 
statesProbDists   = statesProbDists(rndgen);
% binResolvedScoresPRE = binResolvedScoresPRE(rndgen);



% calculating congruence of shuffled PBEs (pooled time swapped) with the
% cross-validation models (the same models trained on actual PBEs)
HMMprctile_pts    = HMMcongruence_crossvalid_pts(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'pooledTimeSwap', subDirectory, sm); % pts: pooled time swap

% back to the original order
HMMprctile_pts    = HMMprctile_pts(rndgen);


end


function [HMMprctile, statesProbDists, unitsSortIdx] = HMMCrossValidAnalysisV2(eventsBinnedfiring, nStates, nShuffles, nFolds, directory, pbePeriod)


subDirectory = fullfile(directory, [pbePeriod '2']);
mkdir(subDirectory)


% randomizing the temporal order of the PBEs
rndgen                     = randperm(size(eventsBinnedfiring, 1)); 
eventsBinnedfiring_rnd     = eventsBinnedfiring(rndgen, :); 



% train cross-validation models build on 80% of PBEs

[CVtransmats, CVlambdas, testEvts, unitsSortIdx] = trainCVmodelsV2(eventsBinnedfiring_rnd, nStates, nFolds);




% calculating congruence of remaining PBEs with the cross-validation models
[HMMprctile, statesProbDists]      = HMMcongruence_crossvalid(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'actual', subDirectory);


% back to the original order
HMMprctile        = HMMprctile(rndgen); 
statesProbDists   = statesProbDists(rndgen);
% binResolvedScoresPRE = binResolvedScoresPRE(rndgen);

unitsSortIdx = unitsSortIdx(:, rndgen);

% 
% 
% % calculating congruence of shuffled PBEs (pooled time swapped) with the
% % cross-validation models (the same models trained on actual PBEs)
% HMMprctile_pts    = HMMcongruence_crossvalid_pts(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'pooledTimeSwap', subDirectory); % pts: pooled time swap
% 
% % back to the original order
% HMMprctile_pts    = HMMprctile_pts(rndgen);
% 

end


