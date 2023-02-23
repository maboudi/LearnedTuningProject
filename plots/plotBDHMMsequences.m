

close all
binDur = 0.02;

%% (1) Example PBEs sorted based on their sequence scores 

%%% POST

directory    = fullfile(mainDir, 'Bayesian', 'POST', 'data');
mkdir(directory)

% seqScores    = BDseqscore.POST.data.zscore;
% maxseqScores = max(seqScores, [], 2);


pBDpHMMind_tmat  = find(HMMSeqScores.aRUN_POST.data(:, 1) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90);
pBDpHMMind_tswap = find(HMMSeqScores.aRUN_POST.data(:, 2) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90);

pBDnHMMind_tmat  = find(HMMSeqScores.aRUN_POST.data(:, 1) < 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90); 
pBDnHMMind_tswap = find(HMMSeqScores.aRUN_POST.data(:, 2) < 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90); 

nBDpHMMind_tmat  = find(HMMSeqScores.aRUN_POST.data(:, 1) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) < 90); 
nBDpHMMind_tswap = find(HMMSeqScores.aRUN_POST.data(:, 2) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) < 90); 



pBD_ptmat_ptswapind  = find(HMMSeqScores.aRUN_POST.data(:, 1) > 90 & HMMSeqScores.aRUN_POST.data(:, 2) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90);
pBD_ptmat_ntswapind  = find(HMMSeqScores.aRUN_POST.data(:, 1) > 90 & HMMSeqScores.aRUN_POST.data(:, 2) < 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90);
pBD_ntmat_ptswapind  = find(HMMSeqScores.aRUN_POST.data(:, 1) < 90 & HMMSeqScores.aRUN_POST.data(:, 2) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) > 90);
nBD_ptmat_ptswapind  = find(HMMSeqScores.aRUN_POST.data(:, 1) > 90 & HMMSeqScores.aRUN_POST.data(:, 2) > 90 & max(BDseqscore.POST.data.prctilescore, [], 2) < 90);
nBD_ntmat_ntswapind  = find(HMMSeqScores.aRUN_POST.data(:, 1) < 90 & HMMSeqScores.aRUN_POST.data(:, 2) < 90 & max(BDseqscore.POST.data.prctilescore, [], 2) < 90);


% plot individual events


%
PBEidx2plot = pBDpHMMind_tmat(1:60);
filename = fullfile(directory , 'pBDpHMMind_tmat');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = pBDpHMMind_tswap(1:60);
filename = fullfile(directory , 'pBDpHMMind_tswap');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



%
PBEidx2plot = pBDnHMMind_tmat(1:60);
filename = fullfile(directory , 'pBDnHMMind_tmat');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = pBDnHMMind_tswap(1:60);
filename = fullfile(directory , 'pBDnHMMind_tswap');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



%
PBEidx2plot = nBDpHMMind_tmat(1:60);
filename = fullfile(directory , 'nBDpHMMind_tmat');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = nBDpHMMind_tswap(1:60);
filename = fullfile(directory , 'nBDpHMMind_tswap');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



% 
PBEidx2plot = pBD_ptmat_ptswapind(1:60);
filename = fullfile(directory , 'pBD_ptmat_ptswapind');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = pBD_ptmat_ntswapind(1:60);
filename = fullfile(directory , 'pBD_ptmat_ntswapind');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = pBD_ntmat_ptswapind(1:60);
filename = fullfile(directory , 'pBD_ntmat_ptswapind');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = nBD_ptmat_ptswapind(1:60);
filename = fullfile(directory , 'nBD_ptmat_ptswapind');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)



PBEidx2plot = nBD_ntmat_ntswapind(1:60);
filename = fullfile(directory , 'nBD_ntmat_ntswapind');

tempPlotBD(POSTbinnedPBEs.data, posteriorProbMatrix.POST.data, BDseqscore.POST.data.prctilescore, begPosition.POST.data, endPosition.POST.data, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, filename)

