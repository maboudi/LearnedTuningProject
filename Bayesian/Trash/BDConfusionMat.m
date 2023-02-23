
% this script is intended to compare the scores from two Bayesian replay
% scoring methods (parts of this code need to modified because I made chnages in the data structure for BD sequence score)

timeswapscore   = BDseqscore.POST.data.weightedCorr.prctilescore; % scores using the first method (weighted correlation and within-PBE time swap)
pfshufflescores = BDseqscore2.POST.data.weightedCorr.prctilescore; % scores using the second method (weighted correlation and unit ID shuffle)

thresh = 80; % the threshold score

tsP_uiP = numel(find(BDseqscore.POST.data.weightedCorr.prctilescore(:, 1) >= thresh & BDseqscore2.POST.data.weightedCorr.prctilescore >= thresh)); % both methods positive
tsP_uiN = numel(find(BDseqscore.POST.data.weightedCorr.prctilescore(:, 1) >= thresh & BDseqscore2.POST.data.weightedCorr.prctilescore < thresh)); % time swap positive and unit Id negative


tsN_uiP = numel(find(BDseqscore.POST.data.weightedCorr.prctilescore(:, 1) < thresh & BDseqscore2.POST.data.weightedCorr.prctilescore >= thresh));
tsN_uiN = numel(find(BDseqscore.POST.data.weightedCorr.prctilescore(:, 1) < thresh & BDseqscore2.POST.data.weightedCorr.prctilescore < thresh));


confMat = [tsP_uiP tsP_uiN;tsN_uiP tsN_uiN];

figure; imagesc(log10(confMat)); colormap('summer'); 
xlabel('unit ID shuffle')
ylabel('within-PBE time swap')

set(gca, 'fontsize', 16,'xtick', [1 2], 'xticklabel', {['>=' num2str(thresh)]; ['<' num2str(thresh)]}, ...
    'ytick', [1 2], 'yticklabel', {['>=' num2str(thresh)]; ['<' num2str(thresh)]})

text(1,1, ['n = ' num2str(tsP_uiP)], 'horizontalalignment', 'center', 'fontsize', 16)

text(2,1, ['n = ' num2str(tsP_uiN)], 'horizontalalignment', 'center', 'fontsize', 16)

text(1,2, ['n = ' num2str(tsN_uiP)], 'horizontalalignment', 'center', 'fontsize', 16)

text(2,2, ['n = ' num2str(tsN_uiN)], 'horizontalalignment', 'center', 'fontsize', 16)


print(gcf, 'timeSwap_UIshuffle_confMat80_data', '-dsvg')


