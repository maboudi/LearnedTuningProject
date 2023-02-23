
figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 18])

BDscoresPRE  = [BDseqscore.PRE.data.weightedCorr.prctilescore(:) BDseqscore.PRE.data.replayScore.prctilescore(:)];
BDscoresRUN  = [BDseqscore.RUN.data.weightedCorr.prctilescore(:) BDseqscore.RUN.data.replayScore.prctilescore(:)];
BDscoresPOST = [BDseqscore.POST.data.weightedCorr.prctilescore(:) BDseqscore.POST.data.replayScore.prctilescore(:)];



nPBEsPRE  = size(BDscoresPRE, 1)/2;
nPBEsRUN  = size(BDscoresRUN, 1)/2;
nPBEsPOST = size(BDscoresPOST, 1)/2;



ax1 = subplot(3,3,4); HMMScores2DPlot(BDscoresPRE, sprintf('%s\n  (n=%d)', 'PRE', nPBEsPRE)); 
ylabel('replay score', 'fontsize', 10)
xlabel('weighted correlation', 'fontsize', 10)
view(ax1, [30 30]);

subplot(3,3,5); HMMScores2DPlot(BDscoresRUN, sprintf('%s\n  (n=%d)', 'RUN', nPBEsRUN)); 
ylabel('replay score', 'fontsize', 10)
xlabel('weighted correlation', 'fontsize', 10)

subplot(3,3,6); HMMScores2DPlot(BDscoresPOST, sprintf('%s\n  (n=%d)', 'POST', nPBEsPOST)); 
ylabel('replay score', 'fontsize', 10)
xlabel('weighted correlation', 'fontsize', 10)



% title
subplot(5,3,[1 2])
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; 'Correlation between the scores from HMM and BD'}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))
filename = fullfile(directory, 'scoreDistributions', 'ScoresCorrelationBwHMMnBD' );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')
