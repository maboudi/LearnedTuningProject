
figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 30])

nPBEsPRE  = length(HMMSeqScores.aRUN_PRE.data);
nPBEsRUN  = length(HMMSeqScores.aRUN_RUN.data);
nPBEsPOST = length(HMMSeqScores.aRUN_POST.data);



subplot(5,3,4); HMMScores2DPlot([HMMSeqScores.aRUN_PRE.data(:, 1) max(BDseqscore.PRE.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'preplay (PRE)', nPBEsPRE)); 
ylabel('HMM tmat-shuffle scores', 'fontsize', 10)

subplot(5,3,5); HMMScores2DPlot([HMMSeqScores.aRUN_RUN.data(:, 1) max(BDseqscore.RUN.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'wake replay(RUN)', nPBEsRUN)); 
ylabel('HMM tmat-shuffle scores', 'fontsize', 10)

subplot(5,3,6); HMMScores2DPlot([HMMSeqScores.aRUN_POST.data(:, 1) max(BDseqscore.POST.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'sleep replay(POST)', nPBEsPOST)); 
ylabel('HMM tmat-shuffle scores', 'fontsize', 10)



subplot(5,3,7); HMMScores2DPlot([HMMSeqScores.aRUN_PRE.data(:, 2) max(BDseqscore.PRE.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'preplay (PRE)', nPBEsPRE)); 
ylabel('HMM tswap scores', 'fontsize', 10)

subplot(5,3,8); HMMScores2DPlot([HMMSeqScores.aRUN_RUN.data(:, 2) max(BDseqscore.RUN.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'wake replay(RUN)', nPBEsRUN)); 
ylabel('HMM tswap scores', 'fontsize', 10)

subplot(5,3,9); HMMScores2DPlot([HMMSeqScores.aRUN_POST.data(:, 2) max(BDseqscore.POST.data.prctilescore, 2)], sprintf('%s\n  (n=%d)', 'sleep replay(POST)', nPBEsPOST)); 
ylabel('HMM tswap scores', 'fontsize', 10)



% title
subplot(5,3,[1 2])
set(gca, 'visible', 'off')

text(0, mean(ylim), {fileinfo.name; 'Correlation between the scores from HMM and BD'}, 'fontsize', 14)


mkdir(fullfile(directory, 'scoreDistributions'))
filename = fullfile(directory, 'scoreDistributions', 'ScoresCorrelationBwHMMnBD' );

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename, 'epsc')


% % save the variables
% 
% filename = fullfile(directory, 'scoreDistributions', 'HMMSeqScores.mat');
% save(filename, 'HMMSeqScores')

