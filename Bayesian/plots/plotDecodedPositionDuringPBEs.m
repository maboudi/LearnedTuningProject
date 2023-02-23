function plotDecodedPositionDuringPBEs(peakDecodedPositionPRE, peakDecodedPositionRUN, peakDecodedPositionPOST, PosBinWidth, fileinfo, subfolder)


figure;
set(gcf, 'position', [0 0 1000 200])

ax1 = subplot(1,3,1);

bins = 0:PosBinWidth:max(peakDecodedPositionPRE);
counts = hist(peakDecodedPositionPRE, bins);
counts = counts/sum(counts);

bar(bins, counts, 'edgecolor', 'none', 'facecolor', 'k')

ylim([min(counts) 1.05*max(counts)])

t = title('PRE', 'fontsize', 10);
set(t, 'horizontalAlignment', 'left')


xlabel('Position on track (cm)', 'fontsize', 10)
ylabel('Ratio of total time bins', 'fontsize', 10)

ax2 = subplot(1,3,2);

bins = 0:PosBinWidth:max(peakDecodedPositionRUN);
counts = hist(peakDecodedPositionRUN, bins);
counts = counts/sum(counts);

bar(bins, counts, 'edgecolor', 'none', 'facecolor', 'k')

ylim([min(counts) 1.05*max(counts)])

t = title('RUN', 'fontsize', 10);
set(t, 'horizontalAlignment', 'left')

xlabel('Position on track (cm)', 'fontsize', 10)
% ylabel('Ratio of PBEs', 'fontsize', 10)


ax3 = subplot(1,3,3);

bins = 0:PosBinWidth:max(peakDecodedPositionPOST);
counts = hist(peakDecodedPositionPOST, bins);
counts = counts/sum(counts);

bar(bins, counts, 'edgecolor', 'none', 'facecolor', 'k')

ylim([min(counts) 1.05*max(counts)])

t = title('POST', 'fontsize', 10);
set(t, 'horizontalAlignment', 'left')

xlabel('Position on track (cm)', 'fontsize', 10)
% ylabel('Ratio of PBEs', 'fontsize', 10)

linkaxes([ax1,ax2,ax3],'xy')

% suptitle('(e) Dist. of decoded location during offline popluation burst events')

filename = fullfile(subfolder, 'PBEsDistOfDecodedPositions');

savepdf(gcf, filename, '-dsvg')
savepdf(gcf, filename, '-dpng')



end