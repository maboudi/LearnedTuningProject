function plotseqscoredist_HMM(seqScore, seqScore_pts, okLegend, currTitle)



hold on

% allSeqScores = [seqScore; seqScore_pts];
% bins = linspace(min(allSeqScores), max(allSeqScores), 100);

bins = linspace(0, 100,100);

% data
h     = calpdf(seqScore, bins);
% pooled time swap
h_pts = calpdf(seqScore_pts, bins);


maxY = max([h h_pts]);


% statistical comparison of data to shuffles

significanceMarker_pts = calSigMarker(seqScore, seqScore_pts); % data with pooled time swap


p_pts = plot(bins, h_pts, 'color', 'k', 'linestyle', ':', 'linewidth', 3);
p     = plot(bins, h, 'color', 'k', 'linestyle', '-', 'linewidth', 3);


text(2, maxY-0.3, significanceMarker_pts, 'fontsize', 10, 'color', [.5 .5 .5])
                       


if okLegend
    
    h_legend = legend([p p_pts],{'data', 'pooled time swap'}, 'Location', 'East');
    
    set(h_legend, 'fontsize', 8)
    legend boxoff 

end

set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.01, 0.01])
xlabel('sequence score', 'fontsize', 10)

% set(gca, 'XTick', bins([1 25 50 75 100]), 'YTick', linspace(0, maxY, 5), ...
%     'XTickLabels', {sprintf('%.1f', bins(1)), sprintf('%.1f', bins(25)), sprintf('%.1f', bins(50)), sprintf('%.1f', bins(75)), sprintf('%.1f', bins(100))}, 'YTickLabels', {'0', '', '', '', sprintf('%.2f', maxY)})


set(gca, 'XTick', bins([1 25 50 75 100]), 'YTick', linspace(0, 1, 5), ...
    'XTickLabels', {'0', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '1'})





title(currTitle, 'fontsize', 10)

% ylim([-0.001 maxY+0.001])
% xlim([bins(1)-0.2 bins(end)])

ylim([-0.05 1])
xlim([-5 100])

axis square

end


function h = calpdf(seqScore, bins)

h = hist(seqScore, bins);
h = h/sum(h);
% h = conv(h, ones(10,1)/10, 'same');
h = cumsum(h);

end

function significanceMarker = calSigMarker(seqScore1, seqScore2)

pval = ranksum(seqScore1, seqScore2, 'tail', 'right');

if pval < 0.001
    significanceMarker = 'p < 0.001';
elseif pval < 0.01
    significanceMarker = 'p < 0.01';
elseif pval < 0.05
    significanceMarker = 'p < 0.05';
else
    significanceMarker = sprintf('%.2f', pval);
end

end
