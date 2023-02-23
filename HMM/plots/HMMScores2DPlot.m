
function corrBWScores = HMMScores2DPlot(HMMprctile, ttl)

edges = {0:10:100 0:10:100};
edges{1}(end) = 101;
edges{2}(end) = 101;

counts       = hist3([HMMprctile(:, 1) HMMprctile(:, 2)], 'Edges', edges);
counts(:,11) = []; 
counts(11,:) = [];

ratios = counts./sum(counts(:));

corrBWScores = corr(HMMprctile(:, 1), HMMprctile(:, 2));

imagesc(ratios); 
text(2,2, {sprintf('rho = %.2f', corrBWScores); sprintf('max ratio = %.2f', max(ratios(:))) }, 'fontsize', 8, 'color', 'w')

set(gca, 'ydir', 'normal', 'XTick', [0.5 10.5], 'YTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTickLabel', {'0', '100'}, 'fontsize', 10)
% xlabel('Bayesian replay score', 'fontsize', 10)
% ylabel('tmat-shuffle scores', 'fontsize', 10)
title(ttl)

colormap('jet')

end