function congruencePlot2(HMMprctile, HMMprctile_pts,  okLegend, curr_title)

hold on

bins = 0:1:100;
% bins = linspace(min([HMMprctile; HMMprctile_pts]), max([HMMprctile; HMMprctile_pts]), 100);

h = hist(HMMprctile*100, bins);
h = h/sum(h)*100;

cumhist = 100-cumsum(h); %%  for the paper
threshBin = find(bins > 50, 1, 'first');
% crosspnt = cumhist(threshBin);


h_pts = hist(HMMprctile_pts*100, bins);
h_pts = h_pts/sum(h_pts)*100;

cumhist_pts = 100-cumsum(h_pts); %% for the paper
crosspnt_pts = cumhist_pts(threshBin);


% statistical comparison of the two distributions

pval = ranksum(HMMprctile, HMMprctile_pts);

if pval < 0.001
    significanceMarker = '***';
elseif pval < 0.01                        
    significanceMarker = '**';
elseif pval < 0.05
    significanceMarker = '*';
else
    significanceMarker = 'ns';
end

% p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9],
% 'linewidth', 4); % the chance line-what would be the chance distribution
% is not established yet


p2 = plot(bins, cumhist_pts, 'color', [140 166 181]/255, 'linewidth', 3);
% line([0 bins(threshBin)], [crosspnt_pts crosspnt_pts], 'color', [140 166 181]/255,'LineStyle', '--', 'linewidth', 2)

p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 3);
% line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

text(bins(threshBin)+10, crosspnt_pts+10, significanceMarker, 'fontsize', 12)


% noEvents = length(HMMprctile);
% line([bins(threshBin) bins(threshBin)], [0 1], 'color', 'k','LineStyle', '--', 'linewidth', 2)



if okLegend
    legend([p1 p2],{'actual', 'pooled time swap'}, 'Location', 'southwest')
    legend boxoff 
end

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01])
xlabel('percentile threshold', 'fontsize', 12)

set(gca, 'XTick', 0:20:100, 'YTick', 0:20:100, ...
    'XTickLabels', {'0', '', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '', '100'})

if okLegend
    ylabel('% sig. PBEs', 'fontsize', 12)
end

title(curr_title, 'fontsize', 12)

ylim([-5 100])
xlim([-5 100])

axis square

end
