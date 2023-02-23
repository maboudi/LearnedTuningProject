function congruencePlot(HMMprctile, HMMprctile_p, HMMprctile_ts, HMMprctile_tc, okLegend, curr_title)

hold on

bins = 0:1:100;

h= hist(HMMprctile*100, bins);
h = h/sum(h)*100;
cumhist = 100-cumsum(h);
% threshBin = find(bins > 50, 1, 'first');
% crosspnt = cumhist(threshBin);

h_p = hist(HMMprctile_p*100, bins);
h_p = h_p/sum(h_p)*100;
cumhist_p = 100-cumsum(h_p);
% crosspnt_p = cumhist_p(threshBin);


h_ts = hist(HMMprctile_ts*100, bins);
h_ts = h_ts/sum(h_ts)*100;
cumhist_ts = 100-cumsum(h_ts);
% crosspnt_ts = cumhist_ts(threshBin);


h_tc = hist(HMMprctile_tc*100, bins);
h_tc = h_tc/sum(h_tc)*100;
cumhist_tc = 100-cumsum(h_tc);
% crosspnt_tc = cumhist_tc(threshBin);

% p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9], 'linewidth', 4);


p2 = plot(bins, cumhist_ts, 'color', [100 180 255]/255, 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt_ts crosspnt_ts], 'color', [100 180 255]/255,'LineStyle', '--', 'linewidth', 2)

p3 = plot(bins, cumhist_p, 'color', [255 180 100]/255, 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt_p crosspnt_p], 'color', [255 180 100]/255,'LineStyle', '--', 'linewidth', 2)

p4 = plot(bins, cumhist_tc, 'color', [255 100 100]/255, 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt_tc crosspnt_tc], 'color', [255 100 100]/255,'LineStyle', '--', 'linewidth', 2)


p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

% 
% noEvents = length(HMMprctile);
% line([bins(threshBin) bins(threshBin)], [0 noEvents], 'color', 'k','LineStyle', '--', 'linewidth', 2)
% 

if okLegend
    legend([p1 p2 p3 p4],{'actual','w PBE time swap','poisson','Pooled time swap'}, 'Location', 'southwest')
    legend boxoff 
end

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01])
xlabel('percentile score', 'fontsize', 12)


set(gca, 'XTick', 0:20:100, 'YTick', 0:20:100, ...
    'XTickLabels', {'0', '', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '', '100'})


if okLegend
    ylabel('% sig. PBEs', 'fontsize', 12)
end

title(curr_title, 'fontsize', 12)


ylim([-5 100])
xlim([-5 100])

% ylim([0 length(HMMprctile)])
axis square

end
