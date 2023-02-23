function pathHist(realData, poissonData, timeswapData, xAxis)


[counts, xbins] = hist(realData, 10);
% counts = counts./sum(counts);


counts_p = hist(poissonData, xbins);
% counts_p = counts_p./sum(counts_p);

counts_ts = hist(timeswapData, xbins);
% counts_ts = counts_ts./sum(counts_ts);

ymax = max([counts counts_p counts_ts])


figure;

rr = bar(xbins, [counts' counts_p' counts_ts']);

% 
% subplot(1,3,1)
% 
% bar(xbins, counts, 'FaceColor', 'k', 'EdgeColor','none')
% ylim([0 ymax])
set(gca, 'fontsize', 20); xlabel(xAxis, 'fontsize', 20); ylabel('Number of Paths', 'fontsize', 20); title('Real', 'fontsize', 20)

% 
% subplot(1,3,2)
% 
% bar(xbins, counts_p, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor','none')
% ylim([0 ymax])
% 
% set(gca, 'fontsize', 20); xlabel(xAxis, 'fontsize', 20);  title('Poisson', 'fontsize', 20)
% 
% 
% subplot(1,3,3)
% 
% bar(xbins, counts_ts, 'FaceColor', 'none', 'EdgeColor','k')
% ylim([0 ymax])
% set(gca, 'fontsize', 20); xlabel(xAxis, 'fontsize', 20);  title('Time Swap', 'fontsize', 20)

end

