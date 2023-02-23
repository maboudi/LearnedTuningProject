function plotseqscoredist(seqScore_ts, seqScore_ui, seqScore_pf, seqScore_posterior, okLegend, currTitle)

hold on

allSeqScores = [seqScore_ts; seqScore_ui; seqScore_pf; seqScore_posterior];
allSeqScores = allSeqScores(:);
bins = linspace(min(allSeqScores), max(allSeqScores), 101);
bincenters = bins(1:end-1) + diff(bins(1:2));

% within-PBE time swap
h_ts    = calpdf(seqScore_ts, bins);

% unit ID shuffle
h_ui   = calpdf(seqScore_ui, bins);

% place field shuffle
h_pf  = calpdf(seqScore_pf, bins);

% pooled time swap
h_posterior = calpdf(seqScore_posterior, bins);


% maxY = max([h h_p h_ts h_pts]);
maxY = 1;

% statistical comparison of data to shuffles
% 
% significanceMarker_p   = calSigMarker(seqScore, seqScore_p); % data with possion
% significanceMarker_ts  = calSigMarker(seqScore, seqScore_ts); % data with time swap (within PBE)
% significanceMarker_pts = calSigMarker(seqScore, seqScore_pts); % data with pooled time swap


% line([0 100], [0 1], 'color', [.7 .7 .7], 'linestyle', ':', 'linewidth', 2)
patch([0 100 100 0], [0 0 1 0],  [.8 .8 .8], 'EdgeColor', 'none', 'FaceAlpha', 0.3)

hold on
p__ts       = plot(bincenters, h_ts, 'color', [0.4940    0.1840    0.5560 0.7], 'linestyle', '-', 'linewidth', 2, 'displayname', 'wPBE time swap');
p_ui        = plot(bincenters, h_ui, 'color', [0.9290    0.6940    0.1250 0.7], 'linestyle', '-', 'linewidth', 2, 'displayname', 'unit identity shuffle');
p_pf        = plot(bincenters, h_pf, 'color', [0.8500    0.3250    0.0980 0.7], 'linestyle', '-', 'linewidth', 2, 'displayname', 'place field c rotation');
p_posterior = plot(bincenters, h_posterior, 'color', [0    0.4470    0.7410 0.7], 'linestyle', '-', 'linewidth', 2, 'displayname', 'posterior c rotation');



%%% added Sep 2020 (KM)
% try
% ThreshPrctile = 95;
% cdfATThresh = h(find(bincenters >= ThreshPrctile, 1, 'first'));
% 
% line([ThreshPrctile ThreshPrctile], [0 1], 'color', [.5 .5 .5], 'linestyle', '-', 'linewidth', 1)
% text(ThreshPrctile+2, 0, '95%', 'color', [.5 .5 .5])
% 
% line([0 100], [cdfATThresh cdfATThresh], 'color', [.5 .5 .5], 'linestyle', '-', 'linewidth', 1)
% text(0, cdfATThresh+0.03, sprintf('%.2f', cdfATThresh), 'color', [.5 .5 .5])
% 
% catch 
%     
%     fprintf('failed')
% end
%%%


% 
% text(2, maxY-0.1, {sprintf('vs pooled time swap  %s', significanceMarker_pts)}, 'fontsize', 6) 
%                    sprintf('vs time swap         %s', significanceMarker_ts), ...
%                    sprintf('vs poisson           %s',
%                    significanceMarker_p)}, ...
                     
% 
% 
% if okLegend
% %     figureSize = get(gcf, 'position');
%     
%     h_legend = legend([p p_pts p_ts p_p],{'data', 'pooled time swap', 'time swap', 'poisson'}, 'Location', 'East');
%     
%     set(h_legend, 'fontsize', 8)
%     legend boxoff 
% 
% end

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.01, 0.01])
xlabel('replay score %', 'fontsize', 14)

% set(gca, 'XTick', bins([1 25 50 75 100]), 'YTick', linspace(0, maxY, 5), ...
%     'XTickLabels', {sprintf('%.1f', bins(1)), sprintf('%.1f', bins(25)), sprintf('%.1f', bins(50)), sprintf('%.1f', bins(75)), sprintf('%.1f', bins(100))}, 'YTickLabels', {'0', '', '', '', sprintf('%.2f', maxY)})


set(gca, 'XTick', bins([1 25 50 75 100]), 'YTick', linspace(0, maxY, 5), ...
    'XTickLabels', {'0', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '1'})





title(currTitle, 'fontsize', 14)

% ylim([-0.001 maxY+0.001])
% xlim([bins(1)-0.2 bins(end)])

ylim([-0.05 1])
xlim([-5 100])

axis square

end


function h = calpdf(seqScore, bins)

h = histc(seqScore, bins);
h(end-1) = h(end-1)+ h(end);
h(end) = [];

h = h/sum(h);
% h = conv(h, ones(10,1)/10, 'same');
h = cumsum(h);

end

function significanceMarker = calSigMarker(seqScore1, seqScore2)

pval = ranksum(seqScore1, seqScore2, 'tail', 'right');

if pval < 0.001
    significanceMarker = '< 0.001';
elseif pval < 0.01
    significanceMarker = '< 0.01';
elseif pval < 0.05
    significanceMarker = '< 0.05';
else
    significanceMarker = sprintf('%.2f', pval);
end

end
