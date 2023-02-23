
function plotHistCorrespondence(data, chance, firstPeriod, secondPeriod, currTitle)


bins = linspace(min([data; chance]), max([data; chance]), 50);
h = hist(data, bins)/length(data);


chance = chance(randi(length(chance), 5000, 1));
h_chance = hist(chance, bins)/length(chance);
pval = ranksum(data, chance, 'tail', 'right');


set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.005, 0.01])

hold on

bar(bins, h_chance, 'FaceColor', 'none', 'EdgeColor', 'k') 
bar(bins, h, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)


h_chance_sm = conv(h_chance, gausswindow(2,5), 'same');
h_sm = conv(h, gausswindow(2,5), 'same');
plot(bins, h_chance_sm, 'color', [0.7 0.7 0.7], 'linewidth', 3)
plot(bins, h_sm, 'k', 'linewidth', 3)

ylim([-0.01 0.2])

ylim1=get(gca,'ylim');
xlim1=get(gca,'xlim');

if pval < 0.001
    text(mean(xlim1),mean(ylim1), 'p < 0.001', 'fontsize', 10)
else
    text(mean(xlim1),mean(ylim1), sprintf('p = %.2f', pval), 'fontsize', 10)
end

xlabel('corres. w/in third period', 'fontsize', 10)
title({currTitle ; '' ; [firstPeriod ' and ' secondPeriod]}, 'fontsize', 10, 'fontweight', 'normal')

end