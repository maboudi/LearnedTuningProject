function plotErrorCdf(data, shuffledGamma, ylabelneeded)
% plotErrorCdf(data, shuffledOrder, shuffledGamma, ylabelneeded)

hold on

% actual data
[h1, bins] = hist(data, 100);
h1c = cumsum(h1)/sum(h1);

dataMed = median(data);
 
curve1 = plot(bins, h1c, 'linewidth', 4, 'color', [150, 150, 255]/255);
line([dataMed dataMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
line([0 dataMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
text(dataMed, 0.03, sprintf('%.1f cm', dataMed), 'color', [150, 150, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')

% % shuffled order
% [h2, bins] = hist(shuffledOrder, 100);
% h2c = cumsum(h2)/sum(h2);
% 
% shuffleMed = median(shuffledOrder);
% 
% curve2 = plot(bins, h2c, 'linewidth', 4, 'color', [150, 200, 255]/255);
% line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
% line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
% text(shuffleMed, 0.07, sprintf('%.1f cm', shuffleMed), 'color', [150, 200, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')
% 

% shuffled Gamma
[h3, bins] = hist(shuffledGamma, 100);
h3c = cumsum(h3)/sum(h3);

shuffleMed = median(shuffledGamma);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(shuffleMed, 0.1, sprintf('%.1f cm', shuffleMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')




xlim([0 max([data; shuffledGamma])])
set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])


xlabel('Decoding Error(cm)', 'fontsize', 16)

if ylabelneeded
    ylabel('Cumulative ratio', 'fontsize', 24)
    
%     legend('actual', 'shuffle state probability', 'shuffle position index', 'Location', 'northwest'); 
%     legend([curve1, curve2, curve3], 'actual', 'shuffled run bins', 'shuffled gamma dist.', 'Location', 'northwest'); 
        legend([curve1, curve3], 'actual', 'shuffled lsPFs', 'Location', 'northwest'); 


    legend boxoff 
end

axis square
end
