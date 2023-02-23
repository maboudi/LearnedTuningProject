
function PlotSeqFreqvsTime(PBEtimes, PBEscores1, PBEscores2, trainPeriod1, trainPeriod2, testPeriod, offPeriods, noEpochs, tbegin, tend)

% noEpochs = 12; % epochs with duration of ~ 15 min

epochs = linspace(tbegin, tend, noEpochs+1);
epochDur = epochs(2) - epochs(1);

epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;

[~, epochIdx] = histc(PBEtimes, epochs);


% model #1
epochMean1 = zeros(noEpochs, 1);
epochErr1 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores1(epochIdx == epoch);
    epochMean1(epoch) = mean(epochPrctls);
    epochErr1(epoch) = std(epochPrctls)/length(epochPrctls);
end





% model #2
epochMean2 = zeros(noEpochs, 1);
epochErr2 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores2(epochIdx == epoch);
    epochMean2(epoch) = mean(epochPrctls);
    epochErr2(epoch) = std(epochPrctls)/length(epochPrctls);
end





hold on

minY = floor(min([epochMean1; epochMean2])*10)/10;
maxY = ceil(max([epochMean1; epochMean2])*10)/10;

for ii = 1:length(offPeriods)
    patch([offPeriods(ii,1) offPeriods(ii,2) offPeriods(ii,2) offPeriods(ii,1)]/3600, [minY minY maxY maxY], ...
           [204 229 255]/255, 'EdgeColor', 'none');
%         set(p,'FaceAlpha', 0.5)
end

h1 = errorbar(epochCenters/3600, epochMean1, epochErr1, 'color', 'k', 'linewidth', 2);

h2 = errorbar(epochCenters/3600, epochMean2, epochErr2, 'color', 'r', 'linewidth', 2);




legend([h1, h2], ['train on ' trainPeriod1], ['train on ' trainPeriod2], 'Location','southeast')
legend boxoff 

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'box', 'off', 'Layer', 'Top')
xlim([0 3])
ylim()

xlabel(['Time in ' testPeriod '(hr)'], 'fontsize', 14)
ylabel('PBE sequence score', 'fontsize', 12)



end

