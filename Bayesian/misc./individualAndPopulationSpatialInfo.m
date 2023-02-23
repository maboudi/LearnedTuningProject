function individualAndPopulationSpatialInfo(spikeStruct, laps, spatialTunings, PF_peaks, sparsity, behavior, speed, runSpeedThresh, posBinSize, binDur, direction, fileinfo, FileBase)


% EXPLAIN THE FUNCTION


noPositionBins = size(spatialTunings, 2);


nlaps               = length(laps);

RunDecodedPositions = cell(nlaps, 1);

decErr              = cell(nlaps, 1);
medianShuffleErr    = cell(nlaps, 1);
diffdecErr          = cell(nlaps, 1); % difference in decoding error between data and shuffle

decErrStats         = zeros(nlaps, 3);
posbinIdx           = cell(nlaps, 1);

noShuffles = 500;


for ii = 1: nlaps
    
    [binnedActiveRun, temp] = binRunData_1D(spikeStruct, laps(ii, 1:2), behavior, speed, runSpeedThresh, binDur, posBinSize, fileinfo);
    posbinIdx{ii}           = temp{1}(:, 1);
    % RUN time bins with speed lower than speed threshold were excluded
    % from decoding and error calculation
    
    
    RunDecodedPositions{ii} = baysDecoder(binnedActiveRun, spatialTunings, binDur); 
    RunDecodedPositions{ii} = RunDecodedPositions{ii} ./ repmat(sum(RunDecodedPositions{ii}, 1), [noPositionBins, 1]); 

    
    [~, peakDec] = max(RunDecodedPositions{ii}, [], 1);
    
%     peakDec21 = peakDec + noPositionBins;
%     peakDec22 = peakDec - noPositionBins;
%     
    Error1 = abs(peakDec' - posbinIdx{ii}) * posBinSize;
%     
%     Error21 = abs(peakDec21' - posbinIdx{ii}) * posBinSize;
%     Error22 = abs(peakDec22' - posbinIdx{ii}) * posBinSize;
%     
%     Error2 = min(Error21, Error22);
    
%     decErr{ii} = min(Error1, Error1);
    decErr{ii} = Error1;

%     decErr{ii} = abs(peakDec' - posbinIdx{ii}) * posBinSize;
    decErrStats(ii, :) = [prctile(decErr{ii}, 25) median(decErr{ii}) prctile(decErr{ii}, 75)];
    
    
    
    [noPosBins, noTimeBins] = size(RunDecodedPositions{ii});
    
    shuffleErr = zeros(noShuffles, noTimeBins);
    
    % method(1): column cycle shuffle. With this method the shuffle error
    % is not uniform across the position
%     
%     for kk = 1:noShuffles
%         
%         
%         randShift = randi(noPosBins-1, 1, noTimeBins);
%         
%         shuffleData  = column_cycle_shuffle(RunDecodedPositions{ii}, randShift);
%         
%         [~, peakDec] = max(shuffleData, [], 1);
%         
%         peakDec21 = peakDec + noPositionBins;
%         peakDec22 = peakDec - noPositionBins;
% 
%         Error1 = abs(peakDec' - posbinIdx{ii}) * posBinSize;
% 
%         Error21 = abs(peakDec21' - posbinIdx{ii}) * posBinSize;
%         Error22 = abs(peakDec22' - posbinIdx{ii}) * posBinSize;
% 
%         Error2 = min(Error21, Error22);
% 
%         shuffleErr(kk, :) = min(Error1, Error2);
% 
% %         shuffleErr(kk, :) = abs(peakDec' - posbinIdx{ii}) * posBinSize;
%     end


    % method(2): time bin shuffle
    for kk = 1:noShuffles
        
        shuffleData = genTimeSwap(RunDecodedPositions(ii));
        [~, peakDec] = max(shuffleData{1}, [], 1);
        
        
%         peakDec21 = peakDec + noPositionBins;
%         peakDec22 = peakDec - noPositionBins;
% 
        Error1 = abs(peakDec' - posbinIdx{ii}) * posBinSize;
% 
%         Error21 = abs(peakDec21' - posbinIdx{ii}) * posBinSize;
%         Error22 = abs(peakDec22' - posbinIdx{ii}) * posBinSize;
% 
%         Error2 = min(Error21, Error22);
% 
%         shuffleErr(kk, :) = min(Error1, Error2);
        shuffleErr(kk, :) = Error1;
    end
    
    
    medianShuffleErr{ii} = median(shuffleErr, 1);
    
    diffdecErr{ii} = decErr{ii}' - medianShuffleErr{ii};
end


cnctDecErr = cell2mat(decErr);
cnctposbinIdx = cell2mat(posbinIdx);
cnctRunDecodedPositions = cell2mat(RunDecodedPositions');


largePosBinSize = 10; % in cm
bulkPosBinsIdx = ceil(cnctposbinIdx/(largePosBinSize/posBinSize)); % using larger position bins (10 cm each)


uniqPosBins = unique(bulkPosBinsIdx);
positioncenters = (1:length(uniqPosBins))*largePosBinSize + largePosBinSize/2;


% data
posBinErrorStat = zeros(length(uniqPosBins), 3);
for tt = 1: length(uniqPosBins)
    
    posBinError = cnctDecErr(bulkPosBinsIdx == uniqPosBins(tt));

    posBinErrorStat(tt, :) = [prctile(posBinError, 25) nanmedian(posBinError) prctile(posBinError, 75)];
    
end


posBinErrorStat(posBinErrorStat(:, 2) > prctile(posBinErrorStat(:, 2), 95), :) = nan;


%shuffle
cnctDecErrShuffle = cell2mat(medianShuffleErr');
posBinErrorStatShuffle = zeros(length(uniqPosBins), 3);
for tt = 1: length(uniqPosBins)
    
    posBinError = cnctDecErrShuffle(bulkPosBinsIdx == uniqPosBins(tt));

    posBinErrorStatShuffle(tt, :) = [prctile(posBinError, 25) nanmedian(posBinError) prctile(posBinError, 75)]; % low quartile, median, and high quartile
end


%difference b/w data and shuffle
cnctDecErrDiff = cell2mat(diffdecErr');
posBinErrorStatDiff = zeros(length(uniqPosBins), 3);
for tt = 1: length(uniqPosBins)
    
    posBinError = cnctDecErrDiff(bulkPosBinsIdx == uniqPosBins(tt));

    posBinErrorStatDiff(tt, :) = [prctile(posBinError, 25) nanmedian(posBinError) prctile(posBinError, 75)]; % low quartile, median, and high quartile
end

openfig(fullfile(FileBase, [fileinfo.name '_placeFields1D_' direction]));
ax1 = gca;
fig1 = get(ax1, 'children');



figure;
set(gcf, 'position', [890 549 1019 476])


% plotting the place fields

ax1 = subplot(2,7, [1 2 8 9]);
copyobj(fig1, ax1)

set(ax1, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')

xlim([0 (noPositionBins+5)*posBinSize+posBinSize/2])

t = title({['(a) Place fields ' direction]; ''}, 'fontsize', 10);
set(t, 'horizontalAlignment', 'right')

xlabel('Position on track (cm)', 'fontsize', 10)

oldPos = ax1.Position;
ax1.Position  = [oldPos(1)-oldPos(3)/4 oldPos(2) oldPos(3) oldPos(4)];

minProb = min(cnctRunDecodedPositions(:));
maxProb = max(cnctRunDecodedPositions(:));
clim = [minProb minProb+0.5*(maxProb - minProb)];


for kk = 1:5
    ax = subplot(2,7,kk+9);
    
    [noPosBins, noTimeBins] = size(RunDecodedPositions{nlaps-(5-kk)});
    
    imagesc((1:noTimeBins)*binDur-binDur/2, (1:noPosBins)*posBinSize-posBinSize/2, RunDecodedPositions{nlaps-(5-kk)}, clim) % 
    set(gca, 'YDir', 'normal')
    colormap('hot')
    
    hold on

    plot((1:noTimeBins)*binDur, posbinIdx{nlaps-(5-kk)}*posBinSize, '--', 'linewidth', 1, 'color', 'y')
    xlim([0 noTimeBins*binDur])
    
    xlabel('Time (sec)', 'fontsize', 10)
    
    if kk == 1
        ylabel({'Position on track (cm)'}, 'fontsize', 10)
        t = title({'(d) Example traversals', ''}, 'fontsize', 10);
        set(t, 'horizontalAlignment', 'center')
    end
%     
%     
%     oldPos = ax.Position;
%     ax.Position = [oldPos(1)+(oldPos(3))/3 oldPos(2) oldPos(3) oldPos(4)];
%  

end



ax2 = subplot(2,7, [3 4]);

% [counts, bins] = hist(sparsity);
% counts = counts./sum(counts);
% 
% 
% bar(bins, counts, 'FaceColor', 'k', 'EdgeColor', 'none')
% 
% xlabel('Sparsity', 'fontsize', 10)
% ylabel('Probability', 'fontsize', 10)
% 
% xlim([0 1])
% 
% t= title({'(b) Spatial information'; ''}, 'fontsize', 10);


PosBinWidth = 10;
bins = 0 : PosBinWidth : noPosBins*posBinSize;
counts = hist(PF_peaks, bins);
% counts = counts./sum(counts);

bar(bins, counts, 'FaceColor', 'k', 'EdgeColor', 'none')

xlabel('Position on track (cm)', 'fontsize', 10)
ylabel('No. of units', 'fontsize', 10)

ylim([0 1.05*max(counts)])

t= title({'(b) Dist. of preferred firing location'; ''}, 'fontsize', 10);
set(t, 'horizontalAlignment', 'center')

set(gca, 'box', 'off')

oldPos = ax2.Position;
ax2.Position = [oldPos(1) oldPos(2) oldPos(3)*0.9 oldPos(4)*0.9];


ax3 = subplot(2,7,[6 7]);


hold on

h1 = plot(positioncenters, posBinErrorStat(:, 2), 'k', 'linewidth', 2);
% plot(positioncenters, posBinErrorStat(:, 1), '--k', 'linewidth', 1)
% plot(positioncenters, posBinErrorStat(:, 3), '--k', 'linewidth', 1)

% 
% h2 = plot(positioncenters, posBinErrorStatShuffle(:, 2), 'color', [.7 .7 .7], 'linewidth', 2);
% plot(positioncenters, posBinErrorStatShuffle(:, 1), 'color', [.7 .7 .7], 'linestyle', '--', 'linewidth', 1)
% plot(positioncenters, posBinErrorStatShuffle(:, 3), 'color', [.7 .7 .7], 'linestyle', '--', 'linewidth', 1)

% 
% h3 = plot(positioncenters, posBinErrorStatDiff(:, 2), 'r', 'linewidth', 2);
% plot(positioncenters, posBinErrorStatDiff(:, 1), '--r', 'linewidth', 1)
% plot(positioncenters, posBinErrorStatDiff(:, 3), '--r', 'linewidth', 1)

hold off

xlim([0 positioncenters(end)+largePosBinSize/2])
ylim([0 max(20, max(posBinErrorStat(:, 2)))])

ylabel('Median decoding error (cm)', 'fontsize', 10)
xlabel('Position on track (cm)')


t= title({'(c) Ensemble decoding of location'; ''}, 'fontsize', 10);
set(t, 'horizontalAlignment', 'center')


set(gca, 'box', 'off')

% legend([h1, h2, h3], {'Data', 'Shuffle', 'Difference'})
% legend([h1, h2], {'Data', 'Shuffle'})
% legend('boxoff')


oldPos = ax3.Position;
ax3.Position = [oldPos(1) oldPos(2) oldPos(3)*0.9 oldPos(4)*0.9];


oldPos = ax3.Position;
ax3.Position  = [oldPos(1)-oldPos(3)/3 oldPos(2) oldPos(3) oldPos(4)];

% 
% ax4 = subplot(2,7,7);
% hold on
% 
% plot(decErrStats(:, 2), 'k', 'linewidth', 3);
% % plot(decErrStats(:, 1), '--k', 'linewidth', 1)
% % plot(decErrStats(:, 3), '--k', 'linewidth', 1)
% 
% xlabel('Laps', 'fontsize', 10)
% ylabel('Decoding error(cm)', 'fontsize', 10)
% 
% 
% oldPos = ax4.Position;
% ax4.Position = [oldPos(1) oldPos(2)-(oldPos(4)-oldPos(2))/10 oldPos(3) oldPos(4)+(oldPos(4)-oldPos(2))/10];
% 
% set(gca, 'box', 'off')



if strcmp(direction, 'LR')
    
    filename = fullfile(FileBase, 'spatialTuning_LR');

elseif strcmp(direction, 'LR')
    
    filename = fullfile(FileBase, 'spatialTuning_RL');
else
    filename = fullfile(FileBase, 'spatialTuning');
end

savepdf(gcf, filename, '-dsvg')
saveas(gcf, filename,  'epsc')
% savefig(gcf, filename)



end




