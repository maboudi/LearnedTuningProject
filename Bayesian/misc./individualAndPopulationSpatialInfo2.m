function individualAndPopulationSpatialInfo2(spikes, laps, spatialTunings, behavior, thetaPeriods, turningPeriods, speed, runSpeedThresh, posBinSize, binDur, direction, fileinfo, FileBase)


% EXPLAIN THE FUNCTION

% spatialTunings = spatialTunings.smoothed; 

noPositionBins = size(spatialTunings, 2);

nlaps               = length(laps);

RunDecodedPositions = cell(nlaps, 1);

peakDec             = cell(nlaps, 1);
decErr              = cell(nlaps, 1);
meanDecErr          = zeros(nlaps, 1);

medianShuffleErr    = cell(nlaps, 1);
diffdecErr          = cell(nlaps, 1); % difference in decoding error between data and shuffle

decErrStats         = zeros(nlaps, 3);
posbinIdx           = cell(nlaps, 1);

noShuffles = 500;


%%%% interpolate to replace missing position data

linearPos = fileinfo.linearPos(:, 1);
tpos = fileinfo.linearPos(:, 2);

% if the pos includes nans interpolate to replace them

nanIdx = find(isnan(linearPos));

linearPos_woNans = linearPos;
linearPos_woNans(isnan(linearPos)) = [];
tpos_woNans = tpos;
tpos_woNans(isnan(linearPos)) = [];

temp = interp1(tpos_woNans, linearPos_woNans, tpos(nanIdx));

linearPos2 = linearPos;
linearPos2(nanIdx) = temp;

fileinfo.linearPos(:, 1) = linearPos2;

%%



laps = laps(randperm(nlaps), :); % randomizing order of the laps for the cross-validation procedure
laps(:, 3) = 1: nlaps;


nFolds = nlaps;

foldSize = floor(nlaps/nFolds);

for fold = 1:nFolds
  
   if fold == nFolds
       decodeLaps = laps((fold-1)*foldSize+1 : nlaps, :); % laps to test (decode based on) the Bayesian encoding models
   else
       decodeLaps = laps((fold-1)*foldSize+1 : fold*foldSize, :);
   end
   encodeLaps = laps(~ismember(1:nlaps, decodeLaps), :); 
   
  
   % calculate the spatial tunings
   spikes = spatialTuning_1D_V2(spikes, encodeLaps, thetaPeriods, turningPeriods, runSpeedThresh, direction, posBinSize, fileinfo); % spatial tuning or Bayesian encoding model
   
   foldSpatialTunings = zeros(numel(spikes), noPositionBins);
   
   for iUnit = 1:numel(spikes)
        foldSpatialTunings(iUnit, :) = spikes(iUnit).CVspatialTuning.(direction);
   end

   foldSpatialTunings(~foldSpatialTunings) = 1e-3;
   
   if fold == 1
      noPositionBins  = size(foldSpatialTunings, 2);
   end
   
   
   for ii = 1: size(decodeLaps, 1)
       
       currLap = decodeLaps(ii, 3);
%         currLap = ii;

       [binnedActiveRun, temp] = binRunData_1D(spikes, decodeLaps(ii, 1:2),runSpeedThresh, turningPeriods, binDur, 0, posBinSize, fileinfo); % need to add theta and turning periods to this ...
       
       posbinIdx{currLap}      = temp(:, 1);
       % RUN time bins with speed lower than speed threshold were excluded
       % from decoding and error calculation


%         RunDecodedPositions{currLap} = baysDecoder(binnedActiveRun{1}, foldSpatialTunings, binDur); 
        RunDecodedPositions{currLap} = baysDecoder(binnedActiveRun, foldSpatialTunings, binDur); 

        RunDecodedPositions{currLap} = RunDecodedPositions{currLap} ./ repmat(sum(RunDecodedPositions{currLap}, 1), [size(RunDecodedPositions{currLap}, 1), 1]); 


        [~, peakDec{currLap}] = max(RunDecodedPositions{currLap}, [], 1);

    %     peakDec21 = peakDec + noPositionBins;
    %     peakDec22 = peakDec - noPositionBins;
    %     
        Error1 = abs(peakDec{currLap}' - posbinIdx{currLap}) * posBinSize - posBinSize/2;
    %     
    %     Error21 = abs(peakDec21' - posbinIdx{ii}) * posBinSize;
    %     Error22 = abs(peakDec22' - posbinIdx{ii}) * posBinSize;
    %     
    %     Error2 = min(Error21, Error22);

    %     decErr{ii} = min(Error1, Error1);
        decErr{currLap} = Error1;
        meanDecErr(currLap) = mean(Error1);
    %     decErr{ii} = abs(peakDec' - posbinIdx{ii}) * posBinSize;
        decErrStats(currLap, :) = [prctile(decErr{currLap}, 25) median(decErr{currLap}) prctile(decErr{currLap}, 75)];



        [noPosBins, noTimeBins] = size(RunDecodedPositions{currLap});

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

            shuffleData = genTimeSwap(RunDecodedPositions(currLap));
            [~, peakDec_shuffle] = max(shuffleData{1}, [], 1);


    %         peakDec21 = peakDec + noPositionBins;
    %         peakDec22 = peakDec - noPositionBins;
    % 
            Error1 = abs(peakDec_shuffle' - posbinIdx{currLap}) * posBinSize;
    % 
    %         Error21 = abs(peakDec21' - posbinIdx{ii}) * posBinSize;
    %         Error22 = abs(peakDec22' - posbinIdx{ii}) * posBinSize;
    % 
    %         Error2 = min(Error21, Error22);
    % 
    %         shuffleErr(kk, :) = min(Error1, Error2);
            shuffleErr(kk, :) = Error1;
        end


        medianShuffleErr{currLap} = median(shuffleErr, 1);

        diffdecErr{currLap} = decErr{currLap}' - medianShuffleErr{currLap};
   end
   
  
end


cnctDecErr = cell2mat(decErr);

cnctposbinIdx    = cell2mat(posbinIdx);
cnctDecposbinIdx = cell2mat(peakDec');

cnctRunDecodedPositions = cell2mat(RunDecodedPositions');



% drawing the decoding error across track position

largePosBinSize = 10; % in cm

bulkPosBinsIdx    = ceil(cnctposbinIdx/(largePosBinSize/posBinSize)); % using larger position bins (10 cm each)
bulkDecPosBinsIdx = ceil(cnctDecposbinIdx/(largePosBinSize/posBinSize)); % using larger position bins (10 cm each)


uniqPosBins = 1:max(unique([bulkPosBinsIdx;bulkDecPosBinsIdx']));
positioncenters = (1:length(uniqPosBins))*largePosBinSize + largePosBinSize/2;


% data
posBinErrorStat = zeros(length(uniqPosBins), 3);
for tt = 1: length(uniqPosBins)
    
    posBinError = cnctDecErr(bulkPosBinsIdx == uniqPosBins(tt));

    posBinErrorStat(tt, :) = [prctile(posBinError, 25) nanmedian(posBinError) prctile(posBinError, 75)];
    
end


% posBinErrorStat(posBinErrorStat(:, 2) > 3*(prctile(posBinErrorStat(:, 2), 75)-prctile(posBinErrorStat(:, 2), 25)), :) = nan; 


%shuffle
cnctDecErrShuffle = cell2mat(medianShuffleErr');
posBinErrorStatShuffle = zeros(length(uniqPosBins), 3);
for tt = 1: length(uniqPosBins)
    
    posBinError = cnctDecErrShuffle(bulkPosBinsIdx == uniqPosBins(tt));

    posBinErrorStatShuffle(tt, :) = [prctile(posBinError, 25) nanmedian(posBinError) prctile(posBinError, 75)]; % low quartile, median, and high quartile
end



% decoding confusion matrix

decConfMat = full(sparse(bulkDecPosBinsIdx, bulkPosBinsIdx, ones(length(bulkPosBinsIdx), 1), length(uniqPosBins), length(uniqPosBins)));
decConfMat = decConfMat./repmat(sum(decConfMat, 2), [1, length(uniqPosBins)]); 

% decConfMat(isnan(decConfMat)) = 0;
% 
% win = gausswindow(2, 5);
% 
% for jj = 1:size(decConfMat, 1)
%     
%     decConfMat(jj, :) = conv(decConfMat(jj, :), win, 'same');
% end
    
    



%difference btw data and shuffle
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


[~, bestLaps] = sort(meanDecErr, 'ascend');
% bestLaps = bestLaps(1:(floor(length(bestLaps))/5):length(bestLaps));
bestLaps = bestLaps(1:5);


for kk = 1:5
    
    ax = subplot(2,7,kk+9);
    
    [noPosBins, noTimeBins] = size(RunDecodedPositions{bestLaps(kk)});
    
    imagesc((1:noTimeBins)*binDur-binDur/2, (1:noPosBins)*posBinSize-posBinSize/2, RunDecodedPositions{bestLaps(kk)}, clim) % 
    
    set(ax, 'YDir', 'normal', 'colormap', hot)
    
    
    hold on

    plot((1:noTimeBins)*binDur-binDur/2, posbinIdx{bestLaps(kk)}*posBinSize, '-', 'linewidth', 1, 'color', 'c')
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



PosBinWidth = 10;
bins = 0 : PosBinWidth : noPosBins*posBinSize;

[peakRates, PF_peaks] = max(spatialTunings, [], 2);
PF_peaks(peakRates < 1) = [];

PF_peaks = PF_peaks*posBinSize;

counts = hist(PF_peaks, bins);
% counts = counts./sum(counts);

bar(bins, counts, 'FaceColor', 'k', 'EdgeColor', 'none')

xlabel('Position on track (cm)', 'fontsize', 10)
ylabel('No. of units', 'fontsize', 10)

ylim([0 1.05*max(counts)])

t= title({'(b) Dist. of preferred firing location'; '(Max FR > 1Hz)'; ''}, 'fontsize', 10);
set(t, 'horizontalAlignment', 'center')

set(gca, 'box', 'off')

oldPos = ax2.Position;
ax2.Position = [oldPos(1) oldPos(2) oldPos(3)*0.9 oldPos(4)*0.9];


ax3 = subplot(2,7,[6 7]);


% hold on
% 
% h1 = plot(positioncenters, posBinErrorStat(:, 2), 'k', 'linewidth', 2);
% % plot(positioncenters, posBinErrorStat(:, 1), '--k', 'linewidth', 1)
% % plot(positioncenters, posBinErrorStat(:, 3), '--k', 'linewidth', 1)
% 
% % 
% % h2 = plot(positioncenters, posBinErrorStatShuffle(:, 2), 'color', [.7 .7 .7], 'linewidth', 2);
% % plot(positioncenters, posBinErrorStatShuffle(:, 1), 'color', [.7 .7 .7], 'linestyle', '--', 'linewidth', 1)
% % plot(positioncenters, posBinErrorStatShuffle(:, 3), 'color', [.7 .7 .7], 'linestyle', '--', 'linewidth', 1)
% 
% % 
% % h3 = plot(positioncenters, posBinErrorStatDiff(:, 2), 'r', 'linewidth', 2);
% % plot(positioncenters, posBinErrorStatDiff(:, 1), '--r', 'linewidth', 1)
% % plot(positioncenters, posBinErrorStatDiff(:, 3), '--r', 'linewidth', 1)
% 
% hold off
% 
% xlim([0 positioncenters(end)+largePosBinSize/2])
% ylim([0 max(20, max(posBinErrorStat(:, 2)))])
% 
% ylabel('Median decoding error (cm)', 'fontsize', 10)
% xlabel('Position on track (cm)')


imagesc(decConfMat); 

set(ax3, 'colormap', flipud(bone), 'YDir', 'normal')
xticks([0.5 length(uniqPosBins)+0.5]);
xticklabels({'0', sprintf('%d', positioncenters(end)+largePosBinSize/2)})

yticks([0.5 length(uniqPosBins)+0.5]);
yticklabels({'0', sprintf('%d', positioncenters(end)+largePosBinSize/2)})


xlabel('Decoded position', 'fontsize', 10)
ylabel('Actual position', 'fontsize', 10)

t= title({'(c) Bayesian decoding confusion matrix'; ''}, 'fontsize', 10);
set(t, 'horizontalAlignment', 'center')


% set(gca, 'box', 'off')


oldPos = ax3.Position;
ax3.Position = [oldPos(1) oldPos(2) oldPos(3)*0.9 oldPos(4)*0.9];


oldPos = ax3.Position;
ax3.Position  = [oldPos(1)-oldPos(3)/3 oldPos(2) oldPos(3) oldPos(4)];

axis square
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

elseif strcmp(direction, 'RL')
    
    filename = fullfile(FileBase, 'spatialTuning_RL');
else
    filename = fullfile(FileBase, 'spatialTuning');
end

savepdf(gcf, filename, '-dsvg')
savepdf(gcf, filename, '-dpng')
% savefig(gcf, filename)



end




