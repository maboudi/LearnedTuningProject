
close all
clc


spikes = spikes_pyr;
nUnits = numel(spikes);

peakPosBins = nan(nUnits, 1);
for iUnit = 1:nUnits
    peakPosBins(iUnit) = spikes(iUnit).peakPosBin.uni;
end

[~, sortIdx] = sort(peakPosBins, 'ascend');


PBEInfo = PBEInfo_replayScores;

pbe = 1612; 
selectTimeBins = [4 8];


displayFirstBin = 1;
displayLastBin  = size(PBEInfo(pbe).posteriorProbMat, 2)-1;

postPr = PBEInfo(pbe).posteriorProbMat(:, displayFirstBin:displayLastBin);


[nPosBins, nTimeBins] = size(postPr);

spikeRaster = PBEInfo(pbe).fr_1msbin(sortIdx, 20*(displayFirstBin-1)+1:20*displayLastBin);
firing_20ms = PBEInfo(pbe).fr_20msbin(sortIdx, displayFirstBin:displayLastBin);


convWin = gausswindow(2,2);

spikeRasterDisp = spikeRaster;
for ii = 1:size(spikeRaster, 2)       
    temp = conv(spikeRaster(:, ii), convWin, 'same'); 
    temp(temp > 0) = 1;
    spikeRasterDisp(:, ii) = temp;  
end
% spikeRasterDisp = spikeRaster;


cm1 = colormap(flipud(copper));
cm2 = colormap('bone');
cm2 = flipud(cm2);



plotheight = 60;
plotwidth  = 30;
fontsize   = 2;

f= figure;
clf(f);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



%% Posterior probability

ax1 = subplot(6,2,[1 2]);

imagesc(1:nTimeBins, 1:nPosBins*2, postPr);
caxis([0 0.05])

set(ax1, 'Ydir', 'normal', 'xtick', [], 'ytick', [0 2*nPosBins], 'YTickLabel', [0 1], 'fontsize', 4)
% xlabel('PBE time bin')
ylabel('norm. position')
title(sprintf('%s: ripple # %d', upper(PBEInfo(pbe).epoch), pbe), 'fontsize', fontsize, 'FontWeight','normal')
ax1.Colormap = cm1;
set(gca, 'box', 'off', 'linewidth', 0.5, 'fontsize', 4, 'TickDir', 'out','TickLength',[0.01, 0.01]) 

oldPos = get(gca, 'Position');
newPos = [oldPos(1) oldPos(2)-0.75*oldPos(4) oldPos(3) 1.5*oldPos(4)];

set(gca, 'position', newPos)



%% spike raster

ax2 = subplot(6,2,[5 6]);

hold on
% imagesc(spikeRasterDisp);

for iUnit = 1:nUnits

    timeStamps = find(spikeRaster(iUnit, :));

    for spk = 1:numel(timeStamps)
        plot([timeStamps(spk) timeStamps(spk)], ...
            [iUnit-0.45 iUnit+0.45], 'color', 'k', 'linewidth', 0.25); % cl(ceil(interp1(19:141, linspace(0,11, 141-19+1), peakPos(iUnit))), :)
    end
end


xlim([0 nTimeBins*20])
ylim([-2 nUnits])

set(ax2, 'yDir', 'normal', ...
    'XTick', [], ...
    'ytick', [1 nUnits], 'fontsize', fontsize)



for ii = 1: nTimeBins-1
   line([ii*20 ii*20], [1 nUnits], 'linestyle', '-', 'color', [0.7 0.7 0.7 0.5], 'linewidth', 0.25);
end

% xlabel('PBE time bin', 'fontsize', fontsize)
ylabel('units', 'fontsize', fontsize)

ax2.Colormap = cm2;

hold on
patch([(selectTimeBins(1)-1)*20 selectTimeBins(1)*20 selectTimeBins(1)*20 (selectTimeBins(1)-1)*20+1], [1 1 nUnits nUnits], 'k', 'FaceColor', 'none', 'facealpha', 0.3, 'edgeColor', '#F9A603', 'EdgeAlpha', 0.9, 'linewidth', 0.5)

hold on
patch([(selectTimeBins(2)-1)*20 selectTimeBins(2)*20 selectTimeBins(2)*20 (selectTimeBins(2)-1)*20+1], [1 1 nUnits nUnits], 'k', 'FaceColor', 'none', 'facealpha', 0.3, 'edgeColor', '#DEB887', 'EdgeAlpha', 0.9, 'linewidth', 0.5)

set(gca, 'box', 'off', 'linewidth', 0.5, 'fontsize', 4, 'TickDir', 'out', 'TickLength', [0.01, 0.01]) 


oldPos = get(gca, 'Position');
newPos = [oldPos(1) oldPos(2)+0.2*oldPos(4) oldPos(3) 1.5*oldPos(4)];

set(gca, 'position', newPos)


%% place fields within representative time bins 

max_N_firingUnits = 10;

for iBin = 1:numel(selectTimeBins)

    currTimeBin = selectTimeBins(iBin);

    firingUnits = find(firing_20ms(:, currTimeBin));
    

    % plot the place fields
    if iBin == 1
        curr_ax = subplot(6,2, [9 11]);
    elseif iBin == 2
        curr_ax = subplot(6,2, [10 12]);
    end

    hold on

    currSpatialTuning = nan(numel(firingUnits), nPosBins);
    for ii = 1: numel(firingUnits)
    
        currSpatialTuning(ii, :) = spikes(sortIdx(firingUnits(ii))).spatialTuning_smoothed.uni + 1e-4;

        if max(currSpatialTuning(ii, :)) > 0.001
            plot(1:nPosBins, currSpatialTuning(ii, :)/max(currSpatialTuning(ii, :))*0.9+ max_N_firingUnits-ii+1, 'linewidth', 0.25, 'Color', 'k')
        else
            plot(1:nPosBins, currSpatialTuning(ii, :)*0.9+max_N_firingUnits-ii+1, 'linewidth', 0.5, 'Color', 'k')
        end 
    end

    plot(1: nPosBins, postPr(:, currTimeBin)/max(postPr(:, currTimeBin))+ max_N_firingUnits+1.5, 'color', 'r', 'linewidth', 0.25);

    ylim([0 max_N_firingUnits+3])

    if iBin == 1
        curr_ax.XColor = '#F9A603';
        curr_ax.YColor = '#F9A603';
    else
        curr_ax.XColor = '#DEB887';
        curr_ax.YColor = '#DEB887';
    end
    

    if iBin == 1
        xlabel('position', 'fontsize', fontsize, 'color', 'k')
        ylabel('participant units', 'fontsize', fontsize, 'Color', 'k')
    end
    set(curr_ax, 'box', 'on', 'xtick', [], 'xticklabel', [], 'ytick', [], 'linewidth', 0.5, 'fontsize', 4, 'TickDir', 'out','TickLength',[0.01, 0.01]) 


    oldPos = get(curr_ax, 'Position');
    newPos = [oldPos(1) oldPos(2)+0.15*oldPos(4) oldPos(3) oldPos(4)];
    
    set(curr_ax, 'position', newPos)

end


%%
% correlation between the place field of acive units within each time bin
% and the decoded posterior within the same time bin 


meanCorr_w_posterior = nan(nTimeBins, 1);

for iBin = 1:nTimeBins

    firingUnits = find(firing_20ms(:, iBin));
    
    currPosterior = postPr(:, iBin);

    currSpatialTuning = nan(numel(firingUnits), nPosBins);
    if numel(firingUnits) > 1
        for ii = 1: numel(firingUnits) 
            currSpatialTuning(ii, :) = spikes(sortIdx(firingUnits(ii))).spatialTuning_smoothed.uni + 1e-4;
        end
        
        allCorrs = corr(currSpatialTuning', currPosterior); %

%         allCorrs = setdiff(allCorrs(:), diag(allCorrs));
        meanCorr_w_posterior(iBin) = nanmean(allCorrs);

    else
        continue
    end

end

ax3 = subplot(6,2,[7 8]);


cm = redblue(); %othercolor('OrRd7'); % 
cm = cm(128:256, :);
imagesc(1:nTimeBins, 1, meanCorr_w_posterior');
caxis([0 0.5])
colormap(ax3, cm)

set(ax3, 'Ydir', 'normal', 'xtick', [], 'ytick', [], 'fontsize', 4, 'linewidth', 0.5)% [10  (nTimeBins-1)*20+10], 'xticklabel', {'1'; sprintf('%d', nTimeBins)}
xlabel('time bins', 'fontsize', fontsize)

oldPos = get(ax3, 'Position');
newPos = [oldPos(1) oldPos(2)+0.8*oldPos(4) oldPos(3) 0.5*oldPos(4)];

set(ax3, 'position', newPos)


