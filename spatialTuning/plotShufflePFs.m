
figure;

x0=0;
y0=0;
width=800;
height=400* nUnits/20;
set(gcf,'units','points','position',[x0,y0,width,height])



ax1 = subplot(1,2,1); plotPlaceFields(spatialTunings_LR, spatialTunings_LRshuffles(:,:,5), 'r')
ax2 = subplot(1,2,2); plotPlaceFields(spatialTunings_RL, spatialTunings_RLshuffles(:,:,5), 'b')


function plotPlaceFields(spatialTunings, spatialTunings_shuffle, cl)


[nUnits, nPosBins] = size(spatialTunings);
posBinSize = 2; % cm

xAxis = (1:nPosBins)*posBinSize; 


[peakRates, peakPos] = max(spatialTunings, [], 2);

[~, sortIdx] = sort(peakPos, 'ascend');

spatialTunings_sorted = spatialTunings(sortIdx, :);
spatialTunings_shuffle_sorted = spatialTunings_shuffle(sortIdx, :);

peakRates_sorted = peakRates(sortIdx);

spatialTunings_sorted_norm  = spatialTunings_sorted./repmat(max(spatialTunings_sorted ,[], 2), [1 nPosBins]); 
spatialTunings_shuffle_norm = spatialTunings_shuffle_sorted./repmat(max(spatialTunings_shuffle_sorted ,[], 2), [1 nPosBins]); 

shuffleRates_sorted = max(spatialTunings_shuffle_sorted ,[], 2);
shuffleRates_sorted = shuffleRates_sorted(sortIdx);

tt = 0;

for jj = 1 : size(spatialTunings_sorted_norm, 1)
    
    if peakRates_sorted(jj) > 1
        
        tt = tt + 1;
%         cl = 'r';
        fill([xAxis fliplr(xAxis)], [0.06*tt+spatialTunings_sorted_norm(jj, :)/20 fliplr(0.06*tt*ones(size(spatialTunings_sorted_norm(jj, :))))], cl,'LineStyle','none')
        hold on
        
        plot(xAxis, 0.06*tt+spatialTunings_sorted_norm(jj, :)/20,'color', 'k','linewidth', 0.5);
        
        
        fill([xAxis fliplr(xAxis)], [0.06*tt+spatialTunings_shuffle_norm(jj, :)/20 fliplr(0.06*tt*ones(size(spatialTunings_shuffle_norm(jj, :))))], [.5 .5 .5],'LineStyle','none')
        hold on
        plot(xAxis, 0.06*tt+spatialTunings_shuffle_norm(jj, :)/20,'color', 'k','linewidth', 0.5, 'linestyle', ':');
        
        
        alpha(0.5)
        
    end

    set(gca, 'YTick', [], 'YTickLabel', [], 'color', 'none', 'YColor', 'none', 'box', 'off')
end

xlim([xAxis(1)-posBinSize/2 xAxis(end)+posBinSize/2])

xlabel('Position on track (cm)', 'fontsize', 10)

h = text(xAxis(1)-5*posBinSize, 0.06*(tt)/2, 'Unit', 'fontsize', 10, 'HorizontalAlignment', 'center');
set(h, 'rotation', 90)

end


