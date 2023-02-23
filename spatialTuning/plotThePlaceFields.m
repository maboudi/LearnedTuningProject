function plotThePlaceFields(placeTunings, spikePositions, spikeUnit, activeUnits, fileinfo, behavior, direction, posBinSize, FileBase)

[noYbins, noXbins] = size(placeTunings(:, :, 1));
noUnits = length(activeUnits);

nCol = 4;
nRo  = ceil(noUnits/nCol);
% nRo = ceil(noXbins/noYbins * nCol);
% posBinSize = 2;

if noUnits > nCol*nRo
    
    rr = noUnits/nCol/nRo;
    nCol = ceil(nCol * rr);
    nRo = ceil(nRo * rr);
    
end



runIdx = find(fileinfo.xyt(:,3) > behavior.time(2,1) & fileinfo.xyt(:,3)< behavior.time(2,2));

runXpose = fileinfo.xyt(runIdx, 1);
runYpose = fileinfo.xyt(runIdx, 2);

graphAspectRatio = (noXbins/noYbins)/(nRo/nCol);

figure;
set(gcf, 'position', [100 100 1200*graphAspectRatio 1200])

for ii = 1:noUnits
   
    subplot(nRo, nCol, ii)
 
    unitSpikePositions = spikePositions(spikeUnit ==  activeUnits(ii), :);

    plot(runXpose-min(runXpose), runYpose-min(runYpose), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])

    hold on

    plot(unitSpikePositions(:, 2)-min(runXpose), unitSpikePositions(:, 1)-min(runYpose), '.', 'markersize', 3, 'color', 'b')
    
    xlim([0 max(runXpose-min(runXpose))])
    ylim([0 max(runYpose-min(runYpose))])
    
    set(gca, 'xtick',[],'ytick',[])
    title(sprintf('u %d', activeUnits(ii)), 'fontsize', 5, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
    
end


savepdf(gcf, fullfile(FileBase, [fileinfo.name '_firingDistribution_' direction]), '-dpng')
% saveas(gcf, fullfile(FileBase, [fileinfo.name '_firingDistribution_' direction]), 'epsc')


figure;
set(gcf, 'position', [100 100 1200*graphAspectRatio 1200])

placeTunings = placeTunings(:,:, activeUnits);

for ii = 1:noUnits
   
    subplot(nRo, nCol, ii)
    
    currTuning = placeTunings(:, :, ii);
    
    plot(runXpose-min(runXpose), runYpose-min(runYpose), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
    
    hold on
    
    imagesc((1:noXbins)*posBinSize, (1:noYbins)*posBinSize, currTuning)
    set(gca, 'YDir', 'normal')
    xlim([0 noXbins*posBinSize])
    ylim([0 noYbins*posBinSize])
    
    alpha(0.7) 
    
    
    set(gca, 'xtick',[],'ytick',[])
    title(sprintf('u %d, %.1f Hz', activeUnits(ii), nanmax(currTuning(:))), 'fontsize', 5, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')
   
end

colormap('jet')


savepdf(gcf, fullfile(FileBase, [fileinfo.name '_2DplaceFields_' direction]), '-dpng')
% saveas(gcf, fullfile(FileBase, [fileinfo.name '_2DplaceFields_' direction]), 'epsc')



figure; 
set(gcf, 'position', [0 0 300 300*noYbins/noXbins])

MUATuning = nanmean(placeTunings, 3);
plot(runXpose-min(runXpose), runYpose-min(runYpose), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])

hold on

imagesc((1:noXbins)*posBinSize, (1:noYbins)*posBinSize, MUATuning)

set(gca, 'YDir', 'normal')

xlim([0 noXbins*posBinSize])
ylim([0 noYbins*posBinSize])

alpha(0.7) 

set(gca, 'xtick',[0 noXbins*posBinSize],'ytick',[0 noYbins*posBinSize])
% xlabel('X(cm)', 'fontsize', 8)
% ylabel('Y(cm)', 'fontsize', 8)
% title('multiunit place field', 'fontsize', 8, 'Units', 'normalized', 'Position', [0 1], 'HorizontalAlignment', 'left')

colormap('jet')

savepdf(gcf, fullfile(FileBase, [fileinfo.name '_2D_MUA_placeField_' direction]), '-dpng')
% saveas(gcf,fullfile(FileBase, [fileinfo.name '_2D_MUA_placeField_' direction]), 'epsc')

end
    
    