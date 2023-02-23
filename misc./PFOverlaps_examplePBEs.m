
spatialTunings_LR = spatialTunings_LR./repmat(sum(spatialTunings_LR, 2), [1, size(spatialTunings_LR, 2)]);
spatialTunings_RL = spatialTunings_RL./repmat(sum(spatialTunings_RL, 2), [1, size(spatialTunings_RL, 2)]);


spatialTunings_LR(~spatialTunings_LR) = 0.001;
spatialTunings_RL(~spatialTunings_RL) = 0.001;


seNum = 6347; % sample event/PBE number
selectTimeBin1 = 7;
selectTimeBin2 = 15;

postPr = posteriorProbMatrix.RUN.data{seNum,1};
[nPosBins, nTimeBins] = size(postPr);


spikeRaster = RUNbinnedPBEs.data{seNum, 1};
nUnits = size(spikeRaster, 1);

convWin = gausswindow(2,2);

spikeRasterDisp = spikeRaster;
for ii = 1:size(spikeRaster, 2)
        
    temp = conv(spikeRaster(:, ii), convWin, 'same'); 
    temp(temp > 0) = 1;
    spikeRasterDisp(:, ii) = temp;
    
end

cm1 = colormap('jet');
cm2 = colormap('bone');
cm2 = flipud(cm2);

figure;
set(gcf, 'units', 'centimeters', 'position', [10 10 20 30])



%% Posterior probability

ax1 = subplot(4,2,[1 2]);

imagesc(postPr);

set(ax1, 'Ydir', 'normal', 'xtick', [1 nTimeBins], 'ytick', [1 nPosBins], 'fontsize', 14)
% xlabel('PBE time bin')
ylabel('position')
title('posterior probability matrix')
ax1.Colormap = cm1;




%% spike raster

ax2 = subplot(4,2,[3 4]);


imagesc(spikeRasterDisp);


set(ax2, 'yDir', 'normal', ...
    'xtick', [10  (nTimeBins-1)*20+10], 'xticklabel', {'1'; sprintf('%d', nTimeBins)}, ...
    'ytick', [1 nUnits], 'fontsize', 14)



for ii = 1: nTimeBins-1
    
   line([ii*20 ii*20], [1 nUnits], 'linestyle', '-', 'color', 'b', 'linewidth', 0.5);
    
end

xlabel('PBE time bin')
ylabel('unit')
title('spike ratser')

ax2.Colormap = cm2;

hold on
patch([(selectTimeBin1-1)*20+1 selectTimeBin1*20 selectTimeBin1*20 (selectTimeBin1-1)*20+1], [1 1 108 108], 'b', 'facealpha', 0.3, 'edgeColor', 'none')

hold on
patch([(selectTimeBin2-1)*20+1 selectTimeBin2*20 selectTimeBin2*20 (selectTimeBin2-1)*20+1], [1 1 108 108], 'g', 'facealpha', 0.3, 'edgeColor', 'none')



%% place fields of exapmple

% selectTimeBin 1

firingUnits = find(RUNbinnedPBEs.data{seNum, 2}(:, selectTimeBin1));

ax3 = subplot(4,2, 5);
ax3.YGrid = 'on';

poss = ax3.Position;

hold on

for ii = 1: numel(firingUnits)
    if max(spatialTunings_LR(firingUnits(ii), :)) > 0.001
        plot(1:nPosBins, spatialTunings_LR(firingUnits(ii), :)/max(spatialTunings_LR(firingUnits(ii), :)) + ii, 'linewidth', 1 )
    else
        plot(1:nPosBins, spatialTunings_LR(firingUnits(ii), :)+ ii, 'linewidth', 1 )
    end
    
end

summedPlaceFields = sum(spatialTunings_LR(firingUnits, :));
prodPlaceFields   = prod(spatialTunings_LR(firingUnits, :));


plot(1: nPosBins, summedPlaceFields/max(summedPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '-');

plot(1: nPosBins, prodPlaceFields/max(prodPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '--');

plot(1: nPosBins, postPr(:, selectTimeBin1)/max(postPr(:, selectTimeBin1))*numel(firingUnits), 'color', [.7 .7 .7], 'linewidth', 1);


set(ax3, 'ytick', 1:numel(firingUnits))
xlabel('position bin')
ylabel({'template 1'; ''; 'firing unit'})

title('\color{blue}example time bin 1')


ax4 = subplot(4,2, 7);
ax4.YGrid = 'on';


hold on

for ii = 1: numel(firingUnits)
    if max(spatialTunings_RL(firingUnits(ii), :)) > 0.001
        plot(1:nPosBins, spatialTunings_RL(firingUnits(ii), :)/max(spatialTunings_RL(firingUnits(ii), :)) + ii, 'linewidth', 1 )
    else
        plot(1:nPosBins, spatialTunings_RL(firingUnits(ii), :) + ii, 'linewidth', 1 )
    end
    
end

summedPlaceFields = sum(spatialTunings_RL(firingUnits, :));
prodPlaceFields   = prod(spatialTunings_RL(firingUnits, :));

plot(1: nPosBins, summedPlaceFields/max(summedPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '-');

plot(1: nPosBins, prodPlaceFields/max(prodPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '--');

plot(1: nPosBins, postPr(:, selectTimeBin1)/max(postPr(:, selectTimeBin1))*numel(firingUnits), 'color', [.7 .7 .7], 'linewidth', 1);




set(ax4, 'ytick', 1:numel(firingUnits))
xlabel('position bin')
ylabel({'template 2'; ''; 'firing unit'})


% selectTimeBin2


firingUnits = find(RUNbinnedPBEs.data{seNum, 2}(:, selectTimeBin2));

ax5 = subplot(4,2, 6);
ax5.YGrid = 'on';
poss = ax5.Position;

hold on

for ii = 1: numel(firingUnits)
    
    if max(spatialTunings_LR(firingUnits(ii), :)) > 0.001
        plot(1:nPosBins, spatialTunings_LR(firingUnits(ii), :)/max(spatialTunings_LR(firingUnits(ii), :)) + ii, 'linewidth', 1 )
    else
        plot(1:nPosBins, spatialTunings_LR(firingUnits(ii), :) + ii, 'linewidth', 1 )
    end
    
end

summedPlaceFields = sum(spatialTunings_LR(firingUnits, :));
prodPlaceFields   = prod(spatialTunings_LR(firingUnits, :));

h1 = plot(1: nPosBins, summedPlaceFields/max(summedPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '-');

h2 = plot(1: nPosBins, prodPlaceFields/max(prodPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '--');

h3 = plot(1: nPosBins, postPr(:, selectTimeBin2)/max(postPr(:, selectTimeBin2))*numel(firingUnits), 'color', [.7 .7 .7], 'linewidth', 1);


legend([h1, h2, h3], 'summation', 'multiplication', 'posterior prob', 'location', 'northout')
ax5.Position = poss;


set(ax5, 'ytick', 1:numel(firingUnits))
xlabel('position bin')
ylabel('firing unit')

title('\color{green}example time bin 2')


ax6 = subplot(4,2, 8);
ax6.YGrid = 'on';

hold on

for ii = 1: numel(firingUnits)
    if max(spatialTunings_RL(firingUnits(ii), :)) > 0.001
        plot(1:nPosBins, spatialTunings_RL(firingUnits(ii), :)/max(spatialTunings_RL(firingUnits(ii), :)) + ii, 'linewidth', 1 )
    else
        plot(1:nPosBins, spatialTunings_RL(firingUnits(ii), :) + ii, 'linewidth', 1 )
    end
    
end

summedPlaceFields = sum(spatialTunings_RL(firingUnits, :));
prodPlaceFields   = prod(spatialTunings_RL(firingUnits, :));

plot(1: nPosBins, summedPlaceFields/max(summedPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '-')

plot(1: nPosBins, prodPlaceFields/max(prodPlaceFields)*numel(firingUnits), 'color', 'k', 'linewidth', 2, 'linestyle', '--')

plot(1: nPosBins, postPr(:, selectTimeBin2)/max(postPr(:, selectTimeBin2))*numel(firingUnits), 'color', [.7 .7 .7], 'linewidth', 1)


set(ax6, 'ytick', 1:numel(firingUnits))
xlabel('position bin')
ylabel('firing unit')



%% save the figure

% saveas(gcf, ['examplePBEPlaceFieldAnalysis_' num2str(seNum)],  'epsc')
print(gcf, ['examplePBEPlaceFieldAnalysis_' num2str(seNum)],  '-dpdf')

close all



