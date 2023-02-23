

% place fields



% RL
[~, peakPositionRL] = max(spatialTunings_RL, [], 2);
peakPositionRL(peakPositionRL == 1) = 0;

% LR
[~, peakPositionLR] = max(spatialTunings_LR, [], 2);
peakPositionLR(peakPositionLR == 1) = 0;

nPosBins = size(spatialTunings_LR, 2);
peakPosition = struct('RL', peakPositionRL, 'LR', peakPositionLR);



% PBE time bins

binnedPBEs.PRE  = PREbinnedPBEs;
binnedPBEs.RUN  = RUNbinnedPBEs;
binnedPBEs.POST = POSTbinnedPBEs;

baseStruct   = struct('PRE', [], 'RUN', [], 'POST', []);
[conctFirings, nUnitsFiring, nCnctTimeBins, fieldRange, medianFieldrange, counts] = deal(baseStruct); 

f = fieldnames(baseStruct);

bins = 0:5:nPosBins+5;
% bins = 0:5:70+5;
binCenters = bins(1:end-1)+2.5;

for ii = 1: length(f)

% pool all of the 20 ms time bins
conctFirings.(f{ii}) = cell2mat(binnedPBEs.(f{ii}).data(:, 2)'); 


% calculate number of units firing in a time bin
conctFirings.(f{ii})(conctFirings.(f{ii})~=0) = 1;
nUnitsFiring.(f{ii}) = sum(conctFirings.(f{ii}), 1);

% exlude all of the silent time bins
conctFirings.(f{ii})(:, nUnitsFiring.(f{ii}) <= 1) = [];
nUnitsFiring.(f{ii})(nUnitsFiring.(f{ii}) <= 1)   = [];


nCnctTimeBins.(f{ii}) = size(conctFirings.(f{ii}), 2);

fieldRange.(f{ii}) = struct('RL', zeros(nCnctTimeBins.(f{ii}), 1), ...
    'LR', zeros(nCnctTimeBins.(f{ii}), 1), 'pooled', []);
medianFieldRange.(f{ii}) = struct('RL', [], 'LR', [], 'pooled', []);

ff = fieldnames(fieldRange.(f{ii}));

for direction = 1:length(ff)
    for tb = 1: nCnctTimeBins.(f{ii})
        try
            
        fieldRange.(f{ii}).(ff{direction})(tb) = range(peakPosition.(ff{direction})(find(conctFirings.(f{ii})(:, tb))));
%         fieldRange.(f{ii}).(ff{direction})(tb) = std(peakPosition.(ff{direction})(find(conctFirings.(f{ii})(:, tb))));

        catch
            fieldRange.(f{ii}).(ff{direction})(tb) = 0;
        end
    end
    
    medianFieldrange.(f{ii}).(ff{direction}) = median(fieldRange.(f{ii}).(ff{direction}));
end


fieldRange.(f{ii}).pooled = [fieldRange.(f{ii}).RL; fieldRange.(f{ii}).LR];
medianFieldRange.(f{ii}).Pooled = median(fieldRange.(f{ii}).pooled);

counts.(f{ii}) = hist(fieldRange.(f{ii}).pooled, bins);
counts.(f{ii})(end) = []; 
counts.(f{ii}) = counts.(f{ii})/sum(counts.(f{ii})); 


end


cnct = [counts.PRE; counts.RUN; counts.POST];
% cnctexcludefirstbin = [counts.PRE(2:end); counts.RUN(2:end); counts.POST(2:end)];


% yAxisBreakStart   = max(cnctexcludefirstbin(:)) * 1.1;
% yAxisBreakEnd = min(cnct(:, 1)) * 0.9;

% 
% widthtogetridof = diff([yAxisBreakStart yAxisBreakEnd]); 
% counts2 = counts;
% for ii = 1: length(f)
%     
%     idx = find(counts.(f{ii}) > yAxisBreakStart);
%     counts2.(f{ii})(idx) = counts.(f{ii})(idx) - widthtogetridof;
%     
% end



figure;
set(gcf, 'units', 'centimeters', 'position', [0 0 9 15])

ax1 = subplot(5,1,1:2);

h1 = bar(binCenters, counts.PRE, 'facecolor', 'b', 'edgecolor', 'none', 'facealpha', .4);
hold on
h2 = bar(binCenters, counts.POST, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .4);
hold on
h3 = bar(binCenters, counts.RUN, 'facecolor', 'none', 'edgecolor', 'k');

% regularYTicks = 0: 0.03: (max(cnct(:))- widthtogetridof);
regularYTicks = 0: 0.03: max(cnct(:));

yTicklabels   = regularYTicks;

% yTicklabels(yTicklabels > yAxisBreakStart) = yTicklabels(yTicklabels > yAxisBreakStart) + widthtogetridof;

yTicklabels = ceil(yTicklabels*100)/100;

set(gca, 'ytick', regularYTicks, 'yticklabel', num2str(yTicklabels'), 'xticklabel', [])

% xlabel('range of pref. position', 'fontsize', 14)
ylabel('probability', 'fontsize', 14)
% 
% rr = xlim;
% patch([rr(1) rr(2) rr(2) rr(1)], ...
%     [yAxisBreakStart-0.01 yAxisBreakStart-0.01 yAxisBreakStart+0.01 yAxisBreakStart+0.01], 'w', ...
% 'EdgeColor', 'none')


line([median(fieldRange.PRE.pooled) median(fieldRange.PRE.pooled)], [0 max(cnct(:))], 'color', 'b', 'linewidth', 2)
line([median(fieldRange.POST.pooled) median(fieldRange.POST.pooled)], [0 max(cnct(:))], 'color', 'r', 'linewidth', 2)
line([median(fieldRange.RUN.pooled) median(fieldRange.RUN.pooled)], [0 max(cnct(:))], 'color', 'k', 'linewidth', 2)



legend([h1, h2, h3], 'PRE', 'POST', 'RUN', 'location', 'northwest', 'box', 'off')


ax2 = subplot(5,1,3:5);

hold on

colors = {'b', 'k', 'r'};
shifts = [-0.2 0 0.2];
for ii = 1:length(f)
    
    uniqNofFiringUnits = unique(nUnitsFiring.(f{ii}));
    numofnFiringUnits  = length(uniqNofFiringUnits);
    
    for jj = 1: numofnFiringUnits
       
        CurrFieldRanges = [fieldRange.(f{ii}).RL(nUnitsFiring.(f{ii}) == uniqNofFiringUnits(jj)); fieldRange.(f{ii}).LR(nUnitsFiring.(f{ii}) == uniqNofFiringUnits(jj))];
        
        h1 = line([prctile(CurrFieldRanges, 25) prctile(CurrFieldRanges, 75)], [uniqNofFiringUnits(jj) uniqNofFiringUnits(jj)]+shifts(ii), 'linewidth', 3, 'color', colors{ii});
        h1.Color = [h1.Color 0.4];
        plot(median(CurrFieldRanges), uniqNofFiringUnits(jj)+shifts(ii), '.', 'markersize', 10, 'color', colors{ii})
    
    end
    
end
        
% legend([h1.PRE h1.POST h1.RUN], 'PRE', 'RUN', 'POST', 'location', 'northwest')
% plot(fieldRange.PRE.pooled, [nUnitsFiring.PRE nUnitsFiring.PRE], '.', 'markersize', 3, 'color', 'b')

set(gca, 'xtick', 0:10:100, 'xticklabel', num2str((0:10:100)'*2))

xlabel('sd of preferred position (cm)', 'fontsize', 14)
ylabel('number of firing units', 'fontsize', 14)

linkaxes([ax1 ax2], 'x')

dim = [0 0.94 0.5 .1];

textOnTop = {'SD of preferred positions of co-firing units in individual time bins'; ...
    '(All time bins with more than one firing unit were included)'; ...
    'Session name: Achilles_10252013'};
annotation('textbox',dim,'String',textOnTop, 'FitBoxToText','on', 'Interpreter', 'latex');

print(gcf, 'filename3', '-dpdf', '-r0')


print(gcf, 'rangeofPreferredLocations', '-dpdf')



