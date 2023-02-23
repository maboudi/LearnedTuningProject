
% Null hypothesis: Distribution of peak positions in place field multiplication in independent of which units are active in each time bin 

spatialTunings_LR2 = spatialTunings_LR ./ repmat(sum(spatialTunings_LR, 2), [1 size(spatialTunings_LR, 2)]);
spatialTunings_RL2 = spatialTunings_RL ./ repmat(sum(spatialTunings_RL, 2), [1 size(spatialTunings_RL, 2)]);

spatialTunings_RL2(~spatialTunings_RL) = 0.001;
spatialTunings_LR2(~spatialTunings_LR) = 0.001;



[nUnits, nPosBins] = size(spatialTunings_LR);


periods = {'PRE'; 'RUN'; 'POST'};
colors  = {'b'; 'k'; 'r'};


binnedPBEs = struct(periods{1}, PREbinnedPBEs, periods{2}, RUNbinnedPBEs, periods{3}, POSTbinnedPBEs);

PFmultiPeakPos = struct(periods{1}, [], periods{2}, [], periods{3}, []);
maxofPFMulti   = struct(periods{1}, [], periods{2}, [], periods{3}, []);
sharpness      = struct(periods{1}, [], periods{2}, [], periods{3}, []);


figure;

for period = 1:numel(periods) 

currentPBEs = binnedPBEs.(periods{period}).data(:, 2);
nEvents       = size(currentPBEs, 1);


innerStruct = struct('data', [], 'sim', []);
PFmultiPeakPos.(periods{period}) = innerStruct; PFmultiPeakPos.(periods{period}).data = cell(nEvents, 1); PFmultiPeakPos.(periods{period}).sim = cell(nEvents, 1);
maxofPFMulti.(periods{period})   = innerStruct; maxofPFMulti.(periods{period}).data   = cell(nEvents, 1); maxofPFMulti.(periods{period}).sim   = cell(nEvents, 1); 
sharpness.(periods{period})      = innerStruct; sharpness.(periods{period}).data      = cell(nEvents, 1); sharpness.(periods{period}).sim      = cell(nEvents, 1);

for pbe = 1: nEvents

    currEvent = currentPBEs{pbe};
    nTimebins = size(currEvent, 2);

    pbeActUnits = find(sum(currEvent, 2));

    % calculate new spatial tunings with shuffled cell IDs
    [spatialTunings_RLshuffles, spatialTunings_LRshuffles] = tuningsUnitIDshuffle_normalize(spatialTunings_RL2, spatialTunings_LR2, pbeActUnits, 1, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'


    maxofPFMulti.(periods{period}).data{pbe}   = zeros(nTimebins, 1);
    sharpness.(periods{period}).data{pbe}      = zeros(nTimebins, 1);
    PFmultiPeakPos.(periods{period}).data{pbe} = zeros(nTimebins, 1);

    maxofPFMulti.(periods{period}).sim{pbe}   = zeros(nTimebins, 1);
    sharpness.(periods{period}).sim{pbe}      = zeros(nTimebins, 1);
    PFmultiPeakPos.(periods{period}).sim{pbe} = zeros(nTimebins, 1);


    for tb = 1:nTimebins

        currBinActiveUnits = find(currEvent(:, tb));

        if numel(currBinActiveUnits) > 0

            % actual unit IDs
            pfMultiplyTwoDirections = prod([spatialTunings_LR2(currBinActiveUnits, :) spatialTunings_RL2(currBinActiveUnits, :)], 1); % concatenated the tunings belonging to two running directions

            [maxofPFMulti.(periods{period}).data{pbe}(tb), peakPosTwodirections] = max(pfMultiplyTwoDirections);
            sharpness.(periods{period}).data{pbe}(tb) = max(pfMultiplyTwoDirections / sum(pfMultiplyTwoDirections));


            peakPos = mod(peakPosTwodirections, nPosBins);
            peakPos(~peakPos) = nPosBins;
            PFmultiPeakPos.(periods{period}).data{pbe}(tb) = peakPos;


            % randomized unit IDs
            pfMultiplyTwoDirections = prod([spatialTunings_LRshuffles(currBinActiveUnits, :) spatialTunings_RLshuffles(currBinActiveUnits, :)], 1); % concatenated the tunings belonging to two running directions

            [maxofPFMulti.(periods{period}).sim{pbe}(tb), peakPosTwodirections] = max(pfMultiplyTwoDirections);
            sharpness.(periods{period}).sim{pbe}(tb) = max(pfMultiplyTwoDirections / sum(pfMultiplyTwoDirections));


            peakPos = mod(peakPosTwodirections, nPosBins);
            peakPos(~peakPos) = nPosBins;
            PFmultiPeakPos.(periods{period}).sim{pbe}(tb) = peakPos;

        else
            PFmultiPeakPos.(periods{period}).data{pbe}(tb) = nan;
            maxofPFMulti.(periods{period}).data{pbe}(tb)   = nan;
            sharpness.(periods{period}).data{pbe}(tb)      = nan;

            PFmultiPeakPos.(periods{period}).sim{pbe}(tb)  = nan;
            maxofPFMulti.(periods{period}).sim{pbe}(tb)    = nan;
            sharpness.(periods{period}).sim{pbe}(tb)       = nan;
        end
    end
end



% plot the histogram of peak positions (in relation to middle of the track) 

subplot(3,3,(period-1)*3+1)

barWidth = 8; % number of position bins 

bins = 0:barWidth: floor(nPosBins/2)+1; 
binCenters = bins + barWidth/2;


% data
currDATA = cell2mat(PFmultiPeakPos.(periods{period}).data);
currDATA = floor(nPosBins/2) - abs(currDATA-floor(nPosBins/2));
countsDATA = hist(currDATA, bins);
% countsDATA(end) = [];
countsDATA = countsDATA/sum(countsDATA);
 
bar(binCenters, countsDATA, 'facecolor', colors{period}, 'facealpha', 0.4, 'edgecolor', 'none')


% shuffle
currSIM = cell2mat(PFmultiPeakPos.(periods{period}).sim);
currSIM = floor(nPosBins/2) - abs(currSIM-floor(nPosBins/2));
countsSIM = hist(currSIM, bins);
% countsSIM(end) = [];
countsSIM = countsSIM/sum(countsSIM);

hold on
bar(binCenters, countsSIM, 'facecolor', 'none', 'edgecolor', 'k')

line([nanmedian(currDATA) nanmedian(currDATA)], [0 max([countsDATA countsSIM])], 'color', colors{period}, 'linewidth', 2)
line([nanmedian(currSIM) nanmedian(currSIM)], [0 max([countsDATA countsSIM])], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

ylim([0 0.3])
set(gca, 'xtick', [4 52], 'xticklabel', {'track ends'; 'track center'}, 'box', 'off')

xlabel('position')
ylabel({periods{period};'';'probability'})

if period == 1
    
    title({'PF multiplication';'peak position'})
end

% plot the histogram of maxofPFMulti (max posterior probability)

subplot(3,3, (period-1)*3+2)


currDATA = cell2mat(maxofPFMulti.(periods{period}).data);
currSIM = cell2mat(maxofPFMulti.(periods{period}).sim);

% bins = 0:0.1:max([currDATA; currSIM]);
bins = linspace(prctile([currDATA; currSIM], 5), prctile([currDATA; currSIM], 95), 10);
binCenters = bins(1:end-1) + (bins(2)-bins(1))/2;

countsDATA = hist(currDATA, bins); countsDATA(end)=[]; %countsDATA = countsDATA/sum(countsDATA);
countsSIM = hist(currSIM, bins); countsSIM(end)=[]; %countsSIM = countsSIM/sum(countsSIM);

countsDATA = log10(countsDATA);
countsSIM = log10(countsSIM);


bar(binCenters, countsDATA, 'faceColor', colors{period}, 'edgeColor', 'none', 'facealpha', 0.5)
hold on
bar(binCenters, countsSIM, 'faceColor', 'none', 'edgeColor', 'k', 'facealpha', 0.5)

line([nanmedian(currDATA) nanmedian(currDATA)], [0 max([countsDATA countsSIM])], 'color', colors{period}, 'linewidth', 2)
line([nanmedian(currSIM) nanmedian(currSIM)], [0 max([countsDATA countsSIM])], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

xlabel('PF multiplication peak prob.')
ylabel('probability')
set(gca, 'box', 'off')

if period == 1
    title({'overlap'; 'of place fields'})
end

% plot the histogram of maxofPFMulti (max posterior probability)

subplot(3,3, period*3)

% bins = 0:0.1:max([sharpness.(periods{period}).data sharpness.(periods{period}).sim]);
bins = 0:0.1:0.7;
binCenters = bins(1:end-1) + 0.05;

% data
currDATA = cell2mat(sharpness.(periods{period}).data);
countsDATA = hist(currDATA, bins); countsDATA(end)=[]; countsDATA = countsDATA/sum(countsDATA);

bar(binCenters, countsDATA, 'faceColor', colors{period}, 'edgeColor', 'none', 'facealpha', 0.5)

% shuffle
currSIM = cell2mat(sharpness.(periods{period}).sim);
countsSIM = hist(currSIM, bins); countsSIM(end)=[]; countsSIM = countsSIM/sum(countsSIM);
hold on
bar(binCenters, countsSIM, 'faceColor', 'none', 'edgeColor', 'k', 'facealpha', 0.5)

line([nanmedian(currDATA) nanmedian(currDATA)], [0 max([countsDATA countsSIM])], 'color', colors{period}, 'linewidth', 2)
line([nanmedian(currSIM) nanmedian(currSIM)], [0 max([countsDATA countsSIM])], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

xlabel('sharpness')
ylabel('probability')
set(gca, 'box', 'off')

if period == 1
   title({'sharpness'; 'of place field multiplication'})
end

end


% nBins = numel(bins)-1;

% if mod(nBins, 2) == 0
%    countsDATA = sum([countsDATA(1:nBins/2); countsDATA(nBins:-1:nBins/2+1)]);
%    
% else
%    countsDATA = [sum([countsDATA(1:floor(nBins/2)); countsDATA(nBins:-1:floor(nBins/2)+2)]) countsDATA(ceil(nBins/2))];
% end

% binCenters = bins(1:ceil(nBins/2));


% medd.data = find(cumsum(countsDATA)/sum(countsDATA) > 0.5, 1, 'first');
% medd.sim  = find(cumsum(countsSIM)/sum(countsSIM) > 0.5, 1, 'first');