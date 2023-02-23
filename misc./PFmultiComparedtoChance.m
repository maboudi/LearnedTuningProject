
% Null hypothesis: Distribution of peak positions in place field multiplication in independent of which units are active in each time bin 


spatialTunings_LR2 = spatialTunings_LR ./ repmat(max(spatialTunings_LR, [], 2), [1 size(spatialTunings_LR, 2)]);
spatialTunings_RL2 = spatialTunings_RL ./ repmat(max(spatialTunings_RL, [], 2), [1 size(spatialTunings_RL, 2)]);

spatialTunings_LR2(~spatialTunings_LR) = 0.001;
spatialTunings_RL2(~spatialTunings_RL) = 0.001;


[nUnits, nPosBins] = size(spatialTunings_LR2);


periods = {'PRE'; 'RUN'; 'POST'};
colors  = {'b'; 'k'; 'r'};


binnedPBEs = struct(periods{1}, PREbinnedPBEs, periods{2}, RUNbinnedPBEs, periods{3}, POSTbinnedPBEs);

PFmultiPeakPos = struct(periods{1}, [], periods{2}, [], periods{3}, []);
maxofPFMulti   = struct(periods{1}, [], periods{2}, [], periods{3}, []);
sharpness      = struct(periods{1}, [], periods{2}, [], periods{3}, []);


figure;

for period = 1:numel(periods) 

    
% multiplication of the place fields

cnctActualTimeBins = cell2mat(binnedPBEs.(periods{period}).data(:, 2)');
cnctActualTimeBins(cnctActualTimeBins > 0) = 1; % convert to binary elements: firing or not firing

totNumTimeBins = size(cnctActualTimeBins, 2); % total number of the time bins in the concatenated PBEs 

numActiveUnits = sum(cnctActualTimeBins, 1); % number of active units in individual time bins. We will sample from the same distribution actual and simulated

% %
setofActUnitsDuringPBEs = find(sum(cnctActualTimeBins, 2));
nActUnits = numel(setofActUnitsDuringPBEs);
% %

nTimeBins     = min(totNumTimeBins, totNumTimeBins); % 
rndTBinSubset = randperm(totNumTimeBins, totNumTimeBins);


innerStruct = struct('data', zeros(nTimeBins, 1), 'sim', zeros(nTimeBins, 1));

PFmultiPeakPos.(periods{period}) = innerStruct; 
maxofPFMulti.(periods{period})      = innerStruct;


for ii = 1: nTimeBins
    
    currtb = rndTBinSubset(ii);
    
    currNumActiveUnits = numActiveUnits(currtb);
    
    
    
    if currNumActiveUnits > 0 
        
        
        % active units in the data
        
        currBinActiveUnits = find(cnctActualTimeBins(:, currtb));

        % Find the likeliest postion considering both running directions at the same
        % time
  
        pfMultiplyTwoDirections = prod([spatialTunings_LR2(currBinActiveUnits, :) spatialTunings_RL2(currBinActiveUnits, :)], 1); % concatenated the tunings belonging to two running directions
        
        [maxofPFMulti.(periods{period}).data(ii), peakPosTwodirections] = max(pfMultiplyTwoDirections);
        sharpness.(periods{period}).data(ii) = max(pfMultiplyTwoDirections / sum(pfMultiplyTwoDirections));
        
        peakPos = mod(peakPosTwodirections, nPosBins);
        peakPos(~peakPos) = nPosBins;
        
        PFmultiPeakPos.(periods{period}).data(ii) = peakPos;
        
        
        % active units in a shuffle -- number of units is the
        % same as in data: but the identity of the units were changed
        
%         currBinActiveUnits = setofActUnitsDuringPBEs(randperm(nActUnits, currNumActiveUnits)); 
        currBinActiveUnits = currBinActiveUnits(randperm(currNumActiveUnits));

        pfMultiplyTwoDirections = prod([spatialTunings_LR2(currBinActiveUnits, :) spatialTunings_RL2(currBinActiveUnits, :)], 1);
        
        [maxofPFMulti.(periods{period}).sim(ii), peakPosTwodirections] = max(pfMultiplyTwoDirections); 
        sharpness.(periods{period}).sim(ii) = max(pfMultiplyTwoDirections / sum(pfMultiplyTwoDirections));
        
        peakPos = mod(peakPosTwodirections, nPosBins);
        PFmultiPeakPos.(periods{period}).sim(ii) = peakPos;

    else 
            PFmultiPeakPos.(periods{period}).data(ii) = nan;
            maxofPFMulti.(periods{period}).data(ii)   = nan;
            sharpness.(periods{period}).data(ii)      = nan;
            
            PFmultiPeakPos.(periods{period}).sim(ii)  = nan;
            maxofPFMulti.(periods{period}).sim(ii)    = nan;
            sharpness.(periods{period}).sim(ii)       = nan;
    end
    
end


% plot the histogram of peak positions (in relation to middle of the track) 

subplot(3,3,(period-1)*3+1)

barWidth = 8; % number of position bins 

bins = 0:barWidth: floor(nPosBins/2)+1; 
binCenters = bins + barWidth/2;


% data
currDATA = floor(nPosBins/2) - abs(PFmultiPeakPos.(periods{period}).data-floor(nPosBins/2));
countsDATA = hist(currDATA, bins);
% countsDATA(end) = [];
countsDATA = countsDATA/sum(countsDATA);
 
bar(binCenters, countsDATA, 'facecolor', colors{period}, 'facealpha', 0.4, 'edgecolor', 'none')


% shuffle
currSIM = floor(nPosBins/2) - abs(PFmultiPeakPos.(periods{period}).sim-floor(nPosBins/2));
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

bins = 0:0.1:max([maxofPFMulti.(periods{period}).data maxofPFMulti.(periods{period}).sim]);
binCenters = bins(1:end-1) + 0.05;

% data
currDATA = maxofPFMulti.(periods{period}).data;
countsDATA = hist(currDATA, bins); countsDATA(end)=[]; countsDATA = countsDATA/sum(countsDATA);

bar(binCenters, countsDATA, 'faceColor', colors{period}, 'edgeColor', 'none', 'facealpha', 0.5)

% shuffle
currSIM = maxofPFMulti.(periods{period}).sim;
countsSIM = hist(currSIM, bins); countsSIM(end)=[]; countsSIM = countsSIM/sum(countsSIM);
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
currDATA = sharpness.(periods{period}).data;
countsDATA = hist(currDATA, bins); countsDATA(end)=[]; countsDATA = countsDATA/sum(countsDATA);

bar(binCenters, countsDATA, 'faceColor', colors{period}, 'edgeColor', 'none', 'facealpha', 0.5)

% shuffle
currSIM = sharpness.(periods{period}).sim;
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