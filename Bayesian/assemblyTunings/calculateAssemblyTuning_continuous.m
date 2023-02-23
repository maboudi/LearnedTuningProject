function [assemblyTunings, px, winCenters] = calculateAssemblyTuning_continuous(spikes, epoch, winDur, stepSize, clusterQuality, ifshuffle)


% Explanation goes here, the difference with the previous versions
% posteriors_allPBEs: is the mean posterior calculated for each unit (the
% unit itself was excluded), using all PBEs in the same epoch




% make sure the units are matched between spikes and PBEInfo.fr_20msbin
% if numel(spikes) ~= size(PBEInfo(1).fr_20msbin)
%     error('There is a misatch in number of units between spikes and PBEInfo.fr_20msbin')
% end


nUnits = numel(spikes);
unitIDs = zeros(nUnits, 2);
for iUnit = 1:nUnits
    unitIDs(iUnit, :) = spikes(iUnit).id;
end


% place fields of each unit needed for calculating the posteriors 

if isfield(spikes(1).spatialTuning_smoothed, 'RL')
    
    twoPlaceFieldsFlag = 1;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.RL);
    
    placeFieldsLR = zeros(nUnits, nPosBins);
    placeFieldsRL = zeros(nUnits, nPosBins);    
    for iUnit = 1:nUnits
        placeFieldsRL(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.RL;
        placeFieldsLR(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.LR;
    end
    
    placeFieldsRL(~placeFieldsRL) = 1e-4;
    placeFieldsLR(~placeFieldsLR) = 1e-4;
    
else
    twoPlaceFieldsFlag = 0;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    placeFields = zeros(nUnits, nPosBins);
    for iUnit = 1:nUnits
        placeFields(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
    end
    
    placeFields(~placeFields) = 1e-4;
end
    


% calculate Lratio threshold by pooling across all shanks
Lratio = [];
for ii = 1:numel(clusterQuality)
    
    temp = clusterQuality(ii).Lratio .*(1-eye(size(clusterQuality(ii).Lratio, 1))); % excluding the diagonal
    Lratio = [Lratio; temp(:)];

end
Lratio_thresh = median(Lratio);




frate = zeros(nUnits, 1);
for ii = 1:nUnits
    frate(ii) =  numel(spikes(ii).time >= epoch(1) & spikes(ii).time < epoch(2))/(epoch(2) - epoch(1));
end
activeUnits = find(frate > 0);


if ifshuffle == 1
    
    % shuffle the unit ID of the spiketimes within each group

    
    medfrate = median(frate(frate > 0));

    hiFiringUnits  = find(frate >= medfrate);
    loFiringUnits  = find(frate > 0 & frate < medfrate);

        
    
    spikes_ui = spikes;

    hiFiringUnits_randIdx = hiFiringUnits(randperm(numel(hiFiringUnits)));
    loFiringUnits_randIdx = loFiringUnits(randperm(numel(loFiringUnits)));


    for iUnit = 1:numel(hiFiringUnits)
        spikes_ui(hiFiringUnits(iUnit)).time = spikes(hiFiringUnits_randIdx(iUnit)).time;
    end

    for iUnit = 1:numel(loFiringUnits)
       spikes_ui(loFiringUnits(iUnit)).time = spikes(loFiringUnits_randIdx(iUnit)).time;
    end

    spikes = spikes_ui;
    clear spikes_ui
    
end



% boundaries of the windows

totalT = epoch(2) - epoch(1);
nWins  = floor((totalT - winDur)/stepSize) + 1;

binStarts = epoch(1) + (0:nWins-1)*stepSize;
binEnds   = binStarts + winDur;
windows   = [binStarts; binEnds]'; 
winCenters = binStarts + stepSize/2;



% Divide the whole period to 20 ms time bins

fileInfo.okUnits = 1:nUnits;

binDur = 0.02;
binnedFirings = timeBinning_20(epoch, spikes, binDur, 0, fileInfo);

nBins = size(binnedFirings, 2);
binCenters  = epoch(1) + (1:nBins)*binDur - binDur/2;



% process only the non-silent bins

tempIdx = sum(binnedFirings, 1) >= 2;
binnedFirings = binnedFirings(:, tempIdx);
binCenters = binCenters(tempIdx);


% bins_subset   = randperm(nBins, floor(.1*nBins));
% binnedFirings = binnedFirings(:, bins_subset);
% binCenters    = binCenters(bins_subset);


assemblyTunings = zeros(nUnits, nPosBins, nWins);

for iUnit = 1:numel(activeUnits)

    currUnit = activeUnits(iUnit);
    
    currUnitClu   = unitIDs(currUnit, 2);
    currUnitShank = unitIDs(currUnit, 1);
    
    
    % all units from the other shanks are included
    includedUnits_otherShanks = find(unitIDs(:, 1) ~= currUnitShank); 
    
    
    % among units on the same shank as of the current unit, units with L-ratio less than threshold L-ratio are included 
    
    idx = clusterQuality(currUnitShank).clus == currUnitClu; % row index corresponding to the current unit
    isLessThanMaxLratio = clusterQuality(currUnitShank).Lratio(:, idx) < Lratio_thresh; % idx of the columns with L-ratio below threshold
    acceptedClusters = clusterQuality(currUnitShank).clus(isLessThanMaxLratio); % the clusters corresponding to the qualified columns
    
    includedUnits_sameShanks = find(unitIDs(:,1) == currUnitShank & ismember(unitIDs(:,2), acceptedClusters));
    
    
    otherUnits = [includedUnits_otherShanks; includedUnits_sameShanks];
    otherUnits = sort(otherUnits, 'ascend');
    
    
    otherUnitsFirings = binnedFirings(otherUnits, :);
    unitFirings       = binnedFirings(currUnit, :);
    
    
    if twoPlaceFieldsFlag == 1

        postprobRL = baysDecoder(otherUnitsFirings, placeFieldsRL(otherUnits, :), binDur);  
        postprobLR = baysDecoder(otherUnitsFirings, placeFieldsLR(otherUnits, :), binDur);  


        %%% marginalized over directions %%%

        posteriorProbMatrix = single(postprobRL + postprobLR);
        posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]);    
    else

        postprob = baysDecoder(otherUnitsFirings, placeFields(otherUnits, :), binDur); 
        postprob = single(postprob);
        posteriorProbMatrix = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);

    end
    
    
    tempIdx =  ~isnan(sum(posteriorProbMatrix, 1));
    posteriorProbMatrix  = posteriorProbMatrix(:, tempIdx);
    binCenters2 = binCenters(tempIdx);

    px = mean(posteriorProbMatrix, 2);
    
   
    for iwin = 1:nWins
        
        currWindow  = windows(iwin, :);
        binIncluded = binCenters2 >= currWindow(1) & binCenters2 < currWindow(2);
        
        unitFirings_win = unitFirings(binIncluded);
        posteriorProbMatrix_win = posteriorProbMatrix(:, binIncluded);
        tempIdx = ~isnan(sum(posteriorProbMatrix_win, 1));

        assemblyTunings(currUnit, :, iwin) = calAssemblyTunings(posteriorProbMatrix_win(:, tempIdx), unitFirings_win(tempIdx), px);
    
    end
end


end



function assemblyTunings = calAssemblyTunings(posteriorProbMatrix, unitFirings, p_of_x)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;    
    assemblyTunings   = weightedSummation ./ p_of_x; 
    
 
end