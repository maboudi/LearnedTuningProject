function posteriorExclude = calculateAssemblyTuning_theta(runBinnedFiring, spikes, clusterQuality, binDur, ifshuffle)


% explanation goes here, the difference with the previous versions



nUnits = numel(spikes);
unitIDs = zeros(nUnits, 2);
for iUnit = 1:nUnits
    unitIDs(iUnit, :) = spikes(iUnit).id;
end



% place fields of each unit needed for calculating the posteriors 

% if isfield(spikes(1).CVspatialTuning, 'RL')
if  isfield(spikes(1).spatialTuning_smoothed, 'RL')
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
% Lratio = [];
% for ii = 1:numel(clusterQuality)
%     
%     temp = clusterQuality(ii).Lratio .*(1-eye(size(clusterQuality(ii).Lratio, 1))); % excluding the diagonal
%     Lratio = [Lratio; temp(:)];
% 
% end
% Lratio_thresh = median(Lratio);

Lratio_thresh = 1e-3;

nRunBouts = numel(runBinnedFiring);


% concatenate the PBEs 


concatBinnedFirings = cell2mat(runBinnedFiring');

activeUnits  = find(sum(concatBinnedFirings, 2) > 0); % possibly there are units that didn't fire in this subset of PBEs
nActiveUnits = numel(activeUnits);


if ifshuffle == 1 
   
   % % generate a unit ID shuffle surrogate

    % divide the set of all units into two (or more) parts based on their average firing rates
    % during the PBEs. Then, shuffle the unit IDs within each group separately.

    totalUnitSpikeCount = sum(concatBinnedFirings, 2);

    medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));

    highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
    lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);
    
    
    runBinnedFiring_ui = cell(nRunBouts, 1);
    for i_run = 1:nRunBouts
       
        runBinnedFiring_ui{i_run} = zeros(size(runBinnedFiring{i_run}));
        
        runBinnedFiring_ui{i_run}(highFiringUnits, :) = runBinnedFiring{i_run}(highFiringUnits(randperm(numel(highFiringUnits))), :);
        runBinnedFiring_ui{i_run}(lowFiringUnits, :)  = runBinnedFiring{i_run}(lowFiringUnits(randperm(numel(lowFiringUnits))), :);  
    
    end
    runBinnedFiring = runBinnedFiring_ui;  
    
end

clear concatBinnedFirings


posteriorExclude = cell(nRunBouts, 1);
for i_run = 1:nRunBouts
    
    nTimeBins = size(runBinnedFiring{i_run}, 2);
    posteriorExclude{i_run} = nan(nPosBins, nTimeBins, nUnits);
end


for iUnit = 1:nActiveUnits
    
    currUnit = activeUnits(iUnit);
    
    currUnitClu   = unitIDs(currUnit, 2);
    currUnitShank = unitIDs(currUnit, 1);
    
    
    % all units from the other shanks are included
    includedUnits_otherShanks = find(unitIDs(:, 1) ~= currUnitShank); 
    
    
    

    % among units on the same shank as of the current unit, units with L-ratio less than threshold L-ratio are included 
    
    if ~isempty(clusterQuality) 

        idx = clusterQuality(currUnitShank).clus == currUnitClu; % row index corresponding to the current unit
        isLessThanMaxLratio = clusterQuality(currUnitShank).Lratio(:, idx) < Lratio_thresh; % idx of the columns with L-ratio below threshold
        acceptedClusters = clusterQuality(currUnitShank).clus(isLessThanMaxLratio); % the clusters corresponding to the qualified columns
        
        includedUnits_sameShanks = find(unitIDs(:,1) == currUnitShank & ismember(unitIDs(:,2), acceptedClusters));
        
            
        otherUnits = [includedUnits_otherShanks; includedUnits_sameShanks];
    else

        otherUnits = setdiff(activeUnits, currUnit);
        otherUnits = sort(otherUnits, 'ascend');
    
    end

    
    % posterior probabilities given each direction
    
    for i_run = 1:nRunBouts
        
        binnedFirings_otherUnits = runBinnedFiring{i_run}(otherUnits, :);
        
        if twoPlaceFieldsFlag == 1

            postprobRL = baysDecoder(binnedFirings_otherUnits, placeFieldsRL(otherUnits, :), binDur);  
            postprobLR = baysDecoder(binnedFirings_otherUnits, placeFieldsLR(otherUnits, :), binDur);  

            posteriorProbMatrix = postprobRL + postprobLR;
            posteriorExclude{i_run}(:, :,  currUnit) = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 

        else

            posteriorProbMatrix = baysDecoder(binnedFirings_otherUnits, placeFields(otherUnits, :), binDur); 
            posteriorExclude{i_run}(:, :,  currUnit) = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]);

        end
    end

end


end
