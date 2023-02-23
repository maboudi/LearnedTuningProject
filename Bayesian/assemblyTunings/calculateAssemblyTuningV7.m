function [assemblyTunings, px, spikeCountProbs, nCoactiveUnits_wEachUnit, nCoactiveUnits] = calculateAssemblyTuningV7(PBEInfo, spikes, selectPBEs, clusterQuality, ifshuffle)


% explanation goes here, the difference with the previous versions

binDur = 0.02; % duration of each PBE time bins needed for calculating the posterior probabilities


% make sure the units are matched between spikes and PBEInfo.fr_20msbin
if numel(spikes) ~= size(PBEInfo(1).fr_20msbin)
    error('There is a misatch in number of units between spikes and PBEInfo.fr_20msbin')
end


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
% Lratio = [];
% for ii = 1:numel(clusterQuality)
%     
%     temp = clusterQuality(ii).Lratio .*(1-eye(size(clusterQuality(ii).Lratio, 1))); % excluding the diagonal
%     Lratio = [Lratio; temp(:)];
% 
% end

% Lratio_thresh = median(Lratio);

Lratio_thresh = 1e-3;


% what subset of PBEs are going to be included in the analysis (used in particular for calculating assembly tunings from PBEs with replay scores in a specific range)
if isempty(selectPBEs) 
    selectPBEs = 1:numel(PBEInfo);
end


nPBEs = numel(selectPBEs);
PBEbinnedFirings = cell(nPBEs, 1);
for pbe = 1:numel(selectPBEs)
    PBEbinnedFirings{pbe} = PBEInfo(selectPBEs(pbe)).fr_20msbin; 
end



% concatenate the PBEs 
concatPBEfirings = cell2mat(PBEbinnedFirings');

activeUnits  = find(sum(concatPBEfirings, 2) > 0); % possibly there are units that didn't fire in this subset of PBEs
nActiveUnits = numel(activeUnits);


if ifshuffle == 1 
   
   % % generate a unit ID shuffle surrogate

    % divide the set of all units into two (or more) parts based on their average firing rates
    % during the PBEs. Then, shuffle the unit IDs within each group separately.

    concatPBEfirings_ui = zeros(size(concatPBEfirings));

%     totalUnitSpikeCount = sum(concatPBEfirings, 2);
% 
%     medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));
% 
%     highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
%     lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);
% 
%     concatPBEfirings_ui(highFiringUnits, :) = concatPBEfirings(highFiringUnits(randperm(numel(highFiringUnits))), :);
%     concatPBEfirings_ui(lowFiringUnits, :)  = concatPBEfirings(lowFiringUnits(randperm(numel(lowFiringUnits))), :);  
    
    concatPBEfirings_ui(activeUnits, :) = concatPBEfirings(activeUnits(randperm(nActiveUnits)) , :);
    concatPBEfirings = concatPBEfirings_ui;
    
end

nCoactiveUnits           = sum(concatPBEfirings > 0, 1);
nCoactiveUnits_wEachUnit = cell(nUnits, 1); % number of coactive units with each given unit within each time bin of concatenated PBEs


assemblyTunings = zeros(nUnits, nPosBins);
px              = zeros(nUnits, nPosBins);


maxFiringCount  = max(concatPBEfirings(:)); 
spikeCountProbs = zeros(nUnits, nPosBins, maxFiringCount+1); % 3rd dimension corresponds to [0 1 ... maxFiringCount]


for iUnit = 1:nActiveUnits

    
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
%     otherUnits = setdiff(activeUnits, currUnit);
    otherUnits = sort(otherUnits, 'ascend');

    concatPBEfirings_unit       = concatPBEfirings(currUnit, :);   
    concatPBEfirings_OtherUnits = concatPBEfirings(otherUnits, :); 


    
    % the median number of other units coactive with the unit
    nCoactiveUnits_wEachUnit{currUnit} = sum(concatPBEfirings_OtherUnits > 0, 1);


    
    % posterior probabilities given each direction

    if twoPlaceFieldsFlag == 1
        
        postprobRL = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsRL(otherUnits, :), binDur);  
        postprobLR = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsLR(otherUnits, :), binDur);  


        %%% marginalized over directions %%%

        posteriorProbMatrix = postprobRL + postprobLR;
        posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 

        [assemblyTunings(currUnit, :), px(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit); %#ok<*PFPIE>
        
        
%         if nargout > 2
%             spikeCountProbs(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
%         end

    else

        posteriorProbMatrix = baysDecoder(concatPBEfirings_OtherUnits, placeFields(otherUnits, :), binDur); 
        posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]);

        [assemblyTunings(currUnit, :), px(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
        
        
%         if nargout > 2
%             spikeCountProbs(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
%         end


    end

end


end


function [assemblyTunings, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
%     [assemblyTunings, assemblyTunings_ci, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;
    p_of_x = sum(posteriorProbMatrix, 2)/nTimeBins;
    
    
    assemblyTunings = weightedSummation ./ p_of_x; 
    
    
%     sigma2 = sum(repmat(unitFirings, [nPosBins 1]) .* ((posteriorProbMatrix ./ repmat(p_of_x, [1 nTimeBins])) - assemblyTunings).^2, 2)/ nTimeBins;
%     sigma = sqrt(sigma2);
%     
%     assemblyTunings_ci = assemblyTunings - sigma/sqrt(nTimeBins) * 1.645;
    
    
end


% function  spikeCountProbs = calExpectedTunning(posteriorProbMatrix, unitFirings, maxFiringCount)
% %     function  [spikeCountProbs, nonZeroSpikeCountProb, silenceProb_otherUnitsSilentPosteriors, silenceProb_otherUnitsFiringPosteriors] = calExpectedTunning(posteriorProbMatrix, unitFirings, nOtherUnitsFiring, maxFiringCount)
% 
% % what is the probability of each spikes count (0, 1, 2, etc) given
% % the position
% 
% 
% %     nTimeBins = numel(unitFirings);
%     nPosBins  = size(posteriorProbMatrix, 1);
%     
%     firingRange = 0:maxFiringCount;
%     
%     
%     denom = sum(posteriorProbMatrix, 2);
%     spikeCountProbs = zeros(1, nPosBins, numel(firingRange));
%     
%     for iUnit = firingRange
%         
%         spikeCountProbs(1, :, iUnit+1) = sum(posteriorProbMatrix(:, unitFirings == iUnit), 2)./ denom;
%      
%     end
%     
%     
% %     % and also the probability of the unit firing at least one given the position
% %     nonZeroSpikeCountProb = sum(posteriorProbMatrix(:, unitFirings > 0), 2)./ denom;
% %     
% %     
% %     % the probability of the unit being silent and others being silent too
% %     silenceProb_otherUnitsSilentPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring == 0), 2) ./ denom;
% %     
% %     
% %     % the probability of the unit being silent but others firing 
% %     silenceProb_otherUnitsFiringPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring > 0), 2) ./ denom;
% 
%     
% end
