function [spatialTuning, spatialTuning_ci, spatialTuning_shuffle, spikeCountProbs, px, nonZeroSpikeCountProb, silenceProb_otherUnitsSilentPosteriors, silenceProb_otherUnitsFiringPosteriors, nCoactiveUnits_wEachUnit, nCoactiveUnits] = spatialTuning_PBEsDecodedPositions_V5(PBEbinnedFirings, spatialTuningRL, spatialTuningLR, selectPBEs)


% compared to the previous version, nCoactiveUnits has been changed


binDur = 0.02;

spatialTuningLR(~spatialTuningLR) = 1e-4;
spatialTuningRL(~spatialTuningRL) = 1e-4;


if isempty(selectPBEs)
    concatPBEfirings = cell2mat(PBEbinnedFirings(:, 2)');
else
    concatPBEfirings = cell2mat(PBEbinnedFirings(selectPBEs, 2)');
end

maxFiringCount = max(concatPBEfirings(:));


[nUnits, nTimeBins] = size(concatPBEfirings);
nPosBins  = size(spatialTuningRL, 2);

spatialTuning         = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));
spatialTuning_ci      = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));

spatialTuning_shuffle = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));


px = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));


spikeCountProbs         = struct('RL', zeros(nUnits, nPosBins, maxFiringCount+1), 'LR', zeros(nUnits, nPosBins, maxFiringCount+1), 'integrated', zeros(nUnits, nPosBins, maxFiringCount+1));
nonZeroSpikeCountProb = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));

silenceProb_otherUnitsSilentPosteriors = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));
silenceProb_otherUnitsFiringPosteriors = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));


nCoactiveUnits_wEachUnit = cell(nUnits, 1);


% unit ID shuffle
% divide the units to two (or more) parts based on their firing rates
% during the PBEs. Then, shuffle the unitIDs within each group separately.
% 
% totalUnitSpikeCount = sum(concatPBEfirings, 2);
% 
% medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));
% 
% unitsGroup_1 = find(totalUnitSpikeCount >= medSpikeCount);
% unitsGroup_2 = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);
% 
% concatPBEfirings(unitsGroup_1, :) = concatPBEfirings(unitsGroup_1(randperm(numel(unitsGroup_1))), :);
% concatPBEfirings(unitsGroup_2, :) = concatPBEfirings(unitsGroup_2(randperm(numel(unitsGroup_2))), :);


for ii = 1:nUnits
    ii
    
    
    nCoactiveUnits = sum(concatPBEfirings > 0, 1);
    
    concatPBEfirings_unit       = concatPBEfirings(ii, :); 
    concatPBEfirings_OtherUnits = concatPBEfirings(setdiff(1:nUnits, ii), :); 

    
    % the median number of other units coactive with the unit
    
    nFiringUnitEachBin = sum(concatPBEfirings_OtherUnits > 0, 1);
    
    
    nCoactiveUnits_wEachUnit{ii} = zeros(nTimeBins, 1);
    nCoactiveUnits_wEachUnit{ii}(concatPBEfirings_unit > 0) = nFiringUnitEachBin(concatPBEfirings_unit > 0);
    
    
    % posterior probabilities given each direction
    
    postprobRL = baysDecoder(concatPBEfirings_OtherUnits, spatialTuningRL(setdiff(1:nUnits, ii), :), binDur);  
    postprobLR = baysDecoder(concatPBEfirings_OtherUnits, spatialTuningLR(setdiff(1:nUnits, ii), :), binDur);  
    
    
    
    %%%%%%% RL direction %%%%%%%%%
    
    posteriorProbMatrix = postprobRL;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    
    [spatialTuning.RL(ii, :), spatialTuning_ci.RL(ii, :),  px.RL(ii, :)] = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit);
    spatialTuning_shuffle.RL(ii, :)         = calSpatialTuning(posteriorProbMatrix(:, randperm(nTimeBins)), concatPBEfirings_unit);
    
    [spikeCountProbs.RL(ii,:, :), nonZeroSpikeCountProb.RL(ii, :), silenceProb_otherUnitsSilentPosteriors.RL(ii, :), silenceProb_otherUnitsFiringPosteriors.RL(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, nFiringUnitEachBin, maxFiringCount);
    

    
    %%%%%%% LR direction %%%%%%%%%
    
    posteriorProbMatrix = postprobLR;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    
    [spatialTuning.LR(ii, :), spatialTuning_ci.LR(ii, :), px.LR(ii, :)] = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit);
    spatialTuning_shuffle.LR(ii, :)         = calSpatialTuning(posteriorProbMatrix(:, randperm(nTimeBins)), concatPBEfirings_unit);

    
    [spikeCountProbs.LR(ii,:, :), nonZeroSpikeCountProb.LR(ii, :), silenceProb_otherUnitsSilentPosteriors.LR(ii, :), silenceProb_otherUnitsFiringPosteriors.LR(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, nFiringUnitEachBin, maxFiringCount);

    

    
    %%% marginalized over directions %%%
    
    posteriorProbMatrix = postprobRL + postprobLR;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    [spatialTuning.integrated(ii, :), spatialTuning_ci.integrated(ii, :), px.integrated(ii, :)] = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit);
    spatialTuning_shuffle.integrated(ii, :)                 = calSpatialTuning(posteriorProbMatrix(:, randperm(nTimeBins)), concatPBEfirings_unit);

    
    [spikeCountProbs.integrated(ii,:, :), nonZeroSpikeCountProb.integrated(ii, :), silenceProb_otherUnitsSilentPosteriors.integrated(ii, :), silenceProb_otherUnitsFiringPosteriors.integrated(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, nFiringUnitEachBin, maxFiringCount);

end



end







function [spatialTuning, spatialTuning_ci, p_of_x] = calSpatialTuning(posteriorProbMatrix, unitFirings)

    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;
    p_of_x = sum(posteriorProbMatrix, 2)/nTimeBins;
    
    
    spatialTuning = weightedSummation ./ p_of_x; 
    
    
    sigma2 = sum(repmat(unitFirings, [nPosBins 1]) .* ((posteriorProbMatrix ./ repmat(p_of_x, [1 nTimeBins])) - spatialTuning).^2, 2)/ nTimeBins;
    
    sigma = sqrt(sigma2);
    
    spatialTuning_ci = spatialTuning - sigma/sqrt(nTimeBins) * 1.645;
    
    
end


function  [spikeCountProbs, nonZeroSpikeCountProb, silenceProb_otherUnitsSilentPosteriors, silenceProb_otherUnitsFiringPosteriors] = calExpectedTunning(posteriorProbMatrix, unitFirings, nOtherUnitsFiring, maxFiringCount)
    
% what is the probability of each spikes count (0, 1, 2, etc) given
% the position


%     nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);
    
    firingRange = 0:maxFiringCount;
    
    
    denom = sum(posteriorProbMatrix, 2);
    spikeCountProbs = zeros(1, nPosBins, numel(firingRange));
    
    for ii = firingRange
        
        spikeCountProbs(1, :, ii+1) = sum(posteriorProbMatrix(:, unitFirings == ii), 2)./ denom;
     
    end
    
    
    % and also the probability of the unit firing at least one given the position
    nonZeroSpikeCountProb = sum(posteriorProbMatrix(:, unitFirings > 0), 2)./ denom;
    
    
    % the probability of the unit being silent and others being silent too
    silenceProb_otherUnitsSilentPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring == 0), 2) ./ denom;
    
    
    % the probability of the unit being silent but others firing 
    silenceProb_otherUnitsFiringPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring > 0), 2) ./ denom;

    
end
