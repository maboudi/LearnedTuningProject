function [spatialTuning, expectation, expectation_firingBIns] = spatialTuning_PBEsDecodedPositions_V3(PBEbinnedFirings, spatialTuningRL, spatialTuningLR, selectPBEs)



binDur = 0.02;

spatialTuningLR(~spatialTuningLR) = 1e-4;
spatialTuningRL(~spatialTuningRL) = 1e-4;


if isempty(selectPBEs)
    concatPBEfrings = cell2mat(PBEbinnedFirings(:, 2)');
else
    concatPBEfrings = cell2mat(PBEbinnedFirings(selectPBEs, 2)');
end

maxFiringCount = max(concatPBEfrings(:));


nUnits = size(concatPBEfrings, 1);
nPosBins  = size(spatialTuningRL, 2);

spatialTuning = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));

expectation          = struct('RL', zeros(nUnits, nPosBins, maxFiringCount+1), 'LR', zeros(nUnits, nPosBins, maxFiringCount+1), 'integrated', zeros(nUnits, nPosBins, maxFiringCount+1));
expectation_firingBIns = struct('RL', zeros(nUnits, nPosBins), 'LR', zeros(nUnits, nPosBins), 'integrated', zeros(nUnits, nPosBins));

for ii = 1:nUnits
    ii
    
    concatPBEfirings_unit       = concatPBEfrings(ii, :); 
    concatPBEfirings_OtherUnits = concatPBEfrings(setdiff(1:nUnits, ii), :); 
    
    
%     silentBins = find(sum(concatPBEfirings_OtherUnits, 1) == 0);
    silentBins = []; 
    
    
    postprobRL = baysDecoder(concatPBEfirings_OtherUnits, spatialTuningRL(setdiff(1:nUnits, ii), :), binDur);  
    postprobLR = baysDecoder(concatPBEfirings_OtherUnits, spatialTuningLR(setdiff(1:nUnits, ii), :), binDur);  
    
    
    
    %%%%%%% RL direction %%%%%%%%%
    
    posteriorProbMatrix = postprobRL;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    
    spatialTuning.RL(ii, :) = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit, silentBins);
    
    [expectation.RL(ii,:, :), expectation_firingBIns.RL(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
    

    
    %%%%%%% LR direction %%%%%%%%%
    
    posteriorProbMatrix = postprobLR;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    
    spatialTuning.LR(ii, :) = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit, silentBins);
    
    [expectation.LR(ii,:, :), expectation_firingBIns.LR(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);

    

    
    %%% marginalized over directions %%%
    
    posteriorProbMatrix = postprobRL + postprobLR;
    posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    spatialTuning.integrated(ii, :) = calSpatialTuning(posteriorProbMatrix, concatPBEfirings_unit, silentBins);

    [expectation.integrated(ii,:, :), expectation_firingBIns.integrated(ii, :)] = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);


end



end


function spatialTuning = calSpatialTuning(posteriorProbMatrix, unitFirings, silentBins)

    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings(setdiff(1:nTimeBins, silentBins)), [nPosBins 1]) .* posteriorProbMatrix(:, setdiff(1:nTimeBins, silentBins)), 2);
    spatialTuning = weightedSummation ./ sum(posteriorProbMatrix(:, setdiff(1:nTimeBins, silentBins)), 2); 

end


function  [expectation, expectation_firingBIns] = calExpectedTunning(posteriorProbMatrix, unitFirings, maxFiringCount)

%     nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);
    
    
    firingRange = 0:maxFiringCount;
    
    
    
    denom = sum(posteriorProbMatrix, 2);
    expectation = zeros(1, nPosBins, numel(firingRange));
    
    for ii = firingRange
        
        expectation(1, :, ii+1) = sum(posteriorProbMatrix(:, unitFirings == ii), 2)./ denom;
     
    end
    
    expectation_firingBIns = sum(posteriorProbMatrix(:, unitFirings > 0), 2)./ denom;
    
    
    

end
