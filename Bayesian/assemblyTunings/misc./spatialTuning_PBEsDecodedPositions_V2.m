function spatialTuning = spatialTuning_PBEsDecodedPositions_V2(PBEbinnedFirings, spatialTuningRL, spatialTuningLR, selectPBEs)


binDur = 0.02;

spatialTuningLR(~spatialTuningLR) = 1e-4;
spatialTuningRL(~spatialTuningRL) = 1e-4;

concatPBEfrings = cell2mat(PBEbinnedFirings(selectPBEs, 2)');


nUnits    = size(concatPBEfrings, 1);
nPosBins  = size(spatialTuningRL, 2);

spatialTuning = zeros(nUnits, nPosBins);

for ii = 1:nUnits
    ii
    currPBEfirings = concatPBEfrings(:, concatPBEfrings(ii, :) > 0); % only the time bins during which unit i fired, because those are the time bins that are used at the final step
    nTimeBins = size(currPBEfirings, 2);


    currPBEfirings_unitExcluded = currPBEfirings(setdiff(1:nUnits, ii), :);

    silentBins = find(sum(currPBEfirings_unitExcluded, 1) == 0); % again we need to exluded silent bins when no other unit is firing

    postprobRL = baysDecoder(currPBEfirings_unitExcluded, spatialTuningRL(setdiff(1:nUnits, ii), :), binDur); 
    postprobLR = baysDecoder(currPBEfirings_unitExcluded, spatialTuningLR(setdiff(1:nUnits, ii), :), binDur); 
    
    posteriorProbMatrix_nn = (postprobRL + postprobLR); % ./repmat(sum(postprobRL + postprobLR, 1), [nPosBins, 1])
    posteriorProbMatrix_nn(:, silentBins) = [];
    
    
    summation = sum(repmat(currPBEfirings(ii, setdiff(1:nTimeBins, silentBins)), [nPosBins 1]) .* posteriorProbMatrix_nn, 2); % summation of posterior probabilities over time bins weighted by firing counts within each time bin
%     spatialTuning(ii, :) = summation ./ sum(posteriorProbMatrix_nn, 2); %
%     in the next version of the code
    
    spatialTuning(ii, :) = summation ./ sum(currPBEfirings(ii, setdiff(1:nTimeBins, silentBins))); % calculating weighted average using the summation

end


end