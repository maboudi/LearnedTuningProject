function [spatialTuning_all, unitIDs_sorted] = spatialTuning_PBEsDecodedPositions(PBEbinnedFirings, postPr)



concatPBEPosteriors = cell2mat(postPr');
concatPBEfrings     = cell2mat(PBEbinnedFirings(:, 2)');
concatPBEfrings(:, sum(concatPBEfrings, 1) == 0) = [];

nUnits   = size(concatPBEfrings, 1);
nPosBins = size(concatPBEPosteriors, 1);


firingUnits = find(sum(concatPBEfrings, 2) > 0);

nFiringUnits = numel(firingUnits);

spatialTuning = zeros(nFiringUnits, nPosBins);
maxPOS        = zeros(nFiringUnits, 1);

for ii = 1:nFiringUnits
    unit = firingUnits(ii);
    
    summation = sum(repmat(concatPBEfrings(unit, :), [nPosBins 1]) .* concatPBEPosteriors, 2); % summation of posterior probabilities over time bins weighted by firing counts within each time bin
    
    spatialTuning(ii, :) = summation / sum(concatPBEfrings(unit, :)); % calculating weighted average using the summation
    
    maxPOS(ii, :) = find(spatialTuning(ii, :) == max(spatialTuning(ii, :)));
    
end

[~, sortInd] = sort(maxPOS, 'descend');

% spatialTuning_sort = spatialTuning(sortInd, :); % sorting based on the tunings peak positions

spatialTuning_all = zeros(nUnits, nPosBins);
spatialTuning_all(firingUnits, :) = spatialTuning ./ repmat(max(spatialTuning, [], 2), [1, nPosBins]); % normalizing to the peak positions 

unitIDs_sorted     = firingUnits(sortInd);


end