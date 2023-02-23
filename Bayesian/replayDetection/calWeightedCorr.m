function wCorr = calWeightedCorr(postPr)

[nPosBins, nTimeBins] = size(postPr);


   
% weighted mean

temp = postPr .* repmat((1:nPosBins)', [1, nTimeBins]);
posWeightedMean = sum(temp(:))/sum(postPr(:));

temp = postPr .* repmat((1:nTimeBins), [nPosBins, 1]);
timeWeightedMean = sum(temp(:))/sum(postPr(:));



% weighted covariance 

temp = postPr .* (repmat((1:nPosBins)', [1, nTimeBins])-posWeightedMean) .* (repmat((1:nTimeBins), [nPosBins, 1])-timeWeightedMean);
weightedCov = sum(temp(:))/sum(postPr(:));



% autocovariances

temp = postPr .* (repmat((1:nPosBins)', [1, nTimeBins])-posWeightedMean) .* (repmat((1:nPosBins)', [1, nTimeBins])-posWeightedMean);
posAutoCov = sum(temp(:))/sum(postPr(:));

temp = postPr .* (repmat((1:nTimeBins), [nPosBins, 1])-timeWeightedMean) .* (repmat((1:nTimeBins), [nPosBins, 1])-timeWeightedMean);
timeAutoCov = sum(temp(:))/sum(postPr(:));



% weighted correlation

wCorr = weightedCov / sqrt(posAutoCov * timeAutoCov);


end


