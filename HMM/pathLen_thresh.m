function MaxLen = pathLen_thresh(transmat, thresholds)

% path analysis


overlapThresh = 1;

MaxLen = zeros(1, length(thresholds));

for jj = 1 : length(thresholds)
    
    transProbThresh = thresholds(jj);
    
    [~, pathLens] = pathFinder(transmat, transProbThresh, overlapThresh);

    if any(pathLens)
        MaxLen(jj) = max(pathLens);
    else
        MaxLen(jj) = 0;
    end
    
end

