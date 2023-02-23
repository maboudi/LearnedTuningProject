function [MaxLen, maxPath] = pathLen_thresh(transmat, thresholds)

% path analysis


overlapThresh = 1;

MaxLen = zeros(1, length(thresholds));
maxPath = cell(1, length(thresholds));

for jj = 1 : length(thresholds)
    
    transProbThresh = thresholds(jj);
    
    [seq, pathLens] = pathFinder(transmat, transProbThresh, overlapThresh);

    if any(pathLens)
        [MaxLen(jj),ind] = max(pathLens);
        maxPath(jj) = seq(ind);
        
    else
        MaxLen(jj) = 0;
    end
    
end

