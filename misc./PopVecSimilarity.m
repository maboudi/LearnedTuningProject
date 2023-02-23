function [PopVectSim, sortAstates] = PopVecSimilarity(x, y)

PopVectSim = acos(x' * y ./ (repmat(diag(sqrt(x'*x)), [1 size(y,2)]) .*  repmat(diag(sqrt(y'*y))', [size(x,2) 1])))*180/pi;

PopVectSim2 = 180 - PopVectSim;   
PopVectSim(PopVectSim > 90) = PopVectSim2(PopVectSim > 90);

[~, maxBstates] = min(PopVectSim, [], 2);
[~, sortAstates] = sort(maxBstates);

% PopVectSim = PopVectSim(sortAstates, :);

end