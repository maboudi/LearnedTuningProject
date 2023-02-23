function [tuningsAShuffles, tuningsBShuffles] = tuningsUnitIDshuffle(tuningsA, tuningsB, shuffleInclusionFRthresh, nShuffles)

% Tunings 1 and 2 belong to two different diretions of running. Regarding
% the unit ID shuffling process, we need to exclude the units with no firing or low firing rates. However, the firing rates could be different between the directions; 
% a unit which is silent on traversals in one direction can have a
% non-zeros firing rate on the other direction. Since, in the calculation
% of posterior probaiblities we marginalize over directions (kind of
% summing over the directions), it makes sense to consider the max firing rates
% over the two directions. 


if isempty(tuningsB)
   tuningsB = tuningsA; 
end

maxFRates = max([tuningsA tuningsB], [], 2); 

includedUnits = find(maxFRates > shuffleInclusionFRthresh);
nIncludedUnits = numel(includedUnits);

tuningsAShuffles = repmat(tuningsA, [1 1 nShuffles]);
tuningsBShuffles = repmat(tuningsB, [1 1 nShuffles]);

for ii = 1:nShuffles
    

    newOrder = includedUnits(randperm(nIncludedUnits, nIncludedUnits));
    
    tuningsAShuffles(includedUnits,:, ii) = tuningsA(newOrder, :); 
    tuningsBShuffles(includedUnits,:, ii) = tuningsB(newOrder, :); 

end

end