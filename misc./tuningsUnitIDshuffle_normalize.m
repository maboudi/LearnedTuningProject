function [tuningsAShuffles, tuningsBShuffles, newOrder] = tuningsUnitIDshuffle_normalize(tuningsA, tuningsB, includedUnits, nShuffles, shuffleType)
% [tuningsAShuffles, tuningsBShuffles] = tuningsUnitIDshuffle_normalize(tuningsA, tuningsB, shuffleInclusionFRthresh, nShuffles, shuffleType)

% 
% Tunings 1 and 2 belong to two different directions of running. Regarding
% the unit ID shuffling process, we need to exclude the units with no firing or low firing rates. However, the firing rates could be different between the directions; 
% a unit silent on traversals in one direction can have a
% non-zero firing rate on the other direction. Since in the calculation
% of posterior probaiblities we marginalize over directions (kind of
% summing over the directions), it makes sense to consider the max firing rates
% over the two directions.


if isempty(tuningsB)
   tuningsB = tuningsA; % this is not going to returned by the function anyway ... 
end

nPosBins = size(tuningsA, 2);    


maxFRatesA = max(tuningsA, [], 2);
maxFRatesB = max(tuningsB, [], 2);

% maxFRates = min([maxFRatesA maxFRatesB], [], 2); 


% includedUnitsA  = find(maxFRatesA > shuffleInclusionFRthresh);
% nIncludedUnitsA = numel(includedUnitsA);
% 
% includedUnitsB  = find(maxFRatesB > shuffleInclusionFRthresh);
% nIncludedUnitsB = numel(includedUnitsB);


% % Paremeterizing the place fields in terms of max firing rates and
% % firing distributions with max rate normalized to one; in the shuffling process, only the firing distributions are permuted and the max rates are going to be the same  


tuningsA_norm = tuningsA./repmat(maxFRatesA, [1 nPosBins]);
tuningsB_norm = tuningsB./repmat(maxFRatesB, [1 nPosBins]);


if strcmp(shuffleType, 'same-peak')
        
    tuningsAShuffles = repmat(tuningsA_norm, [1 1 nShuffles]);
    tuningsBShuffles = repmat(tuningsB_norm, [1 1 nShuffles]);
    
    for ii = 1:nShuffles
        
        newOrder = includedUnits(randperm(numel(includedUnits)));        
%         newOrderA = includedUnitsA(randperm(nIncludedUnitsA, nIncludedUnitsA));
%         newOrderB = includedUnitsB(randperm(nIncludedUnitsB, nIncludedUnitsB));
        
        tuningsAShuffles(includedUnits,:, ii) = tuningsA_norm(newOrder, :); 
        tuningsBShuffles(includedUnits,:, ii) = tuningsB_norm(newOrder, :); 

    end
    
    tuningsAShuffles = tuningsAShuffles.*repmat(maxFRatesA, [1 nPosBins]);
    tuningsBShuffles = tuningsBShuffles.*repmat(maxFRatesB, [1 nPosBins]);
    

elseif strcmp(shuffleType, 'regular')
    
    tuningsAShuffles = repmat(tuningsA, [1 1 nShuffles]);
    tuningsBShuffles = repmat(tuningsB, [1 1 nShuffles]);
    
    for ii = 1:nShuffles
        
        newOrder = includedUnits(randperm(numel(includedUnits)));
%         newOrderA = includedUnitsA(randperm(nIncludedUnitsA, nIncludedUnitsA));
%         newOrderB = includedUnitsB(randperm(nIncludedUnitsB, nIncludedUnitsB));
        
        tuningsAShuffles(includedUnits,:, ii) = tuningsA(newOrder, :); 
        tuningsBShuffles(includedUnits,:, ii) = tuningsB(newOrder, :); 

    end
    
end


end

