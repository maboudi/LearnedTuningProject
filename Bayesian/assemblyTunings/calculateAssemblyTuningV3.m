function [assemblyTunings, px, spikeCountProbs, nCoactiveUnits_wEachUnit, nCoactiveUnits] = calculateAssemblyTuningV3(PBEbinnedFirings, placeFieldsRL, placeFieldsLR, selectPBEs, ifshuffle)

% compared to the previous version "assemblyTuningV2", we calculate the
% posterior within time bins in which the unit fired at least once.
% (when we calculate the spike trigerred average of the posterior, it is worthless to calcualte the posteriors in time bins within which the unit is silent)
% jan 13 2021



binDur = 0.02; % duration of each PBE time bins needed for calculating the posterior probabilities

% place fields of each unit needed for calculating the posteriors 
if isempty(placeFieldsLR)
    
    twoPlaceFieldsFlag = 0;
    
    placeFields = placeFieldsRL;
    placeFields(~placeFields) = 1e-4;
    
    nPosBins = size(placeFields, 2);

else    
    twoPlaceFieldsFlag = 1;
    
    placeFieldsLR(~placeFieldsLR) = 1e-4;
    placeFieldsRL(~placeFieldsRL) = 1e-4;  
    
    nPosBins = size(placeFieldsLR, 2);
    
end


% what subset of PBEs are going to be included in the analysis
if isempty(selectPBEs)
    selectPBEs = 1:size(PBEbinnedFirings(:, 2), 1);
end
PBEbinnedFirings = PBEbinnedFirings(selectPBEs, :);


% concatenate the PBEs 
concatPBEfirings    = cell2mat(PBEbinnedFirings(:, 2)');
[nUnits, nTimeBins] = size(concatPBEfirings);
activeUnits         = find(sum(concatPBEfirings, 2) > 0);
nActiveUnits        = numel(activeUnits);


if ifshuffle == 1
   
   % % generate a unit ID shuffle surrogate

    % divide the units to two (or more) parts based on their average firing rates
    % during the PBEs. Then, shuffle the unit IDs within each group separately.

    concatPBEfirings_ui = zeros(size(concatPBEfirings));

    totalUnitSpikeCount = sum(concatPBEfirings, 2);

    medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));

    highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
    lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);

    concatPBEfirings_ui(highFiringUnits, :) = concatPBEfirings(highFiringUnits(randperm(numel(highFiringUnits))), :);
    concatPBEfirings_ui(lowFiringUnits, :)  = concatPBEfirings(lowFiringUnits(randperm(numel(lowFiringUnits))), :);  
    
    concatPBEfirings = concatPBEfirings_ui;
    
end


maxFiringCount = max(concatPBEfirings(:));


nCoactiveUnits           = sum(concatPBEfirings > 0, 1);
nCoactiveUnits_wEachUnit = cell(nUnits, 1); % number of coactive units with each given unit within each time bin of concatenated PBEs


assemblyTunings = zeros(nUnits, nPosBins);
px              = zeros(nUnits, nPosBins);
spikeCountProbs = zeros(nUnits, nPosBins, maxFiringCount+1);


for iUnit = 1:nActiveUnits
    
    iUnit
    currUnit = activeUnits(iUnit);


    concatPBEfirings_unit       = concatPBEfirings(currUnit, :);  
    concatPBEfirings_OtherUnits = concatPBEfirings(setdiff(1:nUnits, currUnit), :); 

    
    % % added Jan 13 2021
     
%     nonSilentBinsIdx = find(concatPBEfirings_unit > 0);
%     
%     concatPBEfirings_unit = concatPBEfirings_unit(nonSilentBinsIdx);
%     concatPBEfirings_OtherUnits = concatPBEfirings_OtherUnits(:, nonSilentBinsIdx);

    % %
    
    
    
    % the median number of other units coactive with the unit
    nCoactiveUnits_wEachUnit{currUnit} = sum(concatPBEfirings_OtherUnits > 0, 1);


    
    % posterior probabilities given each direction

    if twoPlaceFieldsFlag == 1
        
        postprobRL = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsRL(setdiff(1:nUnits, currUnit), :), binDur);  
        postprobLR = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsLR(setdiff(1:nUnits, currUnit), :), binDur);  


%             %%%%%%% RL direction %%%%%%%%%
% 
%             posteriorProbMatrix = postprobRL;
%             posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
% 
% 
%             [assemblyTunings.RL(currUnit, :),  px.RL(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
%             spikeCountProbs.RL(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
% 
%             %%%%%%% LR direction %%%%%%%%%
% 
%             posteriorProbMatrix = postprobLR;
%             posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
% 
%             [assemblyTunings.LR(currUnit, :), px.LR(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
%             spikeCountProbs.LR(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);


        %%% marginalized over directions %%%

        posteriorProbMatrix = postprobRL + postprobLR;
        posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 

        [assemblyTunings(currUnit, :), px(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit); %#ok<*PFPIE>
        
        
        if nargout > 2
            spikeCountProbs(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
        end

    else

        posteriorProbMatrix = baysDecoder(concatPBEfirings_OtherUnits, placeFields(setdiff(1:nUnits, currUnit), :), binDur); 
        posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]);

        [assemblyTunings(currUnit, :), px(currUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
        
        
        if nargout > 2
            spikeCountProbs(currUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
        end


    end

end


end


function [assemblyTunings, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
%     [assemblyTunings, assemblyTunings_ci, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;
    p_of_x = sum(posteriorProbMatrix, 2)/nTimeBins;
    
    
    assemblyTunings = weightedSummation ./ p_of_x; 
    
    
%     sigma2 = sum(repmat(unitFirings, [nPosBins 1]) .* ((posteriorProbMatrix ./ repmat(p_of_x, [1 nTimeBins])) - assemblyTunings).^2, 2)/ nTimeBins;
%     sigma = sqrt(sigma2);
%     
%     assemblyTunings_ci = assemblyTunings - sigma/sqrt(nTimeBins) * 1.645;
    
    
end


function  spikeCountProbs = calExpectedTunning(posteriorProbMatrix, unitFirings, maxFiringCount)
%     function  [spikeCountProbs, nonZeroSpikeCountProb, silenceProb_otherUnitsSilentPosteriors, silenceProb_otherUnitsFiringPosteriors] = calExpectedTunning(posteriorProbMatrix, unitFirings, nOtherUnitsFiring, maxFiringCount)

% what is the probability of each spikes count (0, 1, 2, etc) given
% the position


%     nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);
    
    firingRange = 0:maxFiringCount;
    
    
    denom = sum(posteriorProbMatrix, 2);
    spikeCountProbs = zeros(1, nPosBins, numel(firingRange));
    
    for iUnit = firingRange
        
        spikeCountProbs(1, :, iUnit+1) = sum(posteriorProbMatrix(:, unitFirings == iUnit), 2)./ denom;
     
    end
    
    
%     % and also the probability of the unit firing at least one given the position
%     nonZeroSpikeCountProb = sum(posteriorProbMatrix(:, unitFirings > 0), 2)./ denom;
%     
%     
%     % the probability of the unit being silent and others being silent too
%     silenceProb_otherUnitsSilentPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring == 0), 2) ./ denom;
%     
%     
%     % the probability of the unit being silent but others firing 
%     silenceProb_otherUnitsFiringPosteriors = sum(posteriorProbMatrix(:, unitFirings == 0 & nOtherUnitsFiring > 0), 2) ./ denom;

    
end
