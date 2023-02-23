function [assemblyTunings, assemblyTunings_ui, spikeCountProbs, px, nCoactiveUnits_wEachUnit, nCoactiveUnits] = calculateAssemblyTuning(PBEbinnedFirings, placeFieldsRL, placeFieldsLR, selectPBEs)

% [assemblyTunings, assemblyTunings_ci, assemblyTunings_ui, spikeCountProbs, px, nonZeroSpikeCountProb, silenceProb_otherUnitsSilentPosteriors, silenceProb_otherUnitsFiringPosteriors, nCoactiveUnits_wEachUnit, nCoactiveUnits] = calculateAssemblyTuning(PBEbinnedFirings, placeFieldsRL, placeFieldsLR, selectPBEs)


binDur = 0.02;

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


if isempty(selectPBEs)
    selectPBEs = 1:size(PBEbinnedFirings(:, 2), 1);
end


PBEbinnedFirings = PBEbinnedFirings(selectPBEs, :);
nPBEs = numel(selectPBEs);


concatPBEfirings = cell2mat(PBEbinnedFirings(:, 2)');
[nUnits, nTimeBins] = size(concatPBEfirings);


maxFiringCount = max(concatPBEfirings(:));




nCoactiveUnits_wEachUnit = cell(nUnits, 1);


% % unit ID shuffle

% divide the units to two (or more) parts based on their firing rates
% during the PBEs. Then, shuffle the unitIDs within each group separately.

% concatPBEfirings_ui = zeros(size(concatPBEfirings));

totalUnitSpikeCount = sum(concatPBEfirings, 2);

medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));

highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);

% concatPBEfirings_ui(highFiringUnits, :) = concatPBEfirings(highFiringUnits(randperm(numel(highFiringUnits))), :);
% concatPBEfirings_ui(lowFiringUnits, :)  = concatPBEfirings(lowFiringUnits(randperm(numel(lowFiringUnits))), :);


concatPBEfirings_ui = [];
for iter = 1
    
    PBEbinnedFirings_ui = cell(nPBEs, 1);
    for pbe = 1:nPBEs

        currPBE = PBEbinnedFirings{pbe, 2};

        PBEbinnedFirings_ui{pbe} = zeros(size(currPBE));
        PBEbinnedFirings_ui{pbe}(highFiringUnits, :) = currPBE(highFiringUnits(randperm(numel(highFiringUnits))), :);
        PBEbinnedFirings_ui{pbe}(lowFiringUnits, :)  = currPBE(lowFiringUnits(randperm(numel(lowFiringUnits))), :);

    end

    concatPBEfirings_ui = [concatPBEfirings_ui cell2mat(PBEbinnedFirings_ui')];
end



for iter = 1:2 % for actual dataset and for unit-ID shuffle surrogate
    
    
    if iter == 1
        currentBinnedFirings = concatPBEfirings;
    else
        currentBinnedFirings = concatPBEfirings_ui;
    end
    
    if iter==2 && nargout == 1
        continue
    end


    for iUnit = 1:nUnits

        iUnit 
        nCoactiveUnits = sum(currentBinnedFirings > 0, 1);

        try
            concatPBEfirings_unit       = currentBinnedFirings(iUnit, :); 
            concatPBEfirings_OtherUnits = currentBinnedFirings(setdiff(1:nUnits, iUnit), :); 
        catch
            iUnit
        end


        % the median number of other units coactive with the unit

        nFiringUnitEachBin = sum(concatPBEfirings_OtherUnits > 0, 1);


        nCoactiveUnits_wEachUnit{iUnit} = nans(nTimeBins, 1);
        nCoactiveUnits_wEachUnit{iUnit}(concatPBEfirings_unit > 0) = nFiringUnitEachBin(concatPBEfirings_unit > 0);

        % posterior probabilities given each direction


        if twoPlaceFieldsFlag == 1

            postprobRL = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsRL(setdiff(1:nUnits, iUnit), :), binDur);  
            postprobLR = baysDecoder(concatPBEfirings_OtherUnits, placeFieldsLR(setdiff(1:nUnits, iUnit), :), binDur);  


    %             %%%%%%% RL direction %%%%%%%%%
    % 
    %             posteriorProbMatrix = postprobRL;
    %             posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    % 
    % 
    %             [currAssemblyTunings.RL(iUnit, :),  curr_px.RL(iUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
    %             currSpikeCountProbs.RL(iUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);
    % 
    %             %%%%%%% LR direction %%%%%%%%%
    % 
    %             posteriorProbMatrix = postprobLR;
    %             posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    % 
    %             [currAssemblyTunings.LR(iUnit, :), curr_px.LR(iUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
    %             currSpikeCountProbs.LR(iUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);

            %%% marginalized over directions %%%

            posteriorProbMatrix = postprobRL + postprobLR;
            posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 

            [currAssemblyTunings.integrated(iUnit, :), curr_px.integrated(iUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
            currSpikeCountProbs.integrated(iUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);

        else

            posteriorProbMatrix = baysDecoder(concatPBEfirings_OtherUnits, placeFields(setdiff(1:nUnits, iUnit), :), binDur); 
            posteriorProbMatrix = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]);

            [currAssemblyTunings.integrated(iUnit, :), curr_px.integrated(iUnit, :)] = calAssemblyTunings(posteriorProbMatrix, concatPBEfirings_unit);
            currSpikeCountProbs.integrated(iUnit,:, :) = calExpectedTunning(posteriorProbMatrix, concatPBEfirings_unit, maxFiringCount);


        end


    %         fprintf('% %.1f\n', iUnit/nUnits)
    end

    if iter == 1
        assemblyTunings = currAssemblyTunings;
        spikeCountProbs = currSpikeCountProbs;
        px = curr_px;
    else
        assemblyTunings_ui = currAssemblyTunings;
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
