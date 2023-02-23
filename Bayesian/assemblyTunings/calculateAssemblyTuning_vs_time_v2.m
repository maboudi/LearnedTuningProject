function assemblyTunings = calculateAssemblyTuning_vs_time_v2(PBEInfo, spikes, selectPBEs, clusterQuality, ifshuffle)


% The difference with the previous versions
% posteriors_allPBEs: is the mean posterior calculated for each unit (the
% unit itself was excluded), using all PBEs in the same epoch

binDur = 0.02; % duration of each PBE time bins needed for calculating the posterior probabilities


% make sure the units are matched between spikes and PBEInfo.fr_20msbin
if numel(spikes) ~= size(PBEInfo(1).fr_20msbin)
    error('There is a misatch in number of units between spikes and PBEInfo.fr_20msbin')
end


nUnits = numel(spikes);
unitIDs = zeros(nUnits, 2);
for iUnit = 1:nUnits
    unitIDs(iUnit, :) = spikes(iUnit).id;
end


% place fields of each unit needed for calculating the posteriors 

if isfield(spikes(1).spatialTuning_smoothed, 'RL')
    
    twoPlaceFieldsFlag = 1;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.RL);
    
    placeFieldsLR = zeros(nUnits, nPosBins);
    placeFieldsRL = zeros(nUnits, nPosBins);    
    for iUnit = 1:nUnits
        placeFieldsRL(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.RL;
        placeFieldsLR(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.LR;
    end
    
    placeFieldsRL(~placeFieldsRL) = 1e-4;
    placeFieldsLR(~placeFieldsLR) = 1e-4;
    
else
    twoPlaceFieldsFlag = 0;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    placeFields = zeros(nUnits, nPosBins);
    for iUnit = 1:nUnits
        placeFields(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni;
    end
    
    placeFields(~placeFields) = 1e-4;
end
    


% calculate Lratio threshold by pooling across all shanks
% Lratio = [];
% for ii = 1:numel(clusterQuality)
%     
%     temp = clusterQuality(ii).Lratio .*(1-eye(size(clusterQuality(ii).Lratio, 1))); % excluding the diagonal
%     Lratio = [Lratio; temp(:)];
% 
% end
% Lratio_thresh = median(Lratio);

Lratio_thresh = 1e-3;


% what subset of PBEs are going to be included in the analysis (used in particular for calculating assembly tunings from PBEs with replay scores in a specific range)
tNum_PBEs = numel(PBEInfo);

if isempty(selectPBEs) 
    selectPBEs = 1:numel(PBEInfo);
end


PBEfirings = cell(tNum_PBEs, 1);
for ipbe = 1:tNum_PBEs
    PBEfirings{ipbe} = PBEInfo(ipbe).fr_20msbin; 
end

clear PBEInfo spikes

% concatenate the PBEs 
concatPBEfirings_all = cell2mat(PBEfirings');

activeUnits  = find(sum(concatPBEfirings_all, 2) > 0); % possibly there are units that didn't fire in this subset of PBEs
nActiveUnits = numel(activeUnits);

clear concatPBEfirings_all


if ifshuffle == 1 
   
   % % generate a unit ID shuffle surrogate

    % divide the set of all units into two (or more) parts based on their average firing rates
    % during the PBEs. Then, shuffle the unit IDs within each group separately.
    
    
    
%     totalUnitSpikeCount = sum(concatPBEfirings_all, 2);
% 
%     medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));
%     
%     highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
%     lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);
    
    
    PBEfirings_ui = cell(tNum_PBEs, 1);
    for ipbe = 1:tNum_PBEs
       
        PBEfirings_ui{ipbe} = zeros(size(PBEfirings{ipbe}));
        
        PBEfirings_ui{ipbe}(activeUnits, :) = PBEfirings{ipbe}(activeUnits(randperm(nActiveUnits)) , :);
        
%         PBEfirings_ui{ipbe}(highFiringUnits, :) = PBEfirings{ipbe}(highFiringUnits(randperm(numel(highFiringUnits))), :);
%         PBEfirings_ui{ipbe}(lowFiringUnits, :)  = PBEfirings{ipbe}(lowFiringUnits(randperm(numel(lowFiringUnits))), :);  
    end
    
    PBEfirings = PBEfirings_ui;   
end



nPBE_subsets    = numel(selectPBEs);
assemblyTunings = nan(nUnits, nPosBins, nPBE_subsets);

for iUnit = 1:nActiveUnits
    
    currUnit = activeUnits(iUnit);
    
    currUnitClu   = unitIDs(currUnit, 2);
    currUnitShank = unitIDs(currUnit, 1);
    
    
    % all units from the other shanks are included
    includedUnits_otherShanks = find(unitIDs(:, 1) ~= currUnitShank); 
    
    
    % among units on the same shank as of the current unit, units with L-ratio less than threshold L-ratio are included 
    
    idx = clusterQuality(currUnitShank).clus == currUnitClu; % row index corresponding to the current unit
    isLessThanMaxLratio = clusterQuality(currUnitShank).Lratio(:, idx) < Lratio_thresh; % idx of the columns with L-ratio below threshold
    acceptedClusters = clusterQuality(currUnitShank).clus(isLessThanMaxLratio); % the clusters corresponding to the qualified columns
    
    includedUnits_sameShanks = find(unitIDs(:,1) == currUnitShank & ismember(unitIDs(:,2), acceptedClusters));
    
    
    
    
    otherUnits = [includedUnits_otherShanks; includedUnits_sameShanks];
    otherUnits = sort(otherUnits, 'ascend');
    
    
    PBEfirings_unit       = cell(tNum_PBEs, 1);
    PBEfirings_otherUnits = cell(tNum_PBEs, 1);
    posteriorProbMatrix   = cell(tNum_PBEs, 1);
    
    for ipbe = 1:tNum_PBEs
        
        PBEfirings_unit{ipbe}       = PBEfirings{ipbe}(currUnit, :);
        PBEfirings_otherUnits{ipbe} = PBEfirings{ipbe}(otherUnits, :);
        
        
        if twoPlaceFieldsFlag == 1
        
            postprobRL = baysDecoder(PBEfirings_otherUnits{ipbe}, placeFieldsRL(otherUnits, :), binDur);  
            postprobLR = baysDecoder(PBEfirings_otherUnits{ipbe}, placeFieldsLR(otherUnits, :), binDur);  

            
            %%% marginalized over directions %%%

            posteriorProbMatrix{ipbe} = single(postprobRL + postprobLR);
            posteriorProbMatrix{ipbe} = posteriorProbMatrix{ipbe}./repmat(sum(posteriorProbMatrix{ipbe}, 1), [nPosBins, 1]);    
        else

            postprob = baysDecoder(PBEfirings_otherUnits{ipbe}, placeFields(otherUnits, :), binDur); 
            postprob = single(postprob);
            posteriorProbMatrix{ipbe} = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);

        end
        
    end
    
    clear PBEfirings_otherUnits
    
    
    posteriorProbMatrix_allPBEs = cell2mat(posteriorProbMatrix');
    posteriorProbMatrix_allPBEs = posteriorProbMatrix_allPBEs(:, ~isnan(sum(posteriorProbMatrix_allPBEs, 1)));
    
    px = mean(posteriorProbMatrix_allPBEs, 2);
    clear posteriorProbMatrix_allPBEs
    
    
    for iSub = 1:nPBE_subsets
        
        currPBEs = selectPBEs{iSub};
        
        if ~isempty(currPBEs)
            
            posteriorProbMatrix_sub = posteriorProbMatrix(currPBEs);
            
            posteriorProbMatrix_sub = cell2mat(posteriorProbMatrix_sub');
            idx = ~isnan(sum(posteriorProbMatrix_sub, 1));
            

            concatPBEfirings_unit_sub = PBEfirings_unit(currPBEs);
            concatPBEfirings_unit_sub = cell2mat(concatPBEfirings_unit_sub');

            assemblyTunings(currUnit, :, iSub) = calAssemblyTunings(posteriorProbMatrix_sub(:, idx), concatPBEfirings_unit_sub(idx), px);

        end
    end

end


end


function assemblyTunings = calAssemblyTunings(posteriorProbMatrix, unitFirings, p_of_x)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;    
    assemblyTunings   = weightedSummation ./ p_of_x; 
    

end