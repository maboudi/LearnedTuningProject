function [AT_prefLocation, AT_placeFieldCorr, AT_KLdist, PBEfirings_unit, concatPBEbinCenters] = calculateAssemblyTuning_woAVG(PBEInfo, spikes, clusterQuality, ifshuffle)


% Explanation goes here, the difference with the previous versions
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
    placeFields   = zeros(nUnits, nPosBins); % spatial tunings without considering running directions
    
    for iUnit = 1:nUnits
        placeFieldsRL(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.RL;
        placeFieldsLR(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.LR;
        placeFields(iUnit, :)   = spikes(iUnit).spatialTuning_smoothed.uni;
    end
    
    placeFieldsRL(~placeFieldsRL) = 1e-4;
    placeFieldsLR(~placeFieldsLR) = 1e-4;
    placeFields(~placeFields)     = 1e-4;
    
    
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
Lratio = [];
for ii = 1:numel(clusterQuality)
    
    temp = clusterQuality(ii).Lratio .*(1-eye(size(clusterQuality(ii).Lratio, 1))); % excluding the diagonal
    Lratio = [Lratio; temp(:)];

end
Lratio_thresh = median(Lratio);




% what subset of PBEs are going to be included in the analysis (used in particular for calculating assembly tunings from PBEs with replay scores in a specific range)
ntPBEs = numel(PBEInfo);


PBEfirings    = cell(ntPBEs, 1);
PBEbinCenters = cell(ntPBEs, 1);

for ipbe = 1:ntPBEs
    
    PBEfirings{ipbe}    = PBEInfo(ipbe).fr_20msbin; 
    PBEbinCenters{ipbe} = PBEInfo(ipbe).startT + (1:size(PBEfirings{ipbe}, 2))*binDur - binDur/2;  
end

concatPBEbinCenters = cell2mat(PBEbinCenters');
ntTimeBins = numel(concatPBEbinCenters);



clear PBEInfo spikes

% concatenate the PBEs 
concatPBEfirings_all = cell2mat(PBEfirings');


if ifshuffle == 1 
   
   % % generate a unit ID shuffle surrogate

   
    % divide the set of all units into two (or more) parts based on their average firing rates
    % during the PBEs. Then, shuffle the unit IDs within each group separately.
    
    totalUnitSpikeCount = sum(concatPBEfirings_all, 2);

    medSpikeCount = median(totalUnitSpikeCount(totalUnitSpikeCount > 0));
    
    highFiringUnits = find(totalUnitSpikeCount >= medSpikeCount);
    lowFiringUnits  = find(totalUnitSpikeCount > 0 & totalUnitSpikeCount < medSpikeCount);
    
    
    PBEfirings_ui = cell(ntPBEs, 1);
    for ipbe = 1:ntPBEs
       
        PBEfirings_ui{ipbe} = zeros(size(PBEfirings{ipbe}));
        
        PBEfirings_ui{ipbe}(highFiringUnits, :) = PBEfirings{ipbe}(highFiringUnits(randperm(numel(highFiringUnits))), :);
        PBEfirings_ui{ipbe}(lowFiringUnits, :)  = PBEfirings{ipbe}(lowFiringUnits(randperm(numel(lowFiringUnits))), :);  
    
    end
    PBEfirings = PBEfirings_ui;  
end


activeUnits  = find(sum(concatPBEfirings_all, 2) > 0); % possibly there are units that didn't fire in this subset of PBEs
nActiveUnits = numel(activeUnits);


clear concatPBEfirings_all


PBEfirings_unit = nan(nUnits, ntTimeBins);

AT_prefLocation   = nan(nUnits, ntTimeBins);
AT_placeFieldCorr = nan(nUnits, ntTimeBins);
AT_KLdist         = nan(nUnits, ntTimeBins);

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
    
    
    PBEfirings_unit_pbe   = cell(ntPBEs, 1);
    PBEfirings_otherUnits = cell(ntPBEs, 1);
    posteriorProbMatrix   = cell(ntPBEs, 1);
    
    for ipbe = 1:ntPBEs
        
        PBEfirings_unit_pbe{ipbe}   = PBEfirings{ipbe}(currUnit, :);
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
    
    PBEfirings_unit(currUnit, :) = cell2mat(PBEfirings_unit_pbe');
    PBEfirings_otherUnits = cell2mat(PBEfirings_otherUnits');
    
    
    posteriorProbMatrix_allPBEs = cell2mat(posteriorProbMatrix');
    
    idx = sum(posteriorProbMatrix_allPBEs, 1) > 0 & PBEfirings_unit(currUnit, :) > 0 & sum(PBEfirings_otherUnits, 1) > 0;
    posteriorProbMatrix_allPBEs = posteriorProbMatrix_allPBEs(:, idx);
    
    
    px = mean(posteriorProbMatrix_allPBEs, 2);
    
    
    % pr(n > 0 | x) = pr(x | n > 0) * pr(n > 0) / p(x)
    indivBinTuning = (posteriorProbMatrix_allPBEs * numel(find(PBEfirings_unit(currUnit, :) > 0))/numel(PBEfirings_unit(currUnit, :) > 0))./repmat(px, [1, size(posteriorProbMatrix_allPBEs, 2)]);
    
    
    [~, AT_prefLocation(currUnit, idx)] = max(indivBinTuning);
    AT_placeFieldCorr(currUnit, idx)    = corr(indivBinTuning, placeFields(currUnit, :)'); 
    
    
    % KL distance 
    
    learnedTuning = indivBinTuning;
    learnedTuning = learnedTuning ./ repmat(sum(learnedTuning, 1), [size(learnedTuning, 1) 1]);
    learnedTuning = learnedTuning + eps;
    
    spatialTuning = placeFields(currUnit, :)';
    spatialTuning = spatialTuning ./ repmat(sum(spatialTuning, 1), [size(spatialTuning, 1) 1]);
    spatialTuning = spatialTuning + eps;
    
    AT_KLdist(currUnit, idx) = sum(learnedTuning .* (log(learnedTuning) - log(repmat(spatialTuning, [1 size(learnedTuning, 2)]))), 1);
    
       

    
%     AT_prefLocation2(currUnit).loc    = AT_prefLocation(currUnit, idx);
%     AT_prefLocation2(currUnit).tvec   = concatPBEbinCenters(idx);
%     
%     AT_placeFieldCorr2(currUnit).corr  = AT_placeFieldCorr(currUnit, idx);
%     AT_placeFieldCorr2(currUnit).tvec = concatPBEbinCenters(idx);
    
%     
%     AT_prefLocation2(currUnit).loc    = repelem(AT_prefLocation(currUnit, idx), 1, PBEfirings_unit(currUnit, idx));
%     AT_prefLocation2(currUnit).tvec   = repelem(concatPBEbinCenters(idx), 1, PBEfirings_unit(currUnit, idx));
%     
%     AT_placeFieldCorr2(currUnit).corr  = repelem(AT_placeFieldCorr(currUnit, idx), 1, PBEfirings_unit(currUnit, idx));
%     AT_placeFieldCorr2(currUnit).tvec = repelem(concatPBEbinCenters(idx), 1, PBEfirings_unit(currUnit, idx));
    
    
end
    
end
