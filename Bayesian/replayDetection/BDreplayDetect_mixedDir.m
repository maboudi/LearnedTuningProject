function [PBEInfo, PBEInfo_null, sigMatrix] = BDreplayDetect_mixedDir(PBEInfo, spikes, fileInfo, binDur)


% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 

PBEInfo_null = [];
sigMatrix = [];



temp = spikes.spatialTuning_smoothed;
if isfield(temp, 'RL')
    sessionType = 'two_templates';
else
    sessionType = 'one_template';
end



% spatial Tunings
nUnits = numel(fileInfo.okUnits);

if strcmp(sessionType, 'two_templates')
    
    nPosBins = numel(spikes(1).spatialTuning_smoothed.RL);
    
    spatialTuningRL = zeros(nUnits, nPosBins);
    unsmoothedSpatialTuningRL = zeros(nUnits, nPosBins);

    spatialTuningLR = zeros(nUnits, nPosBins);
    unsmoothedSpatialTuningLR = zeros(nUnits, nPosBins);

    for iunit = 1:nUnits
        currUnit = fileInfo.okUnits(iunit);

        spatialTuningRL(iunit, :)           = spikes(currUnit).spatialTuning_smoothed.RL;
        unsmoothedSpatialTuningRL(iunit, :) = spikes(currUnit).spatialTuning_unsmoothed.RL;

        spatialTuningLR(iunit, :)           = spikes(currUnit).spatialTuning_smoothed.LR;
        unsmoothedSpatialTuningLR(iunit, :) = spikes(currUnit).spatialTuning_unsmoothed.LR;
    end
    
    spatialTuningRL(~spatialTuningRL) = 1e-4; % replace all of zeros in tuning with epsilon or a small value
    spatialTuningLR(~spatialTuningLR) = 1e-4;
    
else
    
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    
    spatialTuning = zeros(nUnits, nPosBins);
    unsmoothedSpatialTuning = zeros(nUnits, nPosBins);

    for iunit = 1:nUnits
        currUnit = fileInfo.okUnits(iunit);

        spatialTuning(iunit, :) = spikes(currUnit).spatialTuning_smoothed.uni;
        unsmoothedSpatialTuning(iunit, :) = spikes(currUnit).spatialTuning_unsmoothed.uni;
    end
    
    spatialTuning(~spatialTuning) = 1e-4;
    
end



nPBEs     = numel(PBEInfo);
nShuffles = 500; %% number of shuffles in each of the methods 


fprintf('\n Calculating replay scores of actuall PBEs .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
for pbe = 1:nPBEs


    currPBE   = PBEInfo(pbe).fr_20msbin; %%% using just the active units
    nTimebins = size(currPBE, 2);
    
    
    %%% Bayesian Sequential Score
    %%% calculate the posterior probabilty matrix without any normalization 
    
    if strcmp(sessionType, 'two_templates')
        postprob.RL = baysDecoder(currPBE, spatialTuningRL, binDur); 
        postprob.LR = baysDecoder(currPBE, spatialTuningLR, binDur);

        %%% Summed probability distribution, every column is normalized to one
        %%% (this is the same method as used in Davidson et al. 2009)

        posteriorProbMat = (postprob.RL + postprob.LR)./repmat(sum(postprob.RL + postprob.LR, 1), [nPosBins, 1]);
        
    else
        
        postprob = baysDecoder(currPBE, spatialTuning, binDur); 
        posteriorProbMat = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);

    end
    
    posteriorProbMat(:, sum(currPBE, 1) == 0) = 0;

    
    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
    

    fitLineDictionary = generateLineDictionary(nPosBins, nTimebins);
    [radonIntegral, weightedCorr, pbeJumpDist, onlyLineElements, begPosition, endPosition] = replayevaluation_polar2(posteriorProbMat, fitLineDictionary);

    if endPosition(2) >= begPosition(2) % positive slope
        coveredLen = abs(min(endPosition(2), nPosBins) - max(begPosition(2), 1))/nPosBins;
        vmaxsign   = 1;
    else
        coveredLen = abs(min(begPosition(2), nPosBins) - max(endPosition(2), 1))/nPosBins;
        vmaxsign   = -1;
    end
    
    
    jumpDist.maxJump     = pbeJumpDist(1);
    jumpDist.normMaxJump = pbeJumpDist(2);
    jumpDist.medianJump  = pbeJumpDist(3);

    
    PBEInfo(pbe).postMat_nonNorm  = postprob;
    PBEInfo(pbe).posteriorProbMat = posteriorProbMat;
    PBEInfo(pbe).radonIntegral    = radonIntegral;
    PBEInfo(pbe).coveredLen       = coveredLen;
    PBEInfo(pbe).weightedCorr     = weightedCorr;
    PBEInfo(pbe).jumpDist         = jumpDist;
    PBEInfo(pbe).onlyLineElements = onlyLineElements;
    PBEInfo(pbe).begPosition      = begPosition;
    PBEInfo(pbe).endPosition      = endPosition;
    PBEInfo(pbe).vmaxsign         = vmaxsign;
    
    
    if mod(pbe, 25) == 0
        if flag == 1
            fprintf(repmat('\b', [1, 4]))
        else
            flag = flag + 1;
        end
        fprintf('%.2f', pbe/nPBEs)
        
    elseif pbe == nPBEs
        fprintf(repmat('\b', [1, 4]))
        fprintf('done!')

    end

end





%% replay order score 

if strcmp(sessionType, 'two_templates')
    % if a two template session calculate the replay order score, otherwise skip 

    runDirections = {'RL'; 'LR'};

    for iDir = 1:2

        runDir = runDirections{iDir};

        all_postProb.(runDir) = cell(nPBEs, 1);

        for pbe = 1:nPBEs
           all_postProb.(runDir){pbe} = PBEInfo(pbe).postMat_nonNorm.(runDir); 
        end

    end


    % replay order calculated for actual PBEs

    for pbe = 1:nPBEs
        postProb  = PBEInfo(pbe).postMat_nonNorm;

        sumLR = sum(postProb.LR, 1);
        sumRL = sum(postProb.RL, 1);

        PBEInfo(pbe).replayOrderScore = mean((sumLR - sumRL) ./ (sumLR + sumRL)) * PBEInfo(pbe).vmaxsign; 
    end


    % replay order calculated for pooled time swap PBEs

    replayOrderScore_shuffle = zeros(nPBEs, nShuffles);

    all_postProb_LR = all_postProb.LR;
    all_postProb_RL = all_postProb.RL;

    parfor sn = 1:nShuffles
        postProb_shuffle_LR = genTimeSwap(all_postProb_LR);
        postProb_shuffle_RL = genTimeSwap(all_postProb_RL);

        for pbe = 1:nPBEs

            sumLR = sum(postProb_shuffle_LR{pbe}, 1);
            sumRL = sum(postProb_shuffle_RL{pbe}, 1);

            replayOrderScore_shuffle(pbe, sn) = mean((sumLR - sumRL) ./ (sumLR + sumRL));
        end

    end

    for pbe = 1:nPBEs  
        PBEInfo(pbe).replayOrderScore_prctile = length(find(abs(replayOrderScore_shuffle(pbe, :)) < abs(PBEInfo(pbe).replayOrderScore)))/nShuffles * 100;
    end
end



%% calculating shuffle distributions 

radonIntegral_null = zeros(nPBEs, nShuffles, 4); % third dimension corresponds to different methods of shuffle
coveredLen_null    = zeros(nPBEs, nShuffles, 4);

weightedCorr_null  = zeros(nPBEs, nShuffles, 4);


temp = arrayfun(@(x) struct('maxJump', zeros(nShuffles, 4), 'normMaxJump', zeros(nShuffles, 4), 'medianJump', zeros(nShuffles, 4)), 1:nPBEs, 'UniformOutput',0);
jumpDist_null = horzcat(temp{:}); 
clear temp




%% 1) within-PBE time swap


fprintf('\n Measuring statistical significances by comparing against within-PBE time swap .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
for pbe = 1: nPBEs
    
    
    currPBE = PBEInfo(pbe).fr_20msbin; 
    nTimeBins = size(currPBE, 2);
    
    fitLineDictionary = generateLineDictionary(nPosBins, nTimeBins);
    
    postProb = PBEInfo(pbe).posteriorProbMat;


    radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe    = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe  = zeros(nShuffles, 1);

    jumpDist_null_pbe      = zeros(nShuffles, 3); % each column corresponds to one way of calculating the jump distance
    

    parfor sn = 1: nShuffles
        
        posteriorProbMatrix_shuffle = genTimeSwap(postProb);

        [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_shuffle, fitLineDictionary);
        
        if endPosition_null(2) >= begPosition_null(2)
            coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
        else
            coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
        end

    end
    
    
    radonIntegral_null(pbe, :, 1) = radonIntegral_null_pbe; % concatenating two directions and permuting the dimensions 
    coveredLen_null(pbe, :, 1)    = coveredLen_null_pbe;

    weightedCorr_null(pbe, :, 1)  = weightedCorr_null_pbe;
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 1) = jumpDist_null_pbe(:, ii);
    end


    % replay scores
    
    PBEInfo(pbe).rt_ts  = length(find(radonIntegral_null(pbe, :, 1) < PBEInfo(pbe).radonIntegral))/nShuffles * 100;
%     PBEInfo(pbe).rt_ts.zscore        = (PBEInfo(pbe).radonIntegral - nanmean(radonIntegral_null(pbe, :, 1)))/std(radonIntegral_null(pbe, :, 1));

    PBEInfo(pbe).wc_ts  = length(find(abs(weightedCorr_null(pbe, :, 1)) < abs(PBEInfo(pbe).weightedCorr)))/nShuffles * 100;
%     PBEInfo(pbe).wc_ts.zscore        = (abs(PBEInfo(pbe).weightedCorr) - nanmean(abs(weightedCorr_null(pbe, :, 1)))) / std(abs(weightedCorr_null(pbe, :, 1)));
    
    
    if mod(pbe, 25) == 0
        if flag == 1
            fprintf(repmat('\b', [1, 4]))
        else
            flag = flag + 1;
        end
        fprintf('%.2f', pbe/nPBEs)
        
    elseif pbe == nPBEs
        fprintf(repmat('\b', [1, 4]))
        fprintf('done!')

    end
    
end



%% 2) unit ID shuffle
% There is much redundancy between the procedures for the two shuffle
% methods, can I make it one?


fprintf('\n Measuring statistical significances by comparing against unit ID shuffle .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
for pbe = 1:nPBEs
    
    currPBE = PBEInfo(pbe).fr_20msbin; 
    
    nTimebins = size(currPBE, 2);
    fitLineDictionary = generateLineDictionary(nPosBins, nTimebins);

    
    pbeActUnits = find(sum(currPBE, 2));
    
    
    % calculate new spatial tunings with shuffled cell IDs
    
    
    if strcmp(sessionType, 'two_templates')
        [tuningRLshuffles, tuningLRshuffles] = tuningsUnitIDshuffle_normalize(spatialTuningRL, spatialTuningLR, pbeActUnits, nShuffles, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'
    else
        tuningshuffles = tuningsUnitIDshuffle_normalize(spatialTuning, [], pbeActUnits, nShuffles, 'same-peak'); 
    end

    sumFiring = sum(currPBE, 1);
    idx       = (sumFiring == 0);

    
    radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe    = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe  = zeros(nShuffles, 1);
    jumpDist_null_pbe      = zeros(nShuffles, 3);
    

    if strcmp(sessionType, 'two_templates') 
        
        parfor sn = 1:nShuffles % shuffle index 

            postprobRL = baysDecoder(currPBE, tuningRLshuffles(:, :, sn), binDur); 
            postprobLR = baysDecoder(currPBE, tuningLRshuffles(:, :, sn), binDur);

            posteriorProbMatrix_PFshuffle = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPosBins, 1]);
            posteriorProbMatrix_PFshuffle(:, idx) = 0;
            
            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);
            
            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
            end
        end
            
    else
        
        parfor sn = 1:nShuffles
            
            postprob = baysDecoder(currPBE, tuningshuffles(:, :, sn), binDur);
            posteriorProbMatrix_PFshuffle = postprob./repmat(sum(postprob, 1), [nPosBins, 1]); 
            posteriorProbMatrix_PFshuffle(:, idx) = 0;
            
            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);
            
            
            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
            end
        end
        
    end
       

    radonIntegral_null(pbe, :, 2)  = permute(radonIntegral_null_pbe, [2, 1]); % concatenating two directions and permuting the dimensions 
    coveredLen_null(pbe, :, 2)   = permute(coveredLen_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :, 2) = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 2) = jumpDist_null_pbe(:, ii);
    end


    
    % sequence scores
    
    PBEInfo(pbe).rt_ui  = length(find(radonIntegral_null(pbe, :, 2) < PBEInfo(pbe).radonIntegral))/nShuffles*100;
%     PBEInfo(pbe).rt_ui.zscore        = (PBEInfo(pbe).radonIntegral - nanmean(radonIntegral_null(pbe, :, 2)))/std(radonIntegral_null(pbe, :, 2));

    PBEInfo(pbe).wc_ui  = length(find(abs(weightedCorr_null(pbe, :, 2)) < abs(PBEInfo(pbe).weightedCorr)))/nShuffles * 100;
%     PBEInfo(pbe).wc_ui.zscore        = (abs(PBEInfo(pbe).weightedCorr) - nanmean(abs(weightedCorr_null(pbe, :, 2)))) / std(abs(weightedCorr_null(pbe, :, 2)));

    
    if mod(pbe, 25) == 0
        if flag == 1
            fprintf(repmat('\b', [1, 4]))
        else
            flag = flag + 1;
        end
        fprintf('%.2f', pbe/nPBEs)
        
    elseif pbe == nPBEs
        fprintf(repmat('\b', [1, 4]))
        fprintf('done!')

    end
   
   
end



%% 3) Circular place field shuffle

% generate shuffle place fields

if strcmp(sessionType, 'two_templates')
    tuningRLshuffles = circularPlaceFieldShuffle(unsmoothedSpatialTuningRL, nShuffles);
    tuningRLshuffles(~tuningRLshuffles) = 1e-4;

    tuningLRshuffles = circularPlaceFieldShuffle(unsmoothedSpatialTuningLR, nShuffles);
    tuningLRshuffles(~tuningLRshuffles) = 1e-4;

else
    tuningShuffles   = circularPlaceFieldShuffle(unsmoothedSpatialTuning, nShuffles);
    tuningShuffles(~tuningShuffles) = 1e-4;

end



fprintf('\n Measuring statistical significances by comparing against Circular place field shuffle .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;

for pbe = 1:nPBEs
    
    currPBE = PBEInfo(pbe).fr_20msbin; 
    
    nTimebins = size(currPBE, 2);
    fitLineDictionary = generateLineDictionary(nPosBins, nTimebins);
    
    sumFiring = sum(currPBE, 1);
    idx       = (sumFiring == 0);
    
    
    
    radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe    = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe  = zeros(nShuffles, 1);
    jumpDist_null_pbe      = zeros(nShuffles, 3);

    if strcmp(sessionType, 'two_templates')
        
        parfor sn = 1:nShuffles % shuffle index

            postprobRL = baysDecoder(currPBE, tuningRLshuffles(:, :, sn), binDur); 
            postprobLR = baysDecoder(currPBE, tuningLRshuffles(:, :, sn), binDur);

            posteriorProbMatrix_PFshuffle = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPosBins, 1]);
            posteriorProbMatrix_PFshuffle(:, idx) = 0;

            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);


            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
            end
        end
        
    else
        parfor sn = 1:nShuffles

            postprob = baysDecoder(currPBE, tuningShuffles(:, :, sn), binDur); 
            posteriorProbMatrix_PFshuffle = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);
            posteriorProbMatrix_PFshuffle(:, idx) = 0;

            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);

            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
            end
                

        end 
    end
    
    
    radonIntegral_null(pbe, :, 3) = permute(radonIntegral_null_pbe, [2, 1]); % concatenating two directions and permuting the dimensions 
    coveredLen_null(pbe, :, 3)    = permute(coveredLen_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :, 3)  = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 3) = jumpDist_null_pbe(:, ii);
    end


    % sequence scores
    
    PBEInfo(pbe).rt_pf  = length(find(radonIntegral_null(pbe, :, 3) < PBEInfo(pbe).radonIntegral))/nShuffles*100;
%     PBEInfo(pbe).rt_pf.zscore        = (PBEInfo(pbe).radonIntegral - nanmean(radonIntegral_null(pbe, :, 3)))/std(radonIntegral_null(pbe, :, 3));

    PBEInfo(pbe).wc_pf  = length(find(abs(weightedCorr_null(pbe, :, 3)) < abs(PBEInfo(pbe).weightedCorr)))/nShuffles * 100;
%     PBEInfo(pbe).wc_pf.zscore        = (abs(PBEInfo(pbe).weightedCorr) - nanmean(abs(weightedCorr_null(pbe, :, 3)))) / std(abs(weightedCorr_null(pbe, :, 3)));

    
    if mod(pbe, 25) == 0
        if flag == 1
            fprintf(repmat('\b', [1, 4]))
        else
            flag = flag + 1;
        end
        fprintf('%.2f', pbe/nPBEs)
        
    elseif pbe == nPBEs
        fprintf(repmat('\b', [1, 4]))
        fprintf('done!')
    end
   
end


%% 4) Circular posterior distribution shuffle

fprintf('\n Measuring statistical significances by comparing against circular posterior distribution shuffle .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
for pbe = 1:nPBEs
    
    currPBE = PBEInfo(pbe).fr_20msbin; 
    nTimeBins = size(currPBE, 2);
    
    fitLineDictionary = generateLineDictionary(nPosBins, nTimeBins);
    
    
    radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe    = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe  = zeros(nShuffles, 1);
    jumpDist_null_pbe      = zeros(nShuffles, 3);

    
    postprob = PBEInfo(pbe).posteriorProbMat;
    
    parfor sn = 1:nShuffles % shuffle index
        
        shiftRange = 1:(nPosBins-1);
        
        posteriorProbMatrix_PFshuffle = zeros(size(postprob));
        for ibin = 1:nTimeBins
            posteriorProbMatrix_PFshuffle(:, ibin) = circshift(postprob(:, ibin), shiftRange(randi(numel(shiftRange))), 1);
        end

        
        [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);
        
        if endPosition_null(2) >= begPosition_null(2) % positive slope
            coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
        else
            coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - max(endPosition_null(2), 1))/nPosBins;
        end
    end
    
    
    radonIntegral_null(pbe, :, 4)      = permute(radonIntegral_null_pbe, [2, 1]); 
    coveredLen_null(pbe, :, 4)   = permute(coveredLen_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :, 4) = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 4) = jumpDist_null_pbe(:, ii);
    end


    % sequence scores
    
    PBEInfo(pbe).rt_ds  = length(find(radonIntegral_null(pbe, :, 4) < PBEInfo(pbe).radonIntegral))/nShuffles*100;
%     PBEInfo(pbe).rt_ds.zscore        = (PBEInfo(pbe).radonIntegral - nanmean(radonIntegral_null(pbe, :, 4)))/std(radonIntegral_null(pbe, :, 4));

    PBEInfo(pbe).wc_ds  = length(find(abs(weightedCorr_null(pbe, :, 4)) < abs(PBEInfo(pbe).weightedCorr)))/nShuffles*100;
%     PBEInfo(pbe).wc_ds.zscore        = (abs(PBEInfo(pbe).weightedCorr) - nanmean(abs(weightedCorr_null(pbe, :, 4)))) / std(abs(weightedCorr_null(pbe, :, 4)));
    
    
    if mod(pbe, 25) == 0
        if flag == 1
            fprintf(repmat('\b', [1, 4]))
        else
            flag = flag + 1;
        end
        fprintf('%.2f', pbe/nPBEs)
        
    elseif pbe == nPBEs
        fprintf(repmat('\b', [1, 4]))
        fprintf('done!')

    end
       
end



%% store the null distribution 

PBEInfo_null.radonIntegral = radonIntegral_null;
PBEInfo_null.coveredLen    = coveredLen_null;
PBEInfo_null.weightedCorr  = weightedCorr_null;
PBEInfo_null.jumpDist      = jumpDist_null;



%% two-feature significance matrices (WC vs JD   and    RT vs SLP)

shuffleMethodNames = {'wPBEtimeswap';'unitIDshuffle'; 'circularPFshuffle'; 'circularPosteriorShuffle'};
sigMatrix          = struct('wPBEtimeswap', [], 'unitIDshuffle', [], 'circularPFshuffle', [], 'circularPosteriorShuffle', []);

for shuffleMethod = 1:4


    %%% absolute weighted correlation and jump distance
    
    weightedCorr = [PBEInfo.weightedCorr];

    maxJump     = zeros(nPBEs, 1);
    normMaxJump = zeros(nPBEs, 1);
    medianJump  = zeros(nPBEs, 1);

    for pbe = 1:nPBEs    
        maxJump(pbe)     = PBEInfo(pbe).jumpDist.maxJump;
        normMaxJump(pbe) = PBEInfo(pbe).jumpDist.normMaxJump;
        medianJump(pbe)  = PBEInfo(pbe).jumpDist.medianJump;
    end


    maxJump_surr     = zeros(nPBEs, nShuffles);
    normMaxJump_surr = zeros(nPBEs, nShuffles);
    medianJump_Surr  = zeros(nPBEs, nShuffles);

    for pbe = 1: nPBEs 
        maxJump_surr(pbe, :)     = jumpDist_null(pbe).maxJump(:, shuffleMethod);
        normMaxJump_surr(pbe, :) = jumpDist_null(pbe).normMaxJump(:, shuffleMethod);
        medianJump_Surr(pbe, :)  = jumpDist_null(pbe).medianJump(:, shuffleMethod);
    end

    
    % calculating the matrices of PBE counts passing pairs of threshold on the two
    % metrics

    wcthresholds = 0:0.1:0.9;
    jdthresholds = 0.1:0.1:1;

    nwcthreshs   = length(wcthresholds);
    njdthreshs   = length(jdthresholds);



    % DATA

    % using max jump distance
    twoDMatrix1 = calRScoreJdistMatrix(weightedCorr, maxJump, wcthresholds, jdthresholds);

    % using normalized max jump distance
    twoDMatrix2 = calRScoreJdistMatrix(weightedCorr, normMaxJump, wcthresholds, jdthresholds);

    % using median jump distance
    twoDMatrix3 = calRScoreJdistMatrix(weightedCorr, medianJump, wcthresholds, jdthresholds);



    % SURROGATES

    twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
    twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
    twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

    for sn = 1: nShuffles

        % max jump distance
        twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), maxJump_surr(:, sn), wcthresholds, jdthresholds);

        % normalized max jump distance
        twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), normMaxJump_surr(:, sn), wcthresholds, jdthresholds);

        % median Jump distance
        twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), medianJump_Surr(:, sn), wcthresholds, jdthresholds);

    end


    % SIGNIFICANCE MATRIX

    sigMatrix.(shuffleMethodNames{shuffleMethod}) = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs), 'rt_slp', []);


    for ii = 1: nwcthreshs
       for jj = 1: njdthreshs

           nSurrGreater = length(find(twoDMatrix_surr1(ii, jj, :) >= twoDMatrix1(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_maxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

           nSurrGreater = length(find(twoDMatrix_surr2(ii, jj, :) >= twoDMatrix2(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_normMaxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);

           nSurrGreater = length(find(twoDMatrix_surr3(ii, jj, :) >= twoDMatrix3(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_medianJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);

       end
    end


    
    %% significance matrix2: replay score and coveredLen

    
    radonIntegral = [PBEInfo.radonIntegral];
    coveredLen    = [PBEInfo.coveredLen];

    
    %%% significance matrix

    rsthresholds         = 0:0.1:0.9;
    coveredLenthresholds = 0:0.1:0.9; 

    nrsthreshs         = length(rsthresholds);
    ncoveredLenthreshs = length(coveredLenthresholds);


    twoDMatrix = calRScoreJdistMatrix2(radonIntegral, coveredLen, rsthresholds, coveredLenthresholds);


    
    % SURROGATES

    twoDMatrix_surr = zeros(nrsthreshs, ncoveredLenthreshs, nShuffles); % max jump distance

    for sn = 1: nShuffles   
        twoDMatrix_surr(:, :, sn) = calRScoreJdistMatrix2(radonIntegral_null(:, sn, shuffleMethod), coveredLen_null(:, sn, shuffleMethod), rsthresholds, coveredLenthresholds);
    end


    sigMatrix.(shuffleMethodNames{shuffleMethod}).rt_slp = zeros(nrsthreshs, ncoveredLenthreshs);
    for ii = 1: nrsthreshs
       for jj = 1: ncoveredLenthreshs

           nSurrGreater = length(find(twoDMatrix_surr(ii, jj, :) >= twoDMatrix(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).rt_slp(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

       end
    end


end



end



function tuningShuffles = circularPlaceFieldShuffle(unsmoothedSpatialTuning, nShuffles)

        
[nUnits, nPosBins] = size(unsmoothedSpatialTuning);

tuningShuffles = zeros(nUnits, nPosBins, nShuffles);


win = gausswindow(3,5); % smoothing with a standard deviation of 6 cm (assuming that each positon bins is 2cm long)

for sn = 1:nShuffles

    for iUnit = 1:nUnits

        shuffledTuning_unSmoothed = circshift(unsmoothedSpatialTuning(iUnit, :), randi(nPosBins), 2);
        tuningShuffles(iUnit, :, sn) = conv(shuffledTuning_unSmoothed, win, 'same');
    end
end

end



function count = calRScoreJdistMatrix(weightedCorr, jumpDistance, wcthresholds, jdthresholds)

nwcthreshs = length(wcthresholds);
njdthreshs = length(jdthresholds);

count = zeros(nwcthreshs, njdthreshs);

for ii = 1:nwcthreshs
    
    wcthresh = wcthresholds(ii);
    
    for jj = 1: njdthreshs
        
        jdthresh = jdthresholds(jj);
        
        count(ii, jj) = length(find(abs(weightedCorr) > wcthresh & jumpDistance < jdthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end

end


function count = calRScoreJdistMatrix2(radonIntegral, coveredLen, rsthresholds, coveredLenthresholds)

nrsthreshs = length(rsthresholds);
ncoveredLenthreshs = length(coveredLenthresholds);

count = zeros(nrsthreshs, ncoveredLenthreshs);

for ii = 1:nrsthreshs
    
    rsthresh = rsthresholds(ii);
    
    for jj = 1: ncoveredLenthreshs
        
        coveredLenthresh = coveredLenthresholds(jj);
        
        count(ii, jj) = length(find(radonIntegral > rsthresh & coveredLen > coveredLenthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end


end
