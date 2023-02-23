function [PBEInfo, PBEInfo_null, sigMatrix] = BDreplayDetect_separateDirs(PBEInfo, spikes, fileInfo, binDur, calReplayScore)


% replay score separately calculated for each running direction
% replay order score is calculated based on Davidson&Kloosterman, Neuron,
% 2009

% Statistical significance of replay order score can be done by shuffling
% the time bins among the PBEs. The sets of PBEs used in shuffling could be
% either all PBEs or only those with significant replay scores



temp = spikes.spatialTuning_smoothed;
if ~isfield(temp, 'RL')
    fprintf('There is only one running direction, use BDreplayDetect_mixedDir instead')
    return
end

nUnits = numel(fileInfo.okUnits);



% spatial tunings

nPosBins = numel(spikes(1).spatialTuning_smoothed.RL);
runDirections = {'RL'; 'LR'};

for iDir = 1:2
    
    runDir = runDirections{iDir};
    
    spatialTuning.(runDir) = zeros(nUnits, nPosBins);
    unsmoothedSpatialTuning.(runDir) = zeros(nUnits, nPosBins);

    for iunit = 1:nUnits
        currUnit = fileInfo.okUnits(iunit);

        spatialTuning.(runDir)(iunit, :)           = spikes(currUnit).spatialTuning_smoothed.(runDir);
        unsmoothedSpatialTuning.(runDir)(iunit, :) = spikes(currUnit).spatialTuning_unsmoothed.(runDir);

    end

    spatialTuning.(runDir)(~spatialTuning.(runDir)) = 1e-4; % replace all of zeros in tuning with epsilon or a small value
end




nPBEs     = numel(PBEInfo);
nShuffles = 200; %% number of shuffles in each of the methods 


fprintf('\n Calculating replay scores of actuall PBEs .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
for pbe = 1:nPBEs


    currPBE   = PBEInfo(pbe).fr_20msbin; %%% using just the active units
    nTimebins = size(currPBE, 2);
    
    fitLineDictionary = generateLineDictionary(nPosBins, nTimebins);
    
    
    for iDir = 1:2
        
        runDir = runDirections{iDir};
    
        postprob.(runDir) = baysDecoder(currPBE, spatialTuning.(runDir), binDur); 
        

        posteriorProbMat.(runDir) = postprob.(runDir) ./ repmat(sum(postprob.(runDir), 1), [nPosBins, 1]); %% normalized posterior probabilty matrix for each direction
        posteriorProbMat.(runDir)(:, sum(currPBE, 1) == 0) = 0;


        %%% calculate the best fitting line and the replay score and replay
        %%% type for the current event

        [radonIntegral.(runDir), weightedCorr.(runDir), pbeJumpDist, onlyLineElements.(runDir), begPosition.(runDir), endPosition.(runDir)] = replayevaluation_polar2(posteriorProbMat.(runDir), fitLineDictionary);
        

        if endPosition.(runDir)(2) > begPosition.(runDir)(2) % positive slope
            coveredLen.(runDir) = abs(min(endPosition.(runDir)(2), nPosBins) - max(begPosition.(runDir)(2), 1))/nPosBins;
            vmaxsign.(runDir)   = 1;
        elseif  endPosition.(runDir)(2) < begPosition.(runDir)(2)
            coveredLen.(runDir) = abs(min(begPosition.(runDir)(2), nPosBins) - min(endPosition.(runDir)(2), 1))/nPosBins;
            vmaxsign.(runDir)   = -1;
        end 

        
        jumpDist.(runDir).maxJump     = pbeJumpDist(1); % need to correct the code in replayevaluation function so I won't need the assignings here
        jumpDist.(runDir).normMaxJump = pbeJumpDist(2);
        jumpDist.(runDir).medianJump  = pbeJumpDist(3);
        
    end
    
    
    % doing the same for summation of posterior over the directions
    
    posteriorProbMat.mix = (postprob.RL + postprob.LR)./repmat(sum(postprob.RL + postprob.LR, 1), [nPosBins, 1]);
    [radonIntegral.mix, weightedCorr.mix, ~, ~, begPosition.mix, endPosition.mix] = replayevaluation_polar2(posteriorProbMat.mix, fitLineDictionary);
    
    if endPosition.mix(2) >= begPosition.mix(2) % positive slope
        coveredLen.mix = abs(min(endPosition.mix(2), nPosBins) - max(begPosition.mix(2), 1))/nPosBins;
        vmaxsign.mix   = 1;
    else
        coveredLen.mix = abs(min(begPosition.mix(2), nPosBins) - min(endPosition.mix(2), 1))/nPosBins;
        vmaxsign.mix   = -1;
    end 
    
    
    
    
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
    
    PBEInfo(pbe).replayOrderScore = mean((sumLR - sumRL) ./ (sumLR + sumRL)) * PBEInfo(pbe).vmaxsign.mix; 

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


%% testing statistical significance of the replay scores by comparing against shuffle distributions

if ~calReplayScore
    PBEInfo_null = [];
     sigMatrix = [];
    return
end


for iDir = 1:2
    runDir = runDirections{iDir};
    
    radonIntegral_null.(runDir) = zeros(nPBEs, nShuffles, 4); % third dimension corresponds to different methods of shuffle
    coveredLen_null.(runDir)    = zeros(nPBEs, nShuffles, 4);

    weightedCorr_null.(runDir)  = zeros(nPBEs, nShuffles, 4);


    temp = arrayfun(@(x) struct('maxJump', zeros(nShuffles, 4), 'normMaxJump', zeros(nShuffles, 4), 'medianJump', zeros(nShuffles, 4)), 1:nPBEs, 'UniformOutput',0);
    jumpDist_null.(runDir) = horzcat(temp{:}); 
    clear temp
end




%% 1) within-PBE time swap

fprintf('\n Measuring statistical significances by comparing against within-PBE time swap .. \n')
fprintf('Fraction of PBEs processed: ')
flag = 0;
  

for pbe = 1: nPBEs
    
    currPBE = PBEInfo(pbe).fr_20msbin; 
    nTimeBins = size(currPBE, 2);
    fitLineDictionary = generateLineDictionary(nPosBins, nTimeBins);
    
    
    for iDir = 1:2
        
        runDir = runDirections{iDir};
   
        
        postProb = PBEInfo(pbe).posteriorProbMat.(runDir);


        radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
        coveredLen_null_pbe    = zeros(nShuffles, 1);

        weightedCorr_null_pbe  = zeros(nShuffles, 1);

        jumpDist_null_pbe      = zeros(nShuffles, 3); % each column corresponds to one way of calculating the jump distance


        parfor sn = 1: nShuffles

            posteriorProbMatrix_shuffle = genTimeSwap(postProb); % within-PBE time swap

            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_shuffle, fitLineDictionary);

            if endPosition_null(2) >= begPosition_null(2)
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - min(endPosition_null(2), 1))/nPosBins;
            end

        end


        
        radonIntegral_null.(runDir)(pbe, :, 1) = radonIntegral_null_pbe; % concatenating two directions and permuting the dimensions 
        coveredLen_null.(runDir)(pbe, :, 1)    = coveredLen_null_pbe;

        weightedCorr_null.(runDir)(pbe, :, 1)  = weightedCorr_null_pbe;


        f = fieldnames(jumpDist_null);
        for ii = 1: length(f)
            jumpDist_null.(runDir)(pbe).(f{ii})(:, 1) = jumpDist_null_pbe(:, ii);
        end


        
        % replay scores
        
        precentileScore = length(find(radonIntegral_null.(runDir)(pbe, :, 1) < PBEInfo(pbe).radonIntegral.(runDir)))/nShuffles * 100;
        zscore          = (PBEInfo(pbe).radonIntegral.(runDir) - nanmean(radonIntegral_null.(runDir)(pbe, :, 1)))/std(radonIntegral_null.(runDir)(pbe, :, 1));
        rt_ts.(runDir)  = [precentileScore zscore];
        
        
        precentileScore = length(find(abs(weightedCorr_null.(runDir)(pbe, :, 1)) < abs(PBEInfo(pbe).weightedCorr.(runDir))))/nShuffles * 100;
        zscore          = (abs(PBEInfo(pbe).weightedCorr.(runDir)) - nanmean(abs(weightedCorr_null.(runDir)(pbe, :, 1)))) / std(abs(weightedCorr_null.(runDir)(pbe, :, 1)));
        wc_ts.(runDir)  = [precentileScore zscore];
        
    end
    
    PBEInfo(pbe).rt_ts = rt_ts;
    PBEInfo(pbe).wc_ts = wc_ts;
    
    
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
    
    sumFiring = sum(currPBE, 1);
    idx       = (sumFiring == 0);

    pbeActUnits = find(sum(currPBE, 2));
    [tuningshuffles.RL, tuningshuffles.LR] = tuningsUnitIDshuffle_normalize(spatialTuning.RL, spatialTuning.LR, pbeActUnits, nShuffles, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'
        
    
    for iDir = 1:2
        
        runDir = runDirections{iDir};

        % calculate new spatial tunings with shuffled cell IDs

        radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
        coveredLen_null_pbe    = zeros(nShuffles, 1);

        weightedCorr_null_pbe  = zeros(nShuffles, 1);
        jumpDist_null_pbe      = zeros(nShuffles, 3);

        
        curr_tuning_shuffles = tuningshuffles.(runDir);
        parfor sn = 1:nShuffles % shuffle index 

            postprob = baysDecoder(currPBE, curr_tuning_shuffles(:, :, sn), binDur);

            posteriorProbMatrix_PFshuffle = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);
            posteriorProbMatrix_PFshuffle(:, idx) = 0;

            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);

            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - min(endPosition_null(2), 1))/nPosBins;
            end
        end


        radonIntegral_null.(runDir)(pbe, :, 2)  = radonIntegral_null_pbe; % concatenating two directions and permuting the dimensions 
        coveredLen_null.(runDir)(pbe, :, 2)   = coveredLen_null_pbe;
        
        weightedCorr_null.(runDir)(pbe, :, 2) = weightedCorr_null_pbe;


        f = fieldnames(jumpDist_null);
        for ii = 1: length(f)
            jumpDist_null.(runDir)(pbe).(f{ii})(:, 2) = jumpDist_null_pbe(:, ii);
        end

        
        % replay scores
        
        precentileScore = length(find(radonIntegral_null.(runDir)(pbe, :, 2) < PBEInfo(pbe).radonIntegral.(runDir)))/nShuffles * 100;
        zscore          = (PBEInfo(pbe).radonIntegral.(runDir) - nanmean(radonIntegral_null.(runDir)(pbe, :, 2)))/std(radonIntegral_null.(runDir)(pbe, :, 2));
        rt_ui.(runDir)  = [precentileScore zscore];
        
        
        precentileScore = length(find(abs(weightedCorr_null.(runDir)(pbe, :, 2)) < abs(PBEInfo(pbe).weightedCorr.(runDir))))/nShuffles * 100;
        zscore          = (abs(PBEInfo(pbe).weightedCorr.(runDir)) - nanmean(abs(weightedCorr_null.(runDir)(pbe, :, 2)))) / std(abs(weightedCorr_null.(runDir)(pbe, :, 2)));
        wc_ui.(runDir)  = [precentileScore zscore];
        
    end
    
    PBEInfo(pbe).rt_ts = rt_ui;
    PBEInfo(pbe).wc_ts = wc_ui;
    
    
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

for iDir = 1:2
    runDir = runDirections{iDir};
    
    tuningshuffles.(runDir) = circularPlaceFieldShuffle(unsmoothedSpatialTuning.(runDir), nShuffles);
    tuningRLshuffles.(runDir)(~tuningRLshuffles.(runDir)) = 1e-4;
    
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
    
    
    for iDir = 1:2
        
        runDir = runDirections{iDir};

        radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
        coveredLen_null_pbe    = zeros(nShuffles, 1);

        weightedCorr_null_pbe  = zeros(nShuffles, 1);
        jumpDist_null_pbe      = zeros(nShuffles, 3);

        
        curr_tuning_shuffles  = tuningshuffles.(runDir);
        parfor sn = 1:nShuffles % shuffle index

            postprob = baysDecoder(currPBE, curr_tuning_shuffles(:, :, sn), binDur); 

            posteriorProbMatrix_PFshuffle = postprob./repmat(sum(postprob, 1), [nPosBins, 1]);
            posteriorProbMatrix_PFshuffle(:, idx) = 0;

            [radonIntegral_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle, fitLineDictionary);


            if endPosition_null(2) >= begPosition_null(2) % positive slope
                coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPosBins) - max(begPosition_null(2), 1))/nPosBins;
            else
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - min(endPosition_null(2), 1))/nPosBins;
            end
        end

        
        radonIntegral_null.(runDir)(pbe, :, 3) = radonIntegral_null_pbe; % concatenating two directions and permuting the dimensions 
        coveredLen_null.(runDir)(pbe, :, 3)    = coveredLen_null_pbe;

        weightedCorr_null.(runDir)(pbe, :, 3)  = weightedCorr_null_pbe;


        f = fieldnames(jumpDist_null);
        for ii = 1: length(f)
            jumpDist_null.(runDir)(pbe).(f{ii})(:, 3) = jumpDist_null_pbe(:, ii);
        end
    

        % replay scores
        
        precentileScore = length(find(radonIntegral_null.(runDir)(pbe, :, 3) < PBEInfo(pbe).radonIntegral.(runDir)))/nShuffles * 100;
        zscore          = (PBEInfo(pbe).radonIntegral.(runDir) - nanmean(radonIntegral_null.(runDir)(pbe, :, 3)))/std(radonIntegral_null.(runDir)(pbe, :, 3));
        rt_pf.(runDir)  = [precentileScore zscore];
        
        
        precentileScore = length(find(abs(weightedCorr_null.(runDir)(pbe, :, 3)) < abs(PBEInfo(pbe).weightedCorr.(runDir))))/nShuffles * 100;
        zscore          = (abs(PBEInfo(pbe).weightedCorr.(runDir)) - nanmean(abs(weightedCorr_null.(runDir)(pbe, :, 3)))) / std(abs(weightedCorr_null.(runDir)(pbe, :, 3)));
        wc_pf.(runDir)  = [precentileScore zscore];
        
    end
    
    PBEInfo(pbe).rt_pf = rt_pf;
    PBEInfo(pbe).wc_pf = wc_pf;
    
    
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
    
    
    for iDir = 1:2
        runDir = runDirections{iDir};
        
        radonIntegral_null_pbe = zeros(nShuffles, 1); % for individual PBEs
        coveredLen_null_pbe    = zeros(nShuffles, 1);

        weightedCorr_null_pbe  = zeros(nShuffles, 1);
        jumpDist_null_pbe      = zeros(nShuffles, 3);

        
        postprob = PBEInfo(pbe).posteriorProbMat.(runDir);

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
                coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPosBins) - min(endPosition_null(2), 1))/nPosBins;
            end
        end


        radonIntegral_null.(runDir)(pbe, :, 4) = radonIntegral_null_pbe; 
        coveredLen_null.(runDir)(pbe, :, 4)    = coveredLen_null_pbe;

        weightedCorr_null.(runDir)(pbe, :, 4)  = permute(weightedCorr_null_pbe, [2, 1]);


        f = fieldnames(jumpDist_null);
        for ii = 1: length(f)
            jumpDist_null.(runDir)(pbe).(f{ii})(:, 4) = jumpDist_null_pbe(:, ii);
        end


        % replay scores
        
        precentileScore = length(find(radonIntegral_null.(runDir)(pbe, :, 4) < PBEInfo(pbe).radonIntegral.(runDir)))/nShuffles * 100;
        zscore          = (PBEInfo(pbe).radonIntegral.(runDir) - nanmean(radonIntegral_null.(runDir)(pbe, :, 4)))/std(radonIntegral_null.(runDir)(pbe, :, 4));
        rt_ds.(runDir)  = [precentileScore zscore];
        
        
        precentileScore = length(find(abs(weightedCorr_null.(runDir)(pbe, :, 4)) < abs(PBEInfo(pbe).weightedCorr.(runDir))))/nShuffles * 100;
        zscore          = (abs(PBEInfo(pbe).weightedCorr.(runDir)) - nanmean(abs(weightedCorr_null.(runDir)(pbe, :, 4)))) / std(abs(weightedCorr_null.(runDir)(pbe, :, 4)));
        wc_ds.(runDir)  = [precentileScore zscore];
        
    end
    
    PBEInfo(pbe).rt_ds = rt_ds;
    PBEInfo(pbe).wc_ds = wc_ds;
    
    
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

    
    for iDir = 1:2
        
        runDir = runDirections{iDir};
        
        %%% absolute weighted correlation and jump distance
        
        weightedCorr = zeros(nPBEs, 1);
        for pbe = 1:nPBEs
            weightedCorr(pbe) = PBEInfo(pbe).weightedCorr.(runDir);
        end

        maxJump     = zeros(nPBEs, 1);
        normMaxJump = zeros(nPBEs, 1);
        medianJump  = zeros(nPBEs, 1);

        for pbe = 1:nPBEs    
            maxJump(pbe)     = PBEInfo(pbe).jumpDist.(runDir).maxJump;
            normMaxJump(pbe) = PBEInfo(pbe).jumpDist.(runDir).normMaxJump;
            medianJump(pbe)  = PBEInfo(pbe).jumpDist.(runDir).medianJump;
        end


        maxJump_surr     = zeros(nPBEs, nShuffles);
        normMaxJump_surr = zeros(nPBEs, nShuffles);
        medianJump_Surr  = zeros(nPBEs, nShuffles);

        for pbe = 1: nPBEs 
            maxJump_surr(pbe, :)     = jumpDist_null(pbe).(runDir).maxJump(:, shuffleMethod);
            normMaxJump_surr(pbe, :) = jumpDist_null(pbe).(runDir).normMaxJump(:, shuffleMethod);
            medianJump_Surr(pbe, :)  = jumpDist_null(pbe).(runDir).medianJump(:, shuffleMethod);
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
            twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null.(runDir)(:, sn, shuffleMethod), maxJump_surr(:, sn), wcthresholds, jdthresholds);

            % normalized max jump distance
            twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null.(runDir)(:, sn, shuffleMethod), normMaxJump_surr(:, sn), wcthresholds, jdthresholds);

            % median Jump distance
            twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null.(runDir)(:, sn, shuffleMethod), medianJump_Surr(:, sn), wcthresholds, jdthresholds);

        end


        % SIGNIFICANCE MATRIX

        sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}) = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs), 'rt_slp', []);


        for ii = 1: nwcthreshs
           for jj = 1: njdthreshs

               nSurrGreater = length(find(twoDMatrix_surr1(ii, jj, :) >= twoDMatrix1(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}).wc_maxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

               nSurrGreater = length(find(twoDMatrix_surr2(ii, jj, :) >= twoDMatrix2(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}).wc_normMaxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);

               nSurrGreater = length(find(twoDMatrix_surr3(ii, jj, :) >= twoDMatrix3(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}).wc_medianJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);

           end
        end



        %% significance matrix2: replay score and coveredLen

        
        radonIntegral = zeros(nPBEs, 1);
        coveredLen    = zeros(nPBEs, 1);
        for pbe = 1:nPBEs
            radonIntegral(pbe) = PBEInfo(pbe).radonIntegral.(runDir);
            coveredLen(pbe)    = PBEInfo(pbe).coveredLen.(runDir);
        end
        

        %%% significance matrix

        rsthresholds         = 0:0.1:0.9;
        coveredLenthresholds = 0:0.1:0.9; 

        nrsthreshs         = length(rsthresholds);
        ncoveredLenthreshs = length(coveredLenthresholds);


        twoDMatrix = calRScoreJdistMatrix2(radonIntegral, coveredLen, rsthresholds, coveredLenthresholds);


        % SURROGATES

        twoDMatrix_surr = zeros(nrsthreshs, ncoveredLenthreshs, nShuffles); % max jump distance

        for sn = 1: nShuffles   
            twoDMatrix_surr(:, :, sn) = calRScoreJdistMatrix2(radonIntegral_null.(runDir)(:, sn, shuffleMethod), coveredLen_null.(runDir)(:, sn, shuffleMethod), rsthresholds, coveredLenthresholds);
        end


        sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}).rt_slp = zeros(nrsthreshs, ncoveredLenthreshs);
        for ii = 1: nrsthreshs
           for jj = 1: ncoveredLenthreshs

               nSurrGreater = length(find(twoDMatrix_surr(ii, jj, :) >= twoDMatrix(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
               sigMatrix.(runDir).(shuffleMethodNames{shuffleMethod}).rt_slp(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

           end
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
