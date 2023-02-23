% This script in intended to calculate ICs 
% from RUN period and their reactivation strengths during sleep periods: PRE and POST. 
% This analysis is based on two time scales: 250ms to detect reactivations of coarse
% time-scale assemblies and 20ms to detect fine time-scale firing regimes. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing spike-times; dividing each period to time bins with an
% amount of overlap between successive time bins


%%%%% training data (from RUN period)

binDur = struct('coarseTS', 0.25, 'fineTS', 0.02);

timeScaleNames = {'coarseTS'; 'fineTS'};
periodNames = {'PRE', 'RUN', 'POST'};

runTraining   = struct('coarseTS', [], 'fineTS', []); % the data that is used for calcualting RUN ICs. No overlap between successive time bins used here
posbinIdx     = struct('coarseTS', [], 'fineTS', []);
linposcenters = struct('coarseTS', [], 'fineTS', []);


posBinSize = 2; % in cm
runSpeedThresh = 0;
binOverlapRatio = 0;

for its = 1:2  
    currTimeScale = timeScaleNames{its};
    [runTraining.(currTimeScale), posbinIdx.(currTimeScale), linposcenters.(currTimeScale)] = binRunData_1D(spikeStruct, behavior.time(2, :), behavior, speed, runSpeedThresh, binDur.(currTimeScale), binOverlapRatio, posBinSize, fileinfo);
end

%%%%% test data (the whole period is binned. Later we add analyses limited to PBEs as well)

qclus = [1 2 3];

binOverlapRatio = struct('coarseTS', 0.5, 'fineTS', 0);

basicStruct  = struct('PRE', [], 'RUN', [], 'POST', []);
binCenters   = struct('coarseTS', basicStruct, 'fineTS', basicStruct);
binnedFiring = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

for iperiod = 1:3
    for its = 1:2 
        
        currPeriod          = periodNames{iperiod};
        currTimeScale       = timeScaleNames{its};
        currBinDur          = binDur.(currTimeScale);
        currBinOverlabRatio = binOverlapRatio.(currTimeScale);
        currStepSize        = (1-currBinOverlabRatio) * currBinDur;
        
        
        binCenters.(currTimeScale).(currPeriod)    = (behavior.time(iperiod, 1):currStepSize:(behavior.time(iperiod, 2)-currBinDur)) + currBinDur/2;
        binnedFiring.(currTimeScale).(currPeriod)  = timeBinning_withOverlap(behavior.time(iperiod, :) , spikeStruct, qclus, currBinDur, currBinOverlabRatio, fileinfo); % the theta periods
    end
end
  


%%%%% mean and standrad deviation of firing rates for the same z-scoring throughout the analysis 

firingMean = struct('coarseTS', basicStruct, 'fineTS', basicStruct);
firingStd  = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

for its = 1:2
    for iperiod = 1:3
        
        currTimeScale = timeScaleNames{its}; 
        currPeriod    = periodNames{iperiod};
        
        currData      = binnedFiring.(currTimeScale).(currPeriod){1,2};

        firingMean.(currTimeScale).(currPeriod) = mean(currData, 2); % in Hz
        firingStd.(currTimeScale).(currPeriod)  = std(currData, 0, 2); % in Hz
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the RUN ICs and IC weights

ICs         = struct('coarseTS', [], 'fineTS', []);
ICweights   = struct('coarseTS', [], 'fineTS', []);
eigenValues = struct('coarseTS', [], 'fineTS', []);
evThresh    = struct('coarseTS', [], 'fineTS', []); % eigen value threshold


for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    
    meanFiringCount = firingMean.(currTimeScale).RUN;
    stdFiringCount  = firingStd.(currTimeScale).RUN;

    [ICs.(currTimeScale), ICweights.(currTimeScale), eigenValues.(currTimeScale), evThresh.(currTimeScale)] = pca_ica_v2(runTraining.(currTimeScale), meanFiringCount, stdFiringCount);
    nICs = size(ICs.(currTimeScale), 1);

    for ii = 1:nICs
        ICweights.(currTimeScale)(ii, :) = ICweights.(currTimeScale)(ii, :)/norm(ICweights.(currTimeScale)(ii, :));
        
        %signs are set such that the highest absolute weight of each pattern is
        % always positive
        
        [~, maxAbsWeightIndex] = max(abs(ICweights.(currTimeScale)(ii, :)));
        maxWeightSign = sign(ICweights.(currTimeScale)(ii, maxAbsWeightIndex));
        ICweights.(currTimeScale)(ii, :) = ICweights.(currTimeScale)(ii,:) * maxWeightSign;
        
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate ICs in jitered data
% I realized that ICA is not actually sensitive to the temporal shuffles
% (time swap or separately for each unit), because it only cares about the
% distributions of firing rates. So, this section of the code can be
% removed.

% generating jittered data

jitterLength.fineTS   = 20;
jitterLength.coarseTS = 5;

nUnits = size(runTraining.coarseTS, 1);
shuffleBinIdx = struct('fineTS', [], 'coarseTS', []);


for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    currJitterLen = jitterLength.(currTimeScale);
    
    nTimeBins = size(runTraining.(currTimeScale), 2);
    
    binIdx = 1:nTimeBins;
 
    loopVariable = repmat(binIdx, [nUnits 1]);
    
    %%% calculating jittered bin-indices
    
    parfor jj= 1:nUnits
        
        alreadyUsedBins = zeros(nTimeBins, 1);
        currentSign = -1;
        
        loop2Variable = loopVariable(jj , :);
        for ii = 1: floor(nTimeBins/2)

            availableBins = setdiff(binIdx, find(alreadyUsedBins));

            bin1 = availableBins(randi(numel(availableBins)));

            flag = 0;

            if find(ismember(setdiff((bin1-currJitterLen):(bin1+currJitterLen), bin1), availableBins))
                while flag == 0
                    currentSign = currentSign*-1;
                    bin2 = bin1 + randi(currJitterLen)*currentSign;

                    if ismember(bin2, availableBins)
                       flag = 1;
                    end
                end

                loop2Variable(bin1) = bin2;
                loop2Variable(bin2) = bin1;

                alreadyUsedBins([bin1 bin2]) = 1;
                
            else
                alreadyUsedBins(bin1) = 1;
            end

        end
        
        loopVariable(jj, :) =  loop2Variable;
    end
    
    shuffleBinIdx.(currTimeScale) = loopVariable;
    
    
end


% calculating surrogate ICs

ICs_s         = struct('coarseTS', [], 'fineTS', []);
ICweights_s   = struct('coarseTS', [], 'fineTS', []);
eigenValues_s = struct('coarseTS', [], 'fineTS', []);
evThresh_s    = struct('coarseTS', [], 'fineTS', []); % eigen value threshold


for its = 1:2
    
    currTimeScale  = timeScaleNames{its};
    currShuffleIdx = shuffleBinIdx.(currTimeScale);
    
    meanFiringCount = firingMean.(currTimeScale).RUN;
    stdFiringCount  = firingStd.(currTimeScale).RUN;
    
    shuffledTrainData = zeros(size(runTraining.(currTimeScale)));
    for iUnit = 1:nUnits 
        shuffledTrainData(iUnit, :) = runTraining.(currTimeScale)(iUnit, currShuffleIdx(iUnit, :));
    end
    
    [ICs_s.(currTimeScale), ICweights_s.(currTimeScale), eigenValues_s.(currTimeScale), evThresh_s.(currTimeScale)] = pca_ica_v2(shuffledTrainData, meanFiringCount, stdFiringCount);
    nICs = size(ICs_s.(currTimeScale), 1);

    for ii = 1:nICs
        ICweights_s.(currTimeScale)(ii, :) = ICweights_s.(currTimeScale)(ii, :)/norm(ICweights_s.(currTimeScale)(ii, :));
        
        %signs are set such that the highest absolute weight of each pattern is
        % always positive
        
        [~, maxAbsWeightIndex] = max(abs(ICweights_s.(currTimeScale)(ii, :)));
        maxWeightSign = sign(ICweights_s.(currTimeScale)(ii, maxAbsWeightIndex));
        ICweights_s.(currTimeScale)(ii, :) = ICweights_s.(currTimeScale)(ii,:) * maxWeightSign;
        
    end    
       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating the (re)activation strengths

activationStrength = struct('coarseTS', [], 'fineTS', []);
activationStrength = structfun(@(x) struct('PRE', [], 'RUN', [], 'POST', []), activationStrength, 'UniformOutPut', false);

overallActivationStrength = struct('coarseTS', [], 'fineTS', []);
overallActivationStrength = structfun(@(x) struct('PRE', [], 'RUN', [], 'POST', []), overallActivationStrength, 'UniformOutPut', false);


for its = 1:2
    activationStrength.(timeScaleNames{its})        = structfun(@(x) struct('original', [], 'smoothed', []), activationStrength.(timeScaleNames{its}), 'UniformOutPut', false);
    overallActivationStrength.(timeScaleNames{its}) = structfun(@(x) struct('original', [], 'smoothed', []), overallActivationStrength.(timeScaleNames{its}), 'UniformOutPut', false);
end


avgWinDur  = 15*60; % running average window in sec

for its = 1:2  
    
    currTimeScale = timeScaleNames{its};
    
    avgWinSize = avgWinDur/binDur.(currTimeScale);
    avgWindow  = ones(avgWinSize, 1)/avgWinSize;
    
    for iperiod = 1:3      
        currPeriod    = periodNames{iperiod};

        meanFiringCount = firingMean.(currTimeScale).(currPeriod);
        stdFiringCount  = firingStd.(currTimeScale).(currPeriod);
        
        activationStrength.(currTimeScale).(currPeriod).original  = calActStrength(ICweights.(currTimeScale), binnedFiring.(currTimeScale).(currPeriod){1,2}, meanFiringCount, stdFiringCount);
        activationStrength.(currTimeScale).(currPeriod).smoothed  = arrayConv(activationStrength.(currTimeScale).(currPeriod).original, avgWindow); % running average
        
        % overall activation strength by summing over all ICs
        overallActivationStrength.(currTimeScale).(currPeriod).original  = mean(activationStrength.(currTimeScale).(currPeriod).original, 1);
        overallActivationStrength.(currTimeScale).(currPeriod).smoothed  = arrayConv(overallActivationStrength.(currTimeScale).(currPeriod).original, avgWindow); % running average
        
    end
end

%% z-score using PRE-distribution
 
for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    
    PREmean = nanmean(activationStrength.(currTimeScale).PRE.smoothed, 2);
    PREstd  = nanstd(activationStrength.(currTimeScale).PRE.smoothed, 0, 2);
    
    PREmeanOverall = nanmean(overallActivationStrength.(currTimeScale).PRE.smoothed, 2);
    PREstdOverall  = nanstd(overallActivationStrength.(currTimeScale).PRE.smoothed, 0, 2);
    
    for iperiod = 1:3
        currPeriod = periodNames{iperiod};

        activationStrength.(currTimeScale).(currPeriod).smoothed        = zscore_KM(activationStrength.(currTimeScale).(currPeriod).smoothed, PREmean, PREstd);
        overallActivationStrength.(currTimeScale).(currPeriod).smoothed = zscore_KM(overallActivationStrength.(currTimeScale).(currPeriod).smoothed, PREmeanOverall, PREstdOverall);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating reactivation strengths within PBEs

binnedFiringPBEs  = struct('PRE', [], 'RUN', [], 'POST', []);

binnedFiringPBEs.PRE  = PREbinnedPBEs.data;
binnedFiringPBEs.RUN  = RUNbinnedPBEs.data;
binnedFiringPBEs.POST = POSTbinnedPBEs.data;


activationStrPBEs = struct('coarseTS', 0.5, 'fineTS', 0);
activationStrPBEs = structfun(@(x) struct('PRE', [], 'RUN', [], 'POST', []), activationStrPBEs, 'UniformOutPut', false);


PBEperiods = struct('PRE', [], 'RUN', [], 'POST', []);
PBEperiods.PRE  = secondaryPBEs_pre;
PBEperiods.RUN  = secondaryPBEs_run;
PBEperiods.POST = secondaryPBEs_post;


for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    for iperiod = 1:3

        meanFiringCount = firingMean.(currTimeScale).(currPeriod);
        stdFiringCount  = firingStd.(currTimeScale).(currPeriod);

        currPeriod = periodNames{iperiod};
        currPBEs = binnedFiringPBEs.(currPeriod);
        
        nPBEs = size(currPBEs, 1);
        activationStrPBEs.(currTimeScale).(currPeriod) = cell(nPBEs, 1);
        
        for pbe = 1:nPBEs
            
            if its == 1
                activationStrPBEs.(currTimeScale).(currPeriod){pbe} = calActStrength(ICweights.(currTimeScale), sum(currPBEs{pbe, 2}, 2), meanFiringCount, stdFiringCount);
            elseif its == 2
                activationStrPBEs.(currTimeScale).(currPeriod){pbe} = calActStrength(ICweights.(currTimeScale), currPBEs{pbe, 2}, meanFiringCount, stdFiringCount);
            end
                   
        end
    end
end


% mean reactivation strengths by averaging over ICs

activationStrPBEs_sumOverICs = struct('coarseTS', 0.5, 'fineTS', 0); % later we are going to investigate how fine- and coarse-time scale ICs reactivate together within individual PBEs
activationStrPBEs_sumOverICs = structfun(@(x) struct('PRE', [], 'RUN', [], 'POST', []), activationStrPBEs_sumOverICs, 'UniformOutPut', false);


for iperiod = 1:3

    currPeriod = periodNames{iperiod};
    currPBEs = binnedFiringPBEs.(currPeriod);

    nPBEs = size(currPBEs, 1);

    activationStrPBEs_sumOverICs.fineTS.(currPeriod)   = zeros(nPBEs, 1);
    activationStrPBEs_sumOverICs.coarseTS.(currPeriod) = zeros(nPBEs, 1);
    
    for pbe = 1:nPBEs
        
        temp = activationStrPBEs.coarseTS.(currPeriod){pbe};
        activationStrPBEs_sumOverICs.coarseTS.(currPeriod)(pbe) = mean(temp(temp~=0));
        
        temp = activationStrPBEs.fineTS.(currPeriod){pbe}(:);
        activationStrPBEs_sumOverICs.fineTS.(currPeriod)(pbe)   = mean(temp(temp~=0));

    end

end


figure;

for iperiod = 1:3
    
    subplot(1,3,iperiod)
    currPeriod = periodNames{iperiod};
    
    plot(activationStrPBEs_sumOverICs.coarseTS.(currPeriod), activationStrPBEs_sumOverICs.fineTS.(currPeriod), '.', 'markersize', 5)
    
    xlabel('mean coarseTS IC reactivation', 'fontsize', 12)
    ylabel('mean fineTS IC reactivation', 'fontsize', 12)
    title(currPeriod, 'fontsize', 12)
    
    xlim([0 30])
    ylim([0 15])
    
    axis square
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigating the degree to which ICs' reactivations and reconstructed positions are associated together when we look within individual time bins during RUN PBEs   
% only for fine-time scale ICs

nICs = size(ICs.fineTS, 1);

nPosBins = size(spatialTunings_LR, 2);

spatialTunings_LR(spatialTunings_LR==0) = 1e-4;
spatialTunings_RL(spatialTunings_RL==0) = 1e-4;

ICAvgPosterior    = struct('PRE', [], 'RUN', [], 'POST', []);
ICpeakPositions   = struct('PRE', [], 'RUN', [], 'POST', []);
ICpositionsortInd = struct('PRE', [], 'RUN', [], 'POST', []);


for iperiod = 1:3
    
    currPeriod           = periodNames{iperiod};
    cnctPBEs             = cell2mat(binnedFiringPBEs.(currPeriod)(:, 2)');
    cnctPBEActivationStr = cell2mat(activationStrPBEs.fineTS.(currPeriod)'); % concatenated reactivation strengths during the PBEs
    
    postPr = struct('LR', [], 'RL', [], 'integrated', []);
    
    postPr.LR = baysDecoder(cnctPBEs, spatialTunings_LR, binDur.fineTS);
    postPr.RL = baysDecoder(cnctPBEs, spatialTunings_RL, binDur.fineTS);
    
    postPr.integrated = (postPr.LR + postPr.RL) ./ repmat(sum(postPr.LR + postPr.RL), [nPosBins, 1]);
    
    postPr.LR = postPr.LR ./ repmat(sum(postPr.LR, 1), [size(postPr.LR, 1) 1]);
    postPr.RL = postPr.RL ./ repmat(sum(postPr.RL, 1), [size(postPr.RL, 1) 1]);

    
    for ic = 1:nICs
        Averagingwts = cnctPBEActivationStr(ic, :);
        ICAvgPosterior.(currPeriod)(:, ic) = sum(repmat(Averagingwts, [nPosBins 1]) .* postPr.integrated , 2) ./ sum(Averagingwts);
    end

    [~, ICpeakPositions.(currPeriod)] = max(ICAvgPosterior.(currPeriod));
    [~, ICpositionsortInd.(currPeriod)] = sort(ICpeakPositions.(currPeriod), 'ascend');
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigating the degree to which ICs' and actual position of the animal are associated
% for both fine- and coarse-time scale ICs but only during RUN

ICtunings = struct('coarseTS', [], 'fineTS', []);
for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    nICs = size(ICs.(currTimeScale), 1);

    uniqPositions = unique(posbinIdx.(currTimeScale)(:, 1));
     
    
    ICtunings.(currTimeScale) = zeros(nICs, numel(uniqPositions));
    for ii = 1: numel(uniqPositions)
        ICtunings.(currTimeScale)(:, ii) = sum(activationStrength.(currTimeScale).RUN.original(:, posbinIdx.(currTimeScale)(:, 1) == uniqPositions(ii)), 2); 
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross-correlograms of ICs' reactivations and temporal biases
% only for fine-time scale PBEs 

nICs = size(ICs.fineTS, 1);

maxLag    = 20; % we set a maximum lag to be able to sum cross-correlation profiles calculated from all PBEs. Otherwise, we will have cross-correlogram with varying length from different PBEs ... 
nShuffles = 3;


% cross-correlations
xcorr_actual = struct('PRE', [], 'RUN', [], 'POST', []); % ICs' cross-correlation, separately calculated within each PBE and then summed
xcorr_actual = structfun(@(x) cell(nICs, nICs), xcorr_actual, 'UniformOutPut', false);

xcorr_null = struct('PRE', [], 'RUN', [], 'POST', []); % ICs' cross-correlation, separately calculated within each PBE and then summed
xcorr_null = structfun(@(x) cell(nICs, nICs), xcorr_null, 'UniformOutPut', false);


% temporal biases calculated based on the cross-correlations
temporalBias = struct('PRE', [], 'RUN', [], 'POST', []);
temporalBias = structfun(@(x) nan(nICs, nICs), temporalBias, 'UniformOutPut', false);

temporalBias_null = struct('PRE', [], 'RUN', [], 'POST', []);
temporalBias_null = structfun(@(x) cell(nICs, nICs), temporalBias_null, 'UniformOutPut', false);


for iperiod = 1:3
    
    iperiod
    
    currPeriod = periodNames{iperiod};
    nPBEs      = size(binnedFiringPBEs.(currPeriod), 1);
    
    for ii = 1:nICs
        ii
        for jj = ii:nICs
            
            %%%%%%%%%%%%%
            %%% ACTUAL %%
            
            % cross-correlations
            xcorr_actual_PBE = zeros(nPBEs, 2*maxLag+1);
            for pbe = 1:nPBEs
                xcorr_actual_PBE(pbe, :) = xcorr(activationStrPBEs.fineTS.(currPeriod){pbe}(ii, :), activationStrPBEs.fineTS.(currPeriod){pbe}(jj, :), maxLag);
            end

            xcorr_actual.(currPeriod){ii, jj} = sum(xcorr_actual_PBE, 1);
            
            % temporal bias
            A = sum(xcorr_actual.(currPeriod){ii, jj}(1:maxLag));
            B = sum(xcorr_actual.(currPeriod){ii, jj}(maxLag+2:2*maxLag+1));
            temporalBias.(currPeriod)(ii, jj) = (A-B)/(A+B);


            %%%%%%%%%%%%%%
            %%% SHUFFLE %%
            
            % cross-correlations(null distributions)
            
            xcorr_null_PBE = zeros(nShuffles, 2*maxLag+1, nPBEs);
            
            for pbe = 1:nPBEs

                s1 = activationStrPBEs.fineTS.(currPeriod){pbe}(ii, :);
                s2 = activationStrPBEs.fineTS.(currPeriod){pbe}(jj, :);

                nTimeBins = size(s1, 2);

                temp = zeros(nShuffles, 2*maxLag+1);
                for ishuffle = 1:nShuffles % why parfor is so slow with this??
                    temp(ishuffle, :) = xcorr(s1, s2(randperm(nTimeBins)), maxLag);
                end

                xcorr_null_PBE(:,:, pbe) = temp;    
            end

            xcorr_null.(currPeriod){ii, jj} = sum(xcorr_null_PBE, 3);
            
            
            % temporal bias (null distributions)
            temporalBias_null.(currPeriod){ii, jj} = zeros(nShuffles, 1);

            for ishuffle = 1:nShuffles

                A = sum(xcorr_null.(currPeriod){ii, jj}(ishuffle, 1:maxLag));
                B = sum(xcorr_null.(currPeriod){ii, jj}(ishuffle, maxLag+2:2*maxLag+1));

                temporalBias_null.(currPeriod){ii, jj}(ishuffle) = (A-B)/(A+B);    
            end  

        end
    end
end


%%%%%%%%%%%%
%% figures %


%% Eigenvalues, IC bases (weights) and spatial tunings

for its = 2

    currTimeScale = timeScaleNames{its};
    currEigValues = eigenValues.(currTimeScale);
    currEVthresh  = evThresh.(currTimeScale);
    currICweights = ICweights.(currTimeScale);

    figure;
    set(gcf, 'position', [2500 700 620 750])

    % eigenValues
    subplot(2,5,[1 2])

    hold on
    plot(1:numel(currEigValues), currEigValues, 'linewidth', 2, 'color', 'k')

    signPCidx  = find(currEigValues > currEVthresh);
    insigPCidx = setdiff(1:numel(currEigValues), signPCidx); 

    plot(signPCidx, currEigValues(signPCidx), 's', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'markersize', 3)
    plot(insigPCidx, currEigValues(insigPCidx), 'o', 'markerfacecolor', 'none', 'markeredgecolor', 'k', 'markersize',1 )
    line([0 numel(currEigValues)], [currEVthresh currEVthresh], 'linestyle','--', 'linewidth', 1, 'color','r')

    hold off

    set(gca, 'fontsize', 10, 'box', 'off', 'linewidth', 1)
    ylabel('eigenvalues', 'fontsize', 12)
    xlabel('principal components', 'fontsize', 12)
    
    title(currTimeScale, 'fontsize', 12)

    % IC bases (weights)
    subplot(2,5, 3:5)
    cm = redblue;

    imagesc(currICweights(:, runTemplate_RL)', [-max(abs(currICweights(:))) max(abs(currICweights(:)))]); colormap(cm); colorbar; 
    set(gca, 'ydir', 'normal', 'fontsize', 10, 'linewidth', 1)
    xlabel('independent components', 'fontsize', 12)
    ylabel('units', 'fontsize', 12)


    % sptial tuning of the ICs
    subplot(2,5,8:10)

    ICtunings_norm = zscore(ICtunings.(currTimeScale)(:, 2:end), [], 2); % excluding the zero position bin and z scoring
    imagesc(1:size(ICtunings_norm, 2)*2, 1:nICs,  ICtunings_norm); colormap(cm); colorbar;
    set(gca, 'ydir', 'normal', 'fontsize', 12)
    xlabel('track position(cm)', 'fontsize', 12)
    ylabel('independent components', 'fontsize', 12)

end


%% ICs tunings based on BD-based reconstrcuted positions


figure;
set(gcf, 'units', 'centimeters', 'position', [100 100 20 6])

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    subplot(1,3,iperiod)
    
    imagesc(1:nICs, (1:nPosBins)*posBinSize , ICAvgPosterior.(currPeriod)(:, ICpositionsortInd.RUN), [0 0.05]); colormap(cm) % , [-0.2 0.2]
    
    set(gca, 'xTickLabels', ICpositionsortInd.(currPeriod))
    
    title(currPeriod, 'fontsize', 12)
    xlabel('ICs', 'fontsize', 10)
    ylabel('track position', 'fontsize', 10)
    
    axis square
end


%% Reactivation profile of ICs during sleep periods (PRE & POST)

for its = 1
    
    currTimeScale = timeScaleNames{its};
    
    nor = 8; % change these if the number of ICs is higher than 32
    noc = 4;
    
    plotICactivationStrengths(activationStrength.(currTimeScale), binCenters.(currTimeScale), behavior, nor, noc, currTimeScale)
end


%% Totoal reactivation strengths by summing over ICs (PRE & POST)

for its = 2
    
    currTimeScale = timeScaleNames{its};
    
    nor = 1;
    noc = 1;

    plotICactivationStrengths(overallActivationStrength.(currTimeScale), binCenters.(currTimeScale), behavior, 1, 1, [currTimeScale '_overall'])
end


%% visualizing spike rasters, posterior probability and reactivation strengh of different ICs within individual PBEs

for iperiod = 1:3
    currentPeriod = periodNames{iperiod};
    
    directory = fullfile(mainDir, 'ICA', currentPeriod);
    if ~exist(directory, 'dir')
        mkdir(directory)
    end

    currBDscores = BDseqscore.(currentPeriod).data.unitIDshuffle.weightedCorr.prctilescore;
    
    
    % PBEs with high BD scores
    PBE2Plot = find(currBDscores > 95);
    PBE2Plot = PBE2Plot(randperm(numel(PBE2Plot), 40));

    filename  = ['ICsActivationDuringPBEs_' currentPeriod 'highBDscores'];
    textOnTop = sprintf('Session: %s\nPeriod: %s', sessionName, currentPeriod);
    
    plotICactivation4PBEs(binnedFiringPBEs.(currentPeriod), activationStrPBEs.fineTS.(currentPeriod), posteriorProbMatrix.(currentPeriod).data, currBDscores, runTemplate_LR, runTemplate_RL, PBE2Plot, directory, filename, textOnTop)


    % randomly selected PBEs
    PBE2Plot = randperm(numel(currBDscores), 40);
    
    filename  = ['ICsActivationDuringPBEs_' currentPeriod 'randomPBEs'];
    textOnTop = sprintf('Session: %s\nPeriod: %s', sessionName, currentPeriod);
    
    plotICactivation4PBEs(binnedFiringPBEs.(currentPeriod), activationStrPBEs.fineTS.(currentPeriod), posteriorProbMatrix.(currentPeriod).data, currBDscores, runTemplate_LR, runTemplate_RL, PBE2Plot, directory, filename, textOnTop)

end 


%% IC cross-correlograms


ICs2Plot  = randperm(nICs, 10);
ICs2Plot = sort(ICs2Plot, 'ascend');
nICs2Plot = numel(ICs2Plot);

lags = -maxLag:maxLag;

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    
    figure;
    for ii = 1:nICs2Plot
        
        flag = 0;
        
       for jj = ii:nICs2Plot

            subplot(nICs2Plot, nICs2Plot, (ii-1)*nICs2Plot+jj)
            set(gca, 'fontsize', 6)
            hold on

            if jj == ii
                bar(lags, xcorr_actual.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}, 'FaceColor', 'b', 'EdgeColor', 'b')
            else
                bar(lags, xcorr_actual.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}, 'FaceColor', 'k', 'EdgeColor', 'k')
            end

            plot(lags, median(xcorr_null.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}, 1), 'color', 'm', 'linewidth', 1)
            plot(lags, prctile(xcorr_null.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}, 95, 1), 'color', 'm', 'linestyle', ':')
            plot(lags, prctile(xcorr_null.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}, 5, 1), 'color', 'm', 'linestyle', ':')


            hold off

            xlim([-20 20])

            if ii == 1
                title(sprintf('IC #%d', ICs2Plot(jj)), 'FontWeight', 'normal');
            end

            if flag == 0 
               xlabel('time bin')
               ylabel({sprintf('IC #%d', ICs2Plot(ii));'';'cross-corr.'})
               flag =1;
            end

            xticks([-20 -10  0  10  20])

            yLimit = ylim;

            pval = length(find(abs(temporalBias_null.(currPeriod){ICs2Plot(ii), ICs2Plot(jj)}) >= abs(temporalBias.(currPeriod)(ICs2Plot(ii), ICs2Plot(jj)))))/(nShuffles+1);
            text(5, yLimit(1)+0.8*range(yLimit), {sprintf('tb= %.2f', temporalBias.(currPeriod)(ICs2Plot(ii), ICs2Plot(jj))); sprintf('pval= %0.2f', pval)}, 'fontsize', 7)

       end
    end
    
    subplot(nICs2Plot, nICs2Plot, (7-1)*nICs2Plot+2)
    text(0, 0, currPeriod, 'fontsize', 16)
    set(gca, 'XColor', 'none', 'YColor', 'none')
    
end


%% ICs temporal biases for each period and correlation between ICs' temporal biases in differen periods


figure; 
set(gcf, 'position', [100 100 900 1200])


% temporal bias in individual periods

colorLimit = [-0.1 0.1];

for iperiod = 1:3
    currPeriod = periodNames{iperiod};
    
    subplot(4,3, iperiod+3)
    
    temp = temporalBias.(currPeriod);
    temp(isnan(temp)) = 0;
    
    imagesc(temp, colorLimit); colormap(cm)
    set(gca, 'linewidth', 2)
    title({currPeriod; '';'IC j'}, 'fontsize', 12, 'fontweight', 'normal')
    ylabel('IC i', 'fontsize', 12)
    axis square
end

h = colorbar;
h.Position(1) = h.Position(1) + 0.08; 
h.TickLabels{1} = ['< ' num2str(colorLimit(1))];
h.TickLabels{end} = ['> ' num2str(colorLimit(2))];



% scatter plots of temporal bias

RUNtemporalBiasAll  = temporalBias.RUN(~isnan(temporalBias.RUN));
PREtemporalBiasAll  = temporalBias.PRE(~isnan(temporalBias.PRE));
POSTtemporalBiasAll = temporalBias.POST(~isnan(temporalBias.POST));


subplot(4,3,7)

set(gca, 'linewidth', 2, 'box', 'off')
plot(RUNtemporalBiasAll, PREtemporalBiasAll, '.k', 'markersize', 3)
x = [ones(length(RUNtemporalBiasAll), 1) RUNtemporalBiasAll];
b= x\PREtemporalBiasAll;
ycalc = x*b;
hold on
plot(RUNtemporalBiasAll, ycalc, 'color', [255 233 64]/255, 'linewidth',2)
xlabel('RUN temporal bias', 'fontsize', 12)
ylabel('PRE temporal bias', 'fontsize', 12)
xlim([-0.3 0.3])
ylim([-0.3 0.3])
grid on
axis square
set(gca, 'linewidth', 2, 'box', 'off')
subplot(4,3,8)

set(gca, 'linewidth', 2, 'box', 'off')
plot(RUNtemporalBiasAll, POSTtemporalBiasAll, '.k', 'markersize', 3)
b= x\POSTtemporalBiasAll;
ycalc = x*b;
hold on
plot(RUNtemporalBiasAll, ycalc, 'color', [64 86 255]/255, 'linewidth',2)
xlabel('RUN temporal bias', 'fontsize', 12)
ylabel('POST temporal bias', 'fontsize', 12)
xlim([-0.3 0.3])
ylim([-0.3 0.3])
grid on
axis square
set(gca, 'linewidth', 2, 'box', 'off')

% calculating correlation 
RUN_PRE_tbCorr  = zeros(nICs, 1);
RUN_POST_tbCorr = zeros(nICs, 1);

RUN_POST_tbCorr_shuffle = zeros(nICs, 1000);
RUN_PRE_tbCorr_shuffle  = zeros(nICs, 1000);

for ii = 1: nICs
    
    PREtbs = [temporalBias.PRE(ii, min(nICs, ii+1): nICs) temporalBias.PRE(1:ii-1, ii)'];
    PREtbs = PREtbs./abs(PREtbs);
    
    RUNtbs = [temporalBias.RUN(ii, min(nICs, ii+1): nICs) temporalBias.RUN(1:ii-1, ii)'];
    RUNtbs = RUNtbs./abs(RUNtbs);
    
    POSTtbs = [temporalBias.POST(ii, min(nICs, ii+1): nICs) temporalBias.POST(1:ii-1, ii)'];
    POSTtbs = POSTtbs./abs(POSTtbs);
    
    
    RUN_PRE_tbCorr(ii)  = abs((PREtbs  * RUNtbs')/norm(PREtbs) /norm(RUNtbs)); % cosine similarity
    RUN_POST_tbCorr(ii) = abs((POSTtbs * RUNtbs')/norm(POSTtbs)/norm(RUNtbs));
    
    for ishuffle = 1:1000
        
        rndgen = randperm(numel(RUNtbs));
        RUN_POST_tbCorr_shuffle(ii, ishuffle) = abs((POSTtbs * RUNtbs(rndgen)')/norm(POSTtbs)/norm(RUNtbs(rndgen)));  
        RUN_PRE_tbCorr_shuffle(ii, ishuffle)  = abs((PREtbs * RUNtbs(rndgen)')/norm(PREtbs)/norm(RUNtbs(rndgen)));
    end
        
end


subplot(4,3,10)

pooled = [RUN_PRE_tbCorr; RUN_POST_tbCorr; RUN_PRE_tbCorr_shuffle(:); RUN_POST_tbCorr_shuffle(:)];
bins = linspace(min(pooled), max(pooled), 20);

hPRE  = hist(RUN_PRE_tbCorr, bins);
hPRE  = hPRE./sum(hPRE);

hPOST = hist(RUN_POST_tbCorr, bins);
hPOST = hPOST./sum(hPOST); 

hPRE_shuffle = hist(RUN_PRE_tbCorr_shuffle(:), bins);
hPRE_shuffle = hPRE_shuffle./sum(hPRE_shuffle);
hPOST_shuffle = hist(RUN_POST_tbCorr_shuffle(:), bins);
hPOST_shuffle = hPOST_shuffle./sum(hPOST_shuffle);


bar(bins, hPRE, 'FaceColor', [255 233 64]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPOST, 'FaceColor', [64 86 255]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
legend('RUN-PRE', 'RUN-POST', 'location', 'northeast', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')
axis square
xlim([-0.1 1])


subplot(4,3,11)

bar(bins, hPRE_shuffle, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPRE, 'FaceColor', [255 233 64]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
title('RUN-PRE', 'fontsize', 12)
legend('shuffle', 'data', 'location', 'northeast', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')
axis square
xlim([-0.1 1])

subplot(4,3,12)

bar(bins, hPOST_shuffle, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6)
hold on
bar(bins, hPOST, 'FaceColor', [64 86 255]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.6)
xlabel({'correlation'; 'with RUN temporal biases'}, 'fontsize', 12)
ylabel('probability', 'fontsize', 12)
title('RUN-POST', 'fontsize', 12)
legend('shuffle', 'data', 'location', 'northeast', 'box', 'off')
set(gca, 'linewidth', 2, 'box', 'off')
axis square
xlim([-0.1 1])

%%%%%%%%%%%
% functions


function activationStrength = calActStrength(ICweights, x, meanCounts, stdCounts)


nTimebins = size(x,2);

silentBins = sum(x,1) == 0;

% x = zscore(x, [], 2);
x = (x - repmat(meanCounts, [1 nTimebins]))./repmat(stdCounts, [1 nTimebins]);
x(isnan(x)) = 0;

nICs = size(ICweights, 1); % rows correspond to weights for each IC

activationStrength = zeros(nICs, nTimebins);

for ii = 1:nICs
    currentICbase = ICweights(ii, :)';
    ProjectorOp = currentICbase * currentICbase';
    ProjectorOp = ProjectorOp - ProjectorOp.*eye(size(ProjectorOp)); % removing the diagonal to prevent high activation strength caused by the isolated activity of a single neuron with high weight in the pattern
    
    for jj = 1:nTimebins
        activationStrength(ii, jj) = x(:, jj)' * ProjectorOp * x(:, jj);
%         activationStrength(ii, jj) = x(:, jj)' * currentICbase;
    end
end

% setting reactivating strength in silent bins zero

activationStrength(:, silentBins) = 0;

end


function smoothedMat = arrayConv(inputMat, avgWindow)

smoothedMat = zeros(size(inputMat));
for ii = 1: size(inputMat, 1)
    smoothedMat(ii, :) = conv(inputMat(ii, :), avgWindow, 'same');
    
    % dealing with the edge effect by replacing the first and last segements of the data (segments equal in length with half length of  avgWindow) with Nans 
    avgWindowlen = length(avgWindow);
    
    smoothedMat(ii, [1:floor(avgWindowlen/2) (end-floor(avgWindowlen/2)):end]) = nan;   
    
end

end


function z = zscore_KM(x, mu, sigma)

z = (x - repmat(mu, [1 size(x, 2)]))./repmat(sigma, [1 size(x, 2)]);

end
