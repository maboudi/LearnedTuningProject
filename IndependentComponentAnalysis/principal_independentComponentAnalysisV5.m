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


posBinSize = 2; % in cm
runSpeedThresh = 0;
binOverlapRatio = 0;

for its = 1:2  
    currTimeScale = timeScaleNames{its};
    runTraining.(currTimeScale) = binRunData_1D(spikeStruct, behavior.time(2, :), behavior, speed, runSpeedThresh, binDur.(currTimeScale), binOverlapRatio, posBinSize, fileinfo);
end


%% for calculating the sptial tuning of the ICs

binOverlapRatio = struct('coarseTS', 0, 'fineTS', 0);

binnedRun4Tunning = struct('coarseTS', [], 'fineTS', []);
posbinIdx     = struct('coarseTS', [], 'fineTS', []);
linposcenters = struct('coarseTS', [], 'fineTS', []);
binCenters4Tuning   = struct('coarseTS', [], 'fineTS', []);

for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    
    currBinDur          = binDur.(currTimeScale);
    currBinOverlapRatio = binOverlapRatio.(currTimeScale);
    currStepSize        = (1-currBinOverlapRatio) * currBinDur;
    
    binCenters4Tuning.(currTimeScale)    = (behavior.time(2, 1):currStepSize:(behavior.time(2, 2)-currBinDur)) + currBinDur/2;
    [binnedRun4Tunning.(currTimeScale), posbinIdx.(currTimeScale), linposcenters.(currTimeScale)] = binRunData_1D(spikeStruct, behavior.time(2, :), behavior, speed, runSpeedThresh, currBinDur, currBinOverlapRatio, posBinSize, fileinfo);
end



%%%%% test data (the whole period is binned. Later we add analyses limited to PBEs as well)
%%

tic
qclus = [1 2 3];

binOverlapRatio = struct('coarseTS', 0.5, 'fineTS', 0);

basicStruct  = struct('PRE', [], 'RUN', [], 'POST', []);
binCenters   = struct('coarseTS', basicStruct, 'fineTS', basicStruct);
binnedFiring = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

binnedFiring_nonOverlap = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

for iperiod = 1:3
    
    for its = 1:2 
        currPeriod          = periodNames{iperiod};
        currTimeScale       = timeScaleNames{its};
        currBinDur          = binDur.(currTimeScale);
        currBinOverlapRatio = binOverlapRatio.(currTimeScale);
        currStepSize        = (1-currBinOverlapRatio) * currBinDur;

        
        binCenters.(currTimeScale).(currPeriod)    = (behavior.time(iperiod, 1):currStepSize:(behavior.time(iperiod, 2)-currBinDur)) + currBinDur/2;
        binnedFiring.(currTimeScale).(currPeriod)  = timeBinning_withOverlap(behavior.time(iperiod, :) , spikeStruct, qclus, currBinDur, currBinOverlapRatio, fileinfo); % the theta periods        
    end
    
    currBinDur = binDur.coarseTS;
    currBinOverlapRatio = 0;
    
    binnedFiring_nonOverlap.coarseTS.(currPeriod) = timeBinning_withOverlap(behavior.time(iperiod, :) , spikeStruct, qclus, currBinDur, currBinOverlapRatio, fileinfo);
    
    binnedFiring_nonOverlap.fineTS.(currPeriod)   = binnedFiring.fineTS.(currPeriod);
end

toc



%%
%%%%% mean and standrad deviation of firing rates for the same z-scoring throughout the analysis 

firingMean = struct('coarseTS', basicStruct, 'fineTS', basicStruct);
firingStd  = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

for its = 1:2
    for iperiod = 1:3
        
        currTimeScale = timeScaleNames{its}; 
        currPeriod    = periodNames{iperiod};
        
        currData      = binnedFiring_nonOverlap.(currTimeScale).(currPeriod){1,2};

        firingMean.(currTimeScale).(currPeriod) = mean(currData, 2); % in Hz
        firingStd.(currTimeScale).(currPeriod)  = std(currData, 0, 2); % in Hz
    end
end


%%
%%% detecting time bins that overlap with PBEs (note that the PBEs that we
%%% defined earlier are not all the PBEs). So, one way is to calcualte the
%%% spike counts within a time bin and compare it with reactivations in those time bins. 

nSpikesCoarseTBs = struct('coarseTS', basicStruct, 'fineTS', basicStruct);
zfiringCoarseTBs = struct('coarseTS', basicStruct, 'fineTS', basicStruct);

for its = 1:2
    currTimeScale = timeScaleNames{its};
    for iperiod = 1:3

        currPeriod = periodNames{iperiod};

        currBinDur = binDur.(currTimeScale);

        nSpikesCoarseTBs.(currTimeScale).(currPeriod) = sum(binnedFiring.(currTimeScale).(currPeriod){1,2}, 1);

        zfiringCoarseTBs.(currTimeScale).(currPeriod) = (nSpikesCoarseTBs.(currTimeScale).(currPeriod) - mean(nSpikesCoarseTBs.(currTimeScale).(currPeriod)))/std(nSpikesCoarseTBs.(currTimeScale).(currPeriod));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the RUN PCs and PC weights
% 
% PCs         = struct('coarseTS', [], 'fineTS', []);
% PCweights   = struct('coarseTS', [], 'fineTS', []);
% eigenValues = struct('coarseTS', [], 'fineTS', []);
% evThresh    = struct('coarseTS', [], 'fineTS', []); % eigen value threshold
% 
% 
% for its = 1:2
%     
%     currTimeScale = timeScaleNames{its};
%     
%     meanFiringCount = firingMean.(currTimeScale).RUN;
%     stdFiringCount  = firingStd.(currTimeScale).RUN;
% 
%     [PCs.(currTimeScale), PCweights.(currTimeScale), eigenValues.(currTimeScale), evThresh.(currTimeScale)] = pca_ica_v3(runTraining.(currTimeScale), 'firingmean', meanFiringCount, 'firingstd', stdFiringCount, 'onlyPCA', 1);
%     nPCs = size(PCs.(currTimeScale), 1);
% 
%     for ii = 1:nPCs
%         PCweights.(currTimeScale)(ii, :) = PCweights.(currTimeScale)(ii, :)/norm(PCweights.(currTimeScale)(ii, :));
%         
%         %signs are set such that the highest absolute weight of each pattern is
%         % always positive
%         
%         [~, maxAbsWeightIndex] = max(abs(PCweights.(currTimeScale)(ii, :)));
%         maxWeightSign = sign(PCweights.(currTimeScale)(ii, maxAbsWeightIndex));
%         PCweights.(currTimeScale)(ii, :) = PCweights.(currTimeScale)(ii,:) * maxWeightSign;
%         
%     end    
% end




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

    [ICs.(currTimeScale), ICweights.(currTimeScale), eigenValues.(currTimeScale), evThresh.(currTimeScale)] = pca_ica_v3(runTraining.(currTimeScale), 'firingmean', meanFiringCount, 'firingstd', stdFiringCount, 'onlyPCA', 1);
    nICs = size(ICs.(currTimeScale), 1);

    for ii = 1:nICs
        ICweights.(currTimeScale)(ii, :) = ICweights.(currTimeScale)(ii, :)/norm(ICweights.(currTimeScale)(ii, :));
        
        %signs are set such that the highest absolute weight of each pattern is
        % always positive
        
%         [~, maxAbsWeightIndex] = max(abs(ICweights.(currTimeScale)(ii, :)));
%         maxWeightSign = sign(ICweights.(currTimeScale)(ii, maxAbsWeightIndex));
%         ICweights.(currTimeScale)(ii, :) = ICweights.(currTimeScale)(ii,:) * maxWeightSign;
        
    end    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate ICs in jitered data


%%% generating temporally shuffled data using two methods; (1) jittering
%%% time bins locally, (2) shuffling the time bin throughout the RUN block

shuffleBinIdx = struct('coarseTS', [], 'fineTS', []);

shuffleMethod = 2;

if shuffleMethod == 1


    % jittering time bins in a limited range independently for each unit


    jitterLength.fineTS   = 20;
    jitterLength.coarseTS = 5;

    nUnits = size(runTraining.coarseTS, 1);
   

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

elseif shuffleMethod == 2

    % Randomly shuffle the time bins across the whole period, again
    % independently for each unit

    nUnits = size(runTraining.coarseTS, 1);
    shuffleBinIdx = struct('fineTS', [], 'coarseTS', []);

    for its = 1:2

        currTimeScale = timeScaleNames{its};
        nTimeBins = size(runTraining.(currTimeScale), 2);

        shuffleBinIdx.(currTimeScale) = zeros(nUnits, nTimeBins);
        for iUnit = 1:nUnits

            shuffleBinIdx.(currTimeScale)(iUnit, :) = randperm(nTimeBins);
        end

    end

end


% calculating surrogate ICs

ICs_s         = struct('coarseTS', [], 'fineTS', []);
ICweights_s   = struct('coarseTS', [], 'fineTS', []);
eigenValues_s = struct('coarseTS', [], 'fineTS', []);
evThresh_s    = struct('coarseTS', [], 'fineTS', []); % eigenvalue threshold


for its = 1:2
    
    currTimeScale  = timeScaleNames{its};
    currShuffleIdx = shuffleBinIdx.(currTimeScale);
    
    meanFiringCount = firingMean.(currTimeScale).RUN;
    stdFiringCount  = firingStd.(currTimeScale).RUN;
    
    shuffledTrainData = zeros(size(runTraining.(currTimeScale)));
    for iUnit = 1:nUnits 
        shuffledTrainData(iUnit, :) = runTraining.(currTimeScale)(iUnit, currShuffleIdx(iUnit, :));
    end
    

    nTopComp = size(ICweights.(currTimeScale), 1);
    [ICs_s.(currTimeScale), ICweights_s.(currTimeScale), eigenValues_s.(currTimeScale), evThresh_s.(currTimeScale)] = pca_ica_v3(shuffledTrainData, 'firingmean', meanFiringCount, 'firingstd', stdFiringCount, 'nComp', nTopComp);

    for ii = 1:size(ICs_s.(currTimeScale), 1)
        ICweights_s.(currTimeScale)(ii, :) = ICweights_s.(currTimeScale)(ii, :)/norm(ICweights_s.(currTimeScale)(ii, :));
        
        %signs are set such that the highest absolute weight of each pattern is
        % always positive
        
        [~, maxAbsWeightIndex] = max(abs(ICweights_s.(currTimeScale)(ii, :)));
        maxWeightSign = sign(ICweights_s.(currTimeScale)(ii, maxAbsWeightIndex));
        ICweights_s.(currTimeScale)(ii, :) = ICweights_s.(currTimeScale)(ii,:) * maxWeightSign;
        
    end    
       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate (re)activation strengths of ICs

initStruct = struct('coarseTS', [], 'fineTS', []);
initStruct = structfun(@(x) struct('PRE', [], 'RUN', [], 'POST', []), initStruct, 'UniformOutPut', false);

[activationStrength, overallActivationStrength, ...
activationStrength_s, overallActivationStrength_s] = deal(initStruct);


for its = 1:2

    currTimeScale = timeScaleNames{its};

    activationStrength.(currTimeScale)          = structfun(@(x) struct('original', [], 'smoothed', []), activationStrength.(currTimeScale), 'UniformOutPut', false);
    overallActivationStrength.(currTimeScale)   = structfun(@(x) struct('original', [], 'smoothed', []), overallActivationStrength.(currTimeScale), 'UniformOutPut', false);
    
    activationStrength_s.(currTimeScale)        = structfun(@(x) struct('original', [], 'smoothed', []), activationStrength_s.(currTimeScale), 'UniformOutPut', false);
    overallActivationStrength_s.(currTimeScale) = structfun(@(x) struct('original', [], 'smoothed', []), overallActivationStrength_s.(currTimeScale), 'UniformOutPut', false);

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


        % overall activation strength by summing over all ICs
        activationStrength.(currTimeScale).(currPeriod).original  = calActStrength(ICweights.(currTimeScale), binnedFiring.(currTimeScale).(currPeriod){1,2}, meanFiringCount, stdFiringCount);
        activationStrength.(currTimeScale).(currPeriod).smoothed  = arrayConv(activationStrength.(currTimeScale).(currPeriod).original, avgWindow); % running average
       

        activationStrength_s.(currTimeScale).(currPeriod).original  = calActStrength(ICweights_s.(currTimeScale), binnedFiring.(currTimeScale).(currPeriod){1,2}, meanFiringCount, stdFiringCount);
        activationStrength_s.(currTimeScale).(currPeriod).smoothed  = arrayConv(activationStrength_s.(currTimeScale).(currPeriod).original, avgWindow); % running average




        % overall activation strength by summing over all ICs
        overallActivationStrength.(currTimeScale).(currPeriod).original  = mean(activationStrength.(currTimeScale).(currPeriod).original, 1);
        overallActivationStrength.(currTimeScale).(currPeriod).smoothed  = arrayConv(overallActivationStrength.(currTimeScale).(currPeriod).original, avgWindow); % running average
        
    
        overallActivationStrength_s.(currTimeScale).(currPeriod).original  = mean(activationStrength_s.(currTimeScale).(currPeriod).original, 1);
        overallActivationStrength_s.(currTimeScale).(currPeriod).smoothed  = arrayConv(overallActivationStrength_s.(currTimeScale).(currPeriod).original, avgWindow); % running average
        
    end
end


%% z-score using PRE-distribution
 
for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    
    PREmean           = nanmean(activationStrength.(currTimeScale).PRE.smoothed, 2);
    PREstd            = nanstd(activationStrength.(currTimeScale).PRE.smoothed, 0, 2);

    PREmean_s         = nanmean(activationStrength_s.(currTimeScale).PRE.smoothed, 2);
    PREstd_s          = nanstd(activationStrength_s.(currTimeScale).PRE.smoothed, 0, 2);

    PREmeanOverall    = nanmean(overallActivationStrength.(currTimeScale).PRE.smoothed, 2);
    PREstdOverall     = nanstd(overallActivationStrength.(currTimeScale).PRE.smoothed, 0, 2);

    PREmeanOverall_s  = nanmean(overallActivationStrength_s.(currTimeScale).PRE.smoothed, 2);
    PREstdOverall_s   = nanstd(overallActivationStrength_s.(currTimeScale).PRE.smoothed, 0, 2);

    
    for iperiod = 1:3
        currPeriod = periodNames{iperiod};

        activationStrength.(currTimeScale).(currPeriod).smoothed          = zscore_KM(activationStrength.(currTimeScale).(currPeriod).smoothed, PREmean, PREstd);
        overallActivationStrength.(currTimeScale).(currPeriod).smoothed   = zscore_KM(overallActivationStrength.(currTimeScale).(currPeriod).smoothed, PREmeanOverall, PREstdOverall);
        
        activationStrength_s.(currTimeScale).(currPeriod).smoothed        = zscore_KM(activationStrength_s.(currTimeScale).(currPeriod).smoothed, PREmean_s, PREstd_s);
        overallActivationStrength_s.(currTimeScale).(currPeriod).smoothed = zscore_KM(overallActivationStrength_s.(currTimeScale).(currPeriod).smoothed, PREmeanOverall_s, PREstdOverall_s);

    end
end


%% calculating reactivation strength of ICs during RUN overlapping time periods

activationStrength_tuning = struct('coarseTS', [], 'fineTS', []);

for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    
    meanFiringCount = firingMean.(currTimeScale).RUN;
    stdFiringCount  = firingStd.(currTimeScale).RUN;

    activationStrength_tuning.(currTimeScale) =  calActStrength(ICweights.(currTimeScale), binnedRun4Tunning.(currTimeScale), meanFiringCount, stdFiringCount);
end


%% Measuring the degree to which the reactivations are correlated with firing rates (or occur during PBEs)  


currPeriod = 'PRE';
currTimeScale = 'coarseTS';

nICs  = size(ICweights.(currTimeScale), 1);

figure; 
for ic = 1:nICs
    
    subplot(8,4, ic)
    set(gca, 'tickDir', 'out', 'fontsize', 8)
    
    currActivations = activationStrength.(currTimeScale).(currPeriod).original;
    
    scatter(currActivations(ic, :), zfiringCoarseTBs.(currTimeScale).(currPeriod), 3, 'k', 'filled')
    hold on
    line([20 20], [-1 10], 'color', 'r', 'linestyle', '-')
    hold on
    line([-20 -20], [-1 10], 'color', 'r', 'linestyle', '-')
    
    xlim([-50 500])
    
    if ic >= 21
        xlabel('IC reactivation strength', 'fontsize', 10)
    end
    
    if mod(ic, 4) == 1
        ylabel('MUA(z)', 'fontsize', 10)
    end
    
    
    text(-40, 11, ['IC ' num2str(ic)], 'fontsize', 10, 'fontweight', 'normal')
    if ic == 1
        text(-20, 15, currPeriod, 'fontsize', 18, 'fontweight', 'bold')
    end
        
end


%%
% % measure to what extent the IC peak activations happen during the PBEs

periodNames = {'PRE'; 'RUN'; 'POST'};
popBurstTimeCourse = struct('PRE', sdat_pre, 'RUN', sdat_run, 'POST', sdat_post);

nICs = size(ICweights.coarseTS, 1);

ICpeakReactivationPoints =  cell(nICs, 3); % columns corresponding to PRE, RUN, and POST

MUAactivationCorr = struct('PRE', [], 'RUN', [], 'POST', []);

for iperiod = 1:3
    
    figure;
    
for ic = 1: nICs
    
    
    subplot(4,8,ic)
    currPeriod      = periodNames{iperiod};
    currActStrength = activationStrength.coarseTS.(currPeriod).original(ic, :);
    
%     highActPnts = find(currActStrength > prctile(currActStrength, 99.9));
    highActPnts = find(currActStrength > 20);
    
    
%     ICpeakReactivationPoints(ic).(currPeriod) = floor((binCenters.coarseTS.(currPeriod)(highActPnts)-behavior.time(iperiod, 1))*1000); % in MILISEC
    ICpeakReactivationPoints{ic, iperiod} = floor((binCenters.coarseTS.(currPeriod)(highActPnts)-behavior.time(iperiod, 1))*1000); % in MILISEC
    
    
    muaAThiActivations    = popBurstTimeCourse.(currPeriod)(floor((binCenters.coarseTS.(currPeriod)(highActPnts)-behavior.time(iperiod, 1))*1000));
    actStrAThiActivations = currActStrength(highActPnts);
    
    MUAactivationCorr.(currPeriod)(ic) = corr(muaAThiActivations, actStrAThiActivations');
    
    
    muaATrandom        = popBurstTimeCourse.(currPeriod)(floor((binCenters.coarseTS.(currPeriod)(highActPnts)-behavior.time(iperiod, 1))*1000) + 1000);
    
    
    errorbar(1, mean(muaAThiActivations), std(muaAThiActivations), '-s', 'MarkerSize', 6, ...
        'markerEdgeColor', 'k', 'markerFaceColor', 'k', 'color', 'k', 'linewidth', 2) 
    
    hold on
    
    errorbar(2, mean(muaATrandom), std(muaATrandom), '-s', 'MarkerSize', 6, ...
        'markerEdgeColor', [0.5 0.5 0.5], 'markerFaceColor', [0.5 0.5 0.5], 'color', [0.5 0.5 0.5], 'linewidth', 2) 
    
    xlim([0 3])
    ylim([-1.5 6])
    text(0, 6.5, ['IC ' num2str(ic)], 'fontsize', 10, 'fontweight', 'normal')
    
    if ic == 1
        text(0, 8, currPeriod, 'fontsize', 14, 'fontweight', 'bold')
    end
    
    ylabel('MUA(z)')
    
    set(gca, 'xtick', [1 2], 'xticklabel', [])
    if ic >= 17 
        set(gca, 'xtick', [1 2], 'xticklabel', {'peak activations'; 'random pnts'}, 'xticklabelrotation', 45)
    end
    
    grid on
    
    hold on
    line([0 3], [3 3], 'color', 'r', 'linewidth', 2)
    text(0.2, 3.5, ['\uparrow' num2str(numel(find(muaAThiActivations >= 3)))], 'color', 'g')
    text(0.2, 2.5, ['\downarrow' num2str(numel(find(muaAThiActivations <  3)))], 'color', 'r')
    
    hold on
    line([0 3], [2 2], 'color', 'r', 'linestyle', ':', 'linewidth', 2)
    text(1.2, 2.5, ['\uparrow' num2str(numel(find(muaAThiActivations >= 2)))], 'color', 'k')
    text(1.2, 1.5, ['\downarrow' num2str(numel(find(muaAThiActivations <  2)))], 'color', 'k')
    
    
    
end

end


%% extracting boundaries of events centered on peak IC reactivations, based on population burst

popBurstTimeCourse = struct('PRE', sdat_pre, 'RUN', sdat_run, 'POST', sdat_post);

nICs = size(ICweights.coarseTS, 1);
selectICs = 1:nICs;


hiActEvts = findHiPCactPeriods(ICpeakReactivationPoints, popBurstTimeCourse, selectICs);


figure; 

subplot(3,1,[2 3])

clear currData

currData.PRE  = hiActEvts.PRE(:, 3);
% currData.RUN  = hiActEvts.RUN(:, 3);
currData.POST = hiActEvts.POST(:, 3);

violins = violinplot(currData); %[255 154 0]/255 

violins(1,1).ViolinColor = [35 110 150]/255;
% violins(1,2).ViolinColor = [100 100 100]/255;
violins(1,2).ViolinColor = [247 95 28]/255;


for ii=1:2

    violins(1,ii).ShowData = 0;
    violins(1,ii).ViolinPlot.LineWidth = 2;
    violins(1,ii).ViolinAlpha = 0.6;
end

ylabel('max MUA', 'fontsize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating reactivation strengths within PBEs

binnedFiringPBEs  = struct('PRE', [], 'RUN', [], 'POST', []);

binnedFiringPBEs.PRE  = PREbinnedPBEs.data;
binnedFiringPBEs.RUN  = RUNbinnedPBEs.data;
binnedFiringPBEs.POST = POSTbinnedPBEs.data;


activationStrPBEs = initStruct;

for its = 1:2

    currTimeScale = timeScaleNames{its};

    activationStrPBEs.(currTimeScale) = structfun(@(x) struct('actual', [], 'surrogate', []), activationStrPBEs.(currTimeScale), 'UniformOutPut', false);
    
end


PBEperiods = struct('PRE', [], 'RUN', [], 'POST', []);
PBEperiods.PRE  = secondaryPBEs_pre;
PBEperiods.RUN  = secondaryPBEs_run;
PBEperiods.POST = secondaryPBEs_post;


for its = 1:2
    
    currTimeScale = timeScaleNames{its};
    for iperiod = 1:3
        
        currPeriod = periodNames{iperiod};
        currPBEs = binnedFiringPBEs.(currPeriod);
        
        meanFiringCount = firingMean.(currTimeScale).(currPeriod);
        stdFiringCount  = firingStd.(currTimeScale).(currPeriod);

        
        nPBEs = size(currPBEs, 1);

        activationStrPBEs.(currTimeScale).(currPeriod).actual = cell(nPBEs, 1);
        activationStrPBEs.(currTimeScale).(currPeriod).surrogate = cell(nPBEs, 1);

        
        for pbe = 1:nPBEs
            
            if its == 1
                activationStrPBEs.(currTimeScale).(currPeriod).actual{pbe}    = calActStrength(ICweights.(currTimeScale), sum(currPBEs{pbe, 2}, 2), meanFiringCount, stdFiringCount);
                activationStrPBEs.(currTimeScale).(currPeriod).surrogate{pbe} = calActStrength(ICweights_s.(currTimeScale), sum(currPBEs{pbe, 2}, 2), meanFiringCount, stdFiringCount);
            
            elseif its == 2
                activationStrPBEs.(currTimeScale).(currPeriod).actual{pbe}    = calActStrength(ICweights.(currTimeScale), currPBEs{pbe, 2}, meanFiringCount, stdFiringCount);
                activationStrPBEs.(currTimeScale).(currPeriod).surrogate{pbe} = calActStrength(ICweights_s.(currTimeScale), currPBEs{pbe, 2}, meanFiringCount, stdFiringCount);

            end
                   
        end
    end
end




%% ICs activation scores

activationStrPBEs_scores = struct('summedStr', [], 'meanStr', [], 'maxStr', [], 'nActBins', []);

activationStrPBEs_scores.summedStr = initStruct;
activationStrPBEs_scores.meanStr   = initStruct;
activationStrPBEs_scores.maxStr    = initStruct;
activationStrPBEs_scores.nActBins  = initStruct;

scoreNames    = {'summedStr', 'meanStr', 'maxStr', 'nActBins'};
dataTypeNames = {'actual', 'surrogate'};

for its = 1:2

    currTimeScale = timeScaleNames{its};
    
    for iscore = 1:4
        activationStrPBEs_scores.(scoreNames{iscore}).(currTimeScale) = structfun(@(x) struct('actual', [], 'surrogate', []), activationStrPBEs_scores.(scoreNames{iscore}).(currTimeScale), 'UniformOutPut', false);
    end
end



for iperiod = 1:3

    currPeriod = periodNames{iperiod};
    currPBEs = binnedFiringPBEs.(currPeriod);

    nPBEs = size(currPBEs, 1);
    
    for its = 1:2
        
        currTimeScale = timeScaleNames{its};

        for iscore = 1:4
            
            activationStrPBEs_scores.(scoreNames{iscore}).(currTimeScale).(currPeriod) = struct('actual', [], 'surrogate', []);
            for itype = 1:2
                
                if itype == 1
                    activationStrengthThresh = prctile(activationStrength.(currTimeScale).(currPeriod).original(:), 99.99);
                else
                    activationStrengthThresh = prctile(activationStrength_s.(currTimeScale).(currPeriod).original(:), 99.99);
                end


                currDataType = dataTypeNames{itype};
            
                currActStrs = activationStrPBEs.(currTimeScale).(currPeriod).(currDataType);

                activationStrPBEs_scores.(scoreNames{iscore}).(currTimeScale).(currPeriod).(currDataType) = zeros(nPBEs, 1);

                for pbe = 1:nPBEs
                
                    PBEActStr = currActStrs{pbe};

                    switch iscore

                        case 1 
                            tempValue = sum(PBEActStr(PBEActStr > 0));
                        case 2
                            tempValue = mean(PBEActStr(PBEActStr > 0));
                        case 3

                            tempValue = max(PBEActStr(PBEActStr > 0));
                            if isempty(tempValue); tempValue = 0; end

                        case 4
                            cmpr2thresh =  PBEActStr > activationStrengthThresh;
                            tempValue   = length(find(sum(cmpr2thresh, 1)));
                            if itype == 2; tempValue = tempValue + 0.5; end
                    end
       
                    activationStrPBEs_scores.(scoreNames{iscore}).(currTimeScale).(currPeriod).(currDataType)(pbe) = tempValue;
                end
            end
        end
    end
end


%% scatter plots of IC activation scores - fineTS vs coarseTS for both actual and surrogate (temporal shuffle) ICs 

figure;

plotScoreNames = {'summed strength'; 'mean strength'; 'max strength'; 'no bins w/ sig. reactivation'};

for iscore = 1:4
    currScore = scoreNames{iscore};
    
    for iperiod = 1:3
        
        currPeriod = periodNames{iperiod};
        ax.(currPeriod) = subplot(4,3, (iscore-1)*3 + iperiod);

        for itype = [1 2]
            currDataType = dataTypeNames{itype};
            
            finTSData    = activationStrPBEs_scores.(scoreNames{iscore}).fineTS.(currPeriod).(currDataType);
            
            if iscore == 4
               coarseTSData = activationStrPBEs_scores.(scoreNames{3}).coarseTS.(currPeriod).(currDataType);
            else
               coarseTSData = activationStrPBEs_scores.(scoreNames{iscore}).coarseTS.(currPeriod).(currDataType);
            end
            
            if itype == 1; color = 'k';else; color = 'r'; end

%             plot(finTSData, coarseTSData, '.', 'color', color, 'markersize', 5)
            scatter(finTSData,coarseTSData, color,'filled','SizeData',3)
            alpha(0.5)
            hold on

        end
        xlabel({plotScoreNames{iscore}; 'fineTS ICs'}, 'fontsize', 8)
        ylabel({plotScoreNames{iscore}; 'coarseTS ICs'}, 'fontsize', 8)
        
        if iscore == 4; ylabel({plotScoreNames{3}; 'coarseTS ICs'}, 'fontsize', 8); end
        if iscore == 1; title(currPeriod, 'fontsize', 8); end


        axis square

    end

    linkaxes([ax.PRE ax.POST], 'xy')
end

%% adding a third dimension based on Bayesian scores to the scatterplot in the last section 


plotScoreNames = {'summed strength'; 'mean strength'; 'max strength'; 'no bins w/ sig. reactivation'};

BDscoreMethods = {'weightedCorr'; 'replayScore'};
BDshuffles     = {'unitIDshuffle'; 'wPBEtimeswap'};

for iBDscore = 1
for iBDshuffle = 1

figure;
for iscore = 1:4
    currScore = scoreNames{iscore};
    
    for iperiod = 1:3
        
        currPeriod = periodNames{iperiod};
        ax.(currPeriod) = subplot(4,3, (iscore-1)*3 + iperiod);

        itype = 1;
        currDataType = dataTypeNames{itype};

        finTSData    = activationStrPBEs_scores.(scoreNames{iscore}).fineTS.(currPeriod).(currDataType);

        if iscore == 4
           coarseTSData = activationStrPBEs_scores.(scoreNames{3}).coarseTS.(currPeriod).(currDataType);
        else
           coarseTSData = activationStrPBEs_scores.(scoreNames{iscore}).coarseTS.(currPeriod).(currDataType);
        end
        

        colormap('jet')
        BDscores = BDseqscore.(currPeriod).data.(BDshuffles{iBDshuffle}).(BDscoreMethods{iBDscore}).prctilescore;

        color = BDscores;
        scatter(finTSData,coarseTSData, 3, color,'filled')
        alpha(0.5)
        

        xlabel({plotScoreNames{iscore}; 'fineTS ICs'}, 'fontsize', 8)
        ylabel({plotScoreNames{iscore}; 'coarseTS ICs'}, 'fontsize', 8)
        
        if iscore == 4; ylabel({plotScoreNames{3}; 'coarseTS ICs'}, 'fontsize', 8); end
        if iscore == 1; title(currPeriod, 'fontsize', 8); end

        if iscore == 1 && iperiod == 1; title({BDscoreMethods{iBDscore};'';BDshuffles{iBDshuffle};'';currPeriod}, 'fontsize', 8); end


        axis square

    end

    linkaxes([ax.PRE ax.POST], 'xy')
end

end
end




%% scatterplot - CoarseTS ICs' recativation scores vs Bayesian replay scores
% total spikes fired during a PBE can be plotted as a third dimension

totalSpikeCounts = struct('PRE', [], 'RUN', [], 'POST', []);
periodNames    = {'PRE', 'RUN', 'POST'};

for iperiod = 1:3
    currPeriod = periodNames{iperiod};
    
    currPBEs = binnedFiringPBEs.(currPeriod);
    nPBEs    = size(currPBEs, 1); 
    totalSpikeCounts.(currPeriod) = zeros(nPBEs, 1);
    
    for ipbe = 1:nPBEs
        currPBE = currPBEs{ipbe, 2};
        
        totalSpikeCounts.(currPeriod)(ipbe) = sum(currPBE(:));  
    end
end
        


plotScoreNames = {'summed strength'; 'mean strength'; 'max strength'; 'no bins w/ sig. reactivation'};

BDscoreMethods = {'weightedCorr'; 'replayScore'};
BDshuffles     = {'unitIDshuffle'; 'wPBEtimeswap'};

scoreNames     = {'summedStr', 'meanStr', 'maxStr', 'nActBins'};
timeScaleNames = {'coarseTS'; 'fineTS'};
periodNames    = {'PRE', 'RUN', 'POST'};


its         = 1;

iBDscore    = 2;
iBDshuffle  = 2;

iscore      = 1;

itype       = 1;

rhoPart = struct('PRE', [], 'RUN', [], 'POST', []);


figure;
for iperiod = 1:3
    currPeriod = periodNames{iperiod};
    currDataType = dataTypeNames{itype};
    
    BDscores  = BDseqscore.(currPeriod).data.(BDshuffles{iBDshuffle}).(BDscoreMethods{iBDscore}).zscore;
    ICAscores = activationStrPBEs_scores.(scoreNames{iscore}).(timeScaleNames{its}).(currPeriod).(currDataType);

    
    
    color = totalSpikeCounts.(currPeriod);
    
    rhoPart.(currPeriod) = partialcorr(BDscores,ICAscores, color);
    
    
    BDscore_IQR  = prctile(BDscores, 75) - prctile(BDscores, 25);
    ICAscore_IQR = prctile(ICAscores, 75) - prctile(ICAscores, 25);
    
    
    sortedBDscores = sort(BDscores, 'ascend');
    minBDscore = sortedBDscores(find(sortedBDscores >= prctile(BDscores, 25)-1.5*BDscore_IQR, 1, 'first'));
    maxBDscore = sortedBDscores(find(sortedBDscores <= prctile(BDscores, 75)+1.5*BDscore_IQR, 1, 'last'));
    
    
    
    BDscoreBins = linspace(minBDscore,maxBDscore, 10);
    scoreBinCenters  = BDscoreBins(1:end-1) + diff(BDscoreBins(1:2))/2;
    BDscoreBins(end) = BDscoreBins(end) + 0.01;
    
    
    nBDscoreBins = numel(BDscoreBins)-1;
    [counts, binIdx] = histc(BDscores, BDscoreBins);
    counts(end) = [];
    
    
    actStrengthStats = zeros(nBDscoreBins, 3); % median, first quartile, 3rd quartile, resp.
    
    for ibin = 1:nBDscoreBins
        
        currICscores = ICAscores(binIdx == ibin);
        
        actStrengthStats(ibin, 1) = median(currICscores);
        actStrengthStats(ibin, 2) = prctile(currICscores, 25);
        actStrengthStats(ibin, 3) = prctile(currICscores, 75);
         
    end
 
    
    
    outliers = (BDscores-prctile(BDscores, 75)) > 2 * BDscore_IQR | ...
                (BDscores-prctile(BDscores, 25)) < -2 * BDscore_IQR | ...
                (ICAscores-prctile(ICAscores, 75)) > 2 * ICAscore_IQR | ...
                (ICAscores-prctile(ICAscores, 25)) < -2 * ICAscore_IQR;
    
    
    ax.(currPeriod) = subplot(4,3, [iperiod+3 iperiod+6]);
    
%     scatter(BDscores(~outliers), ICAscores(~outliers), 3, color(~outliers), 'filled')
    scatter(BDscores, ICAscores, 3, color, 'filled')
    
    p = polyfit(BDscores(~isnan(BDscores)), ICAscores(~isnan(BDscores)), 1);
    fitline = p(1)*BDscores + p(2);
    hold on
    plot(BDscores(~isnan(BDscores)), fitline(~isnan(BDscores)), 'k', 'linewidth', 2)
    
    [rho, pval] = corr(BDscores(~isnan(BDscores)), ICAscores(~isnan(BDscores)));
    text(BDscores(end)-1, fitline(end)+200, sprintf('corr=%.2f\npvalue=%.3f', rho, pval), 'fontsize', 8, 'color', [0.5 0.5 0.5])
    
    
    colormap('jet')
    if iperiod == 1; ylabel({plotScoreNames{iscore}; [timeScaleNames{its} 'ICs']}, 'fontsize', 10); end
    xlabel('BD score', 'fontsize', 10)
    title(currPeriod, 'fontsize', 10)
%     ylim([0 1000])
%     hold on
    
    ax2.(currPeriod) = subplot(4,3, iperiod + 9);
    
    hold on
    plot(scoreBinCenters, actStrengthStats(:, 1), 'r', 'linewidth', 3)
    plot(scoreBinCenters, actStrengthStats(:, 2), 'r', 'linewidth', 1)
    plot(scoreBinCenters, actStrengthStats(:, 3), 'r', 'linewidth', 1)
    hold off
    
    if iperiod == 1; ylabel({plotScoreNames{iscore}; 'coarseTS ICs'}, 'fontsize', 10); end
    xlabel('BD score', 'fontsize', 10)
 
end



linkaxes([ax.PRE ax.RUN ax.POST], 'xy')

linkaxes([ax2.PRE ax2.POST], 'y')
linkaxes([ax2.PRE ax2.RUN ax2.POST], 'x')

linkaxes([ax.PRE ax2.PRE], 'x')



%% BD scores and IC score as a function of PBEs total spike counts


its         = 1;

iBDscore    = 1;
iBDshuffle  = 1;

iscore      = 1;

itype       = 1;

figure; 

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    currDataType = dataTypeNames{itype};

    
    BDscores  = BDseqscore.(currPeriod).data.(BDshuffles{iBDshuffle}).(BDscoreMethods{iBDscore}).zscore;
    ICAscores = activationStrPBEs_scores.(scoreNames{iscore}).(timeScaleNames{its}).(currPeriod).(currDataType);

    firingRates = totalSpikeCounts.(currPeriod);
    
    
    
    
    %%%
    ax0.(currPeriod) = subplot(5,3,iperiod+3);
    scatter(BDscores, ICAscores, 3, firingRates, 'filled')
    
    p = polyfit(BDscores(~isnan(BDscores)), ICAscores(~isnan(BDscores)), 1);
    fitline = p(1)*BDscores + p(2);
    hold on
    plot(BDscores(~isnan(BDscores)), fitline(~isnan(BDscores)), 'k', 'linewidth', 2)
    
    [rho, pval] = corr(BDscores(~isnan(BDscores)), ICAscores(~isnan(BDscores)));
    [rhoPart, pvalPart] = partialcorr(BDscores,ICAscores, firingRates);
    
    text(BDscores(end)-1, fitline(end)+200, sprintf('corr=%.2f\npvalue=%.3f', rho, pval), 'fontsize', 8, 'color', [0.5 0.5 0.5])
    
    
    
    colormap('jet')
    
    title(currPeriod, 'fontsize', 10)
    if iperiod == 1; ylabel({plotScoreNames{iscore}; [timeScaleNames{its} 'ICs']}, 'fontsize', 10); end
    xlabel('BD score', 'fontsize', 10)
    
    %%%
    ax1.(currPeriod) = subplot(5,3,iperiod+6);
    scatter(BDscores, firingRates, 3, 'k', 'filled')
    
    p = polyfit(BDscores(~isnan(BDscores)), firingRates(~isnan(BDscores)), 1);
    fitline = p(1)*BDscores + p(2);
    hold on
    plot(BDscores(~isnan(BDscores)), fitline(~isnan(BDscores)), 'b', 'linewidth', 2)
    
    [rho, pval] = corr(BDscores(~isnan(BDscores)), firingRates(~isnan(BDscores)));
    text(BDscores(end)-1, fitline(end)+100, sprintf('corr=%.2f\npvalue=%.3f', rho, pval), 'fontsize', 8, 'color', [0.5 0.5 0.5])
    
    xlabel('BD score', 'fontsize', 10)
    if iperiod == 1; ylabel('PBE spike count', 'fontsize', 10); end

    
    
    %%%
    ax2.(currPeriod) = subplot(5,3,iperiod+9);
    
    scatter(ICAscores, firingRates, 3, 'k', 'filled')
    
    p = polyfit(ICAscores(~isnan(ICAscores)), firingRates(~isnan(ICAscores)), 1);
    fitline = p(1)*ICAscores + p(2);
    hold on
    plot(ICAscores(~isnan(ICAscores)), fitline(~isnan(ICAscores)), 'b', 'linewidth', 2)
    
    [rho, pval] = corr(ICAscores(~isnan(ICAscores)), firingRates(~isnan(ICAscores)));
    text(ICAscores(end)-1, fitline(end)+100, sprintf('corr=%.2f\npvalue=%.3f', rho, pval), 'fontsize', 8, 'color', [0.5 0.5 0.5])
    
    xlabel({plotScoreNames{iscore}; 'coarseTS ICs'}, 'fontsize', 10)
    if iperiod == 1; ylabel('PBE spike count', 'fontsize', 10); end
    
    
    aa = subplot(5,3, iperiod+12);
    
    text(0,0, sprintf('partialcorr=%.2f\npvalue=%.3f', rhoPart, pvalPart), 'fontsize', 12, 'color', [0.5 0.5 0.5])
    set(aa, 'Visible', 'off');
    
end

linkaxes([ax0.PRE ax0.RUN ax0.POST], 'xy')
linkaxes([ax1.PRE ax1.RUN ax1.POST], 'xy')
linkaxes([ax2.PRE ax2.RUN ax2.POST], 'xy')

%% Correlation between IC scores and Bayesian scores 
% Correlations separately can be calclulated for each time scale


BDmethodNames  = {'WC-UI'; 'WC-TS'; 'RT-UI'; 'RT-TS'}; 
BDmethodColors = {'r';[55 55 0]/255; 'g';'b'};

figure;

for iscore = 1:4

    currScore = scoreNames{iscore};

    for iperiod = 1:3

        currPeriod = periodNames{iperiod};
        ax.(currPeriod) = subplot(5,3, (iscore)*3 + iperiod);    

        
% 
%         if iscore == 4
%            coarseTSData = activationStrPBEs_scores.(scoreNames{3}).coarseTS.(currPeriod).(currDataType);
%         else
%            coarseTSData = activationStrPBEs_scores.(scoreNames{iscore}).coarseTS.(currPeriod).(currDataType);
%         end

        for iBDscore = 1:2
            for iBDshuffle = 1:2
                
                BDscores = BDseqscore.(currPeriod).data.(BDshuffles{iBDshuffle}).(BDscoreMethods{iBDscore}).zscore;
                for itype = 1:2

                    currDataType = dataTypeNames{itype};
                    ICAData = activationStrPBEs_scores.(currScore).fineTS.(currPeriod).(currDataType);

                    corrVal = corr(ICAData(~isnan(BDscores) & ~isnan(ICAData)), BDscores(~isnan(BDscores) & ~isnan(ICAData)));
                    
                    if itype == 1; markerColor = BDmethodColors{(iBDscore-1)*2+iBDshuffle}; else; markerColor = [0.7 0.7 0.7]; end
                    plot((iBDscore-1)*2+iBDshuffle, corrVal, 'marker', '.', 'color', markerColor, 'markersize', 20)
                    hold on
                end

            end
        end


        xlim([0 5])
        xticks(1:4)
        xticklabels(BDmethodNames)
        xtickangle(45)

        if iperiod == 1; ylabel({'correlation w/ fineTS ICs'; currScore}, 'fontsize', 8); end
        if iscore  == 4; xlabel({'BD methods'}); end
        if iscore == 1; title(currPeriod, 'fontsize', 10);end
        
        set(gca, 'box', 'off', 'linewidth', 1)
        
    end
    linkaxes([ax.PRE ax.RUN ax.POST], 'y')
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
%         ICtunings.(currTimeScale)(:, ii) = sum(activationStrength.(currTimeScale).RUN.original(:, posbinIdx.(currTimeScale)(:, 1) == uniqPositions(ii)), 2);  
        ICtunings.(currTimeScale)(:, ii) = sum(activationStrength_tuning.(currTimeScale)(:, posbinIdx.(currTimeScale)(:, 1) == uniqPositions(ii)), 2)/numel(find(posbinIdx.(currTimeScale)(:, 1) == uniqPositions(ii)))/binDur.coarseTS; 
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

currDataType = 'actual';


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
                xcorr_actual_PBE(pbe, :) = xcorr(activationStrPBEs.fineTS.(currPeriod).(currDataType){pbe}(ii, :), activationStrPBEs.fineTS.(currPeriod).(currDataType){pbe}(jj, :), maxLag);
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

                s1 = activationStrPBEs.fineTS.(currPeriod).(currDataType){pbe}(ii, :);
                s2 = activationStrPBEs.fineTS.(currPeriod).(currDataType){pbe}(jj, :);

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

for its = 1
    
    currTimeScale = timeScaleNames{its};
    
    nICs          = size(ICs.(currTimeScale), 1);

    
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

    imagesc(currICweights(:, runTemplate_LR)', [-max(abs(currICweights(:))) max(abs(currICweights(:)))]); colormap(cm); colorbar; 
    set(gca, 'ydir', 'normal', 'fontsize', 10, 'linewidth', 1)
    xlabel('principal components', 'fontsize', 12)
    ylabel('units', 'fontsize', 12)


    % sptial tuning of the ICs
    subplot(2,5,8:10)

    ICtunings_norm = zscore(ICtunings.(currTimeScale), [], 2); % excluding the zero position bin and z scoring
    
    imagesc(1:size(ICtunings_norm, 2)*2, 1:nICs,  ICtunings_norm, [-max(ICtunings_norm(:)) max(ICtunings_norm(:))]); colormap(cm); colorbar;
    set(gca, 'ydir', 'normal', 'fontsize', 12)
    xlabel('track position(cm)', 'fontsize', 12)
    
    
    subplot(2,5, 6)
    
    plot(max(ICtunings.(currTimeScale), [], 2), (1:nICs)- 0.5, 'linewidth', 3)
    ylim([0 20])
    set(gca, 'fontsize', 12)

    xlabel('peak activation rate (1/sec)', 'fontsize', 12)
    ylabel('principal components', 'fontsize', 12)
end


%% activations of PCs lap by lap

figure;
ics2plot = 1:20;


actThresh = 0.0001;

for ic = 1: length(ics2plot)
    
    currIC = ics2plot(ic);
    
    subplot(length(ics2plot),1,ic)
    
    currActivations = activationStrength_tuning.coarseTS(currIC, :);
    currActivations(currActivations < actThresh) = nan;
    
    plot(fileinfo.xyt2(:,1), fileinfo.xyt2(:,2) - behavior.time(2,1), '.', 'markersize', 0.5, 'color', [0.5 0.5 0.5])
    
    hold on 
    
    s = scatter(linposcenters.coarseTS(posbinIdx.coarseTS(:,1)), binCenters4Tuning.coarseTS - behavior.time(2,1), currActivations);
    
    s.MarkerEdgeColor = 'none';
    s.MarkerFaceColor = 'b';
    s.MarkerFaceAlpha = 0.5;
    
    set(gca, 'linewidth', 1)
    
    if ic ~= length(ics2plot)
        set(gca, 'xticklabel', [])
    end
    
    if ic == length(ics2plot)
        xlabel('track position(cm)', 'fontsize', 12)
    end
    
    ylabel({['PC ' num2str(currIC)]; ''; 'time(sec)'}, 'fontsize', 12)
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
    
    nor = 5; % change these if the number of ICs is higher than 32
    noc = 4;
    
    plotICactivationStrengths(activationStrength.(currTimeScale), binCenters.(currTimeScale), behavior, nor, noc, currTimeScale)
end


%% Totoal reactivation strengths by summing over ICs (PRE & POST)

for its = 1
    
    currTimeScale = timeScaleNames{its};
    
    nor = 1;
    noc = 1;

    plotICactivationStrengths(overallActivationStrength.(currTimeScale), binCenters.(currTimeScale), behavior, 1, 1, [currTimeScale '_overall'])
end


%% test

spikeStruct.RLtemplateIdx = zeros(size(spikeStruct.unit));
spikeStruct.LRtemplateIdx = zeros(size(spikeStruct.unit));


for ii = 1: numel(runTemplate_RL)
    spikeStruct.RLtemplateIdx(spikeStruct.unit == runTemplate_RL(ii)) = ii;
end


for ii = 1: numel(runTemplate_LR)
    spikeStruct.LRtemplateIdx(spikeStruct.unit == runTemplate_LR(ii)) = numel(runTemplate_LR) - ii + 1;
end


%%

PBEtimes = struct('PRE', secondaryPBEs_pre, 'RUN', secondaryPBEs_run, 'POST', secondaryPBEs_post);


figure; 

duration2Plot = 300; % in sec
noBins2Plot   = (duration2Plot/binDur.coarseTS - binOverlapRatio.coarseTS) / (1 - binOverlapRatio.coarseTS);
startTimes    = [400 50 500];  % in respect to beginning of the period
% startBins     = (startTimes/binDur.coarseTS - binOverlapRatio.coarseTS)./ (1 - binOverlapRatio.coarseTS)+10; % because of the overlap

currStepSize  = (1-binOverlapRatio.coarseTS) * binDur.coarseTS;
startBins     = startTimes/currStepSize + 1;
startBins     = floor(startBins);


start20msBins   = startTimes/binDur.fineTS + 1;
start20msBins   = floor(start20msBins);
no20msBins2Plot = duration2Plot/binDur.fineTS;


%%% Bayesian posterior probability matrices for 20ms time-binned firings


nPosBins = size(spatialTunings_LR, 2);

spatialTunings_LR(spatialTunings_LR==0) = 1e-4;
spatialTunings_RL(spatialTunings_RL==0) = 1e-4;

postPr = struct('PRE', [], 'RUN', [], 'POST', []);


for iperiod = 1:3
    iperiod
    currPeriod           = periodNames{iperiod};
    currFirings          = binnedFiring_nonOverlap.fineTS.(currPeriod){1,2}(:, start20msBins(iperiod)+(0:no20msBins2Plot-1));
    
    sumFiring = sum( currFirings, 1);
    idx       = find(sumFiring == 0);
    
    
    postPr.(currPeriod)  = struct('LR', [], 'RL', [], 'integrated', []);
    
    postPr.(currPeriod).LR = baysDecoder(currFirings, spatialTunings_LR, binDur.fineTS);
    postPr.(currPeriod).RL = baysDecoder(currFirings, spatialTunings_RL, binDur.fineTS);
    
    postPr.(currPeriod).integrated = (postPr.(currPeriod).LR + postPr.(currPeriod).RL) ./ repmat(sum(postPr.(currPeriod).LR + postPr.(currPeriod).RL), [nPosBins, 1]);
    postPr.(currPeriod).integrated(:, idx) = 0;
    
    
    postPr.(currPeriod).LR = postPr.(currPeriod).LR ./ repmat(sum(postPr.(currPeriod).LR, 1), [size(postPr.(currPeriod).LR, 1) 1]);
    postPr.(currPeriod).RL = postPr.(currPeriod).RL ./ repmat(sum(postPr.(currPeriod).RL, 1), [size(postPr.(currPeriod).RL, 1) 1]);
    
    postPr.(currPeriod).LR(:, idx) = 0;
    postPr.(currPeriod).RL(:, idx) = 0;
    
end




popBurstTimeCourse = struct('PRE', sdat_pre, 'RUN', sdat_run, 'POST', sdat_post);

periodNames = {'PRE'; 'RUN'; 'POST'};

for iperiod = 1:3
    
    currPeriod = periodNames{iperiod};
    ax0.(currPeriod) = subplot(12, 3, ((1:6)-1)*3 +iperiod);    
    
    currActivationStrs = activationStrength.coarseTS.(currPeriod).original;
    
    currStartTime = startTimes(iperiod)+behavior.time(iperiod, 1);
    currEndTime   = currStartTime + duration2Plot;
    
    spikes2Plot = find(spikeStruct.t >= currStartTime & spikeStruct.t <= currEndTime & ismember(spikeStruct.qclu, qclus));
    

    currBinCenters = (startTimes(iperiod):currStepSize:(startTimes(iperiod) + duration2Plot -binDur.coarseTS)) + binDur.coarseTS/2;
    
    
    
    plot(currBinCenters-startTimes(iperiod), currActivationStrs(:, floor(startBins(iperiod)+(0:noBins2Plot-1)))', 'linewidth',1.5)
    hold on
    plot(spikeStruct.t(spikes2Plot)-startTimes(iperiod)-behavior.time(iperiod, 1), spikeStruct.LRtemplateIdx(spikes2Plot)+ 50, '.', 'color', [.5 .5 .5 0.5], 'markersize', 7)
    hold on
    plot(spikeStruct.t(spikes2Plot)-startTimes(iperiod)-behavior.time(iperiod, 1), spikeStruct.RLtemplateIdx(spikes2Plot)+ 200, '.', 'color', [.5 .5 .5 0.5], 'markersize', 7)
    
%     scatter(spikeStruct.t(spikes2Plot)-startTimes(iperiod)-behavior.time(iperiod, 1), spikeStruct.LRtemplateIdx(spikes2Plot)+ 50, 5, ICweights.coarseTS(4, spikeStruct.unit(spikes2Plot)), 'filled')
%     hold on
%     scatter(spikeStruct.t(spikes2Plot)-startTimes(iperiod)-behavior.time(iperiod, 1), spikeStruct.RLtemplateIdx(spikes2Plot)+ 200, 5, ICweights.coarseTS(4, spikeStruct.unit(spikes2Plot)), 'filled')
%     
    
    colormap(ax0.(currPeriod), cm)
    ylim([-20 350])
    
    
    if iperiod == 3
       line([0 5], [250 250], 'linewidth', 2, 'color', 'k')
%        line([0 0], [400 500], 'linewidth', 2, 'color', 'k')
    end
    
    title(periodNames{iperiod}, 'fontsize', 12, 'FontWeight', 'normal')
    
    if iperiod == 1; ylabel('ICs reactivation strength', 'fontsize', 12); end
    set(ax0.(currPeriod), 'xcolor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)

    
    
    % ripple LPF 
    ax1.(currPeriod) = subplot(12, 3, (7-1)*3 +iperiod); 
    
    currRipple = rippleLFP(currStartTime*1250+1: currEndTime*1250);
    plot([1:length(currRipple)]/1250, currRipple, 'linewidth', 1.5, 'color', [.5 .5 .5]) %  
    
    ylim([-1000 1000])
    
    if iperiod== 1; yh = ylabel('ripple band', 'fontsize', 8); set(yh, 'rotation', 0);  end
    
    set(ax1.(currPeriod), 'XColor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)

    
    % population burst
    
    ax2.(currPeriod) = subplot(12, 3, (8-1)*3 +iperiod); 
    
    currPopBurst = popBurstTimeCourse.(currPeriod)(startTimes(iperiod)*1000+1: (startTimes(iperiod)+duration2Plot)*1000);
    plot([1:length(currPopBurst)]/1000, currPopBurst, 'linewidth', 1.5, 'color', [.5 .5 .5]) %  
    hold on
    line([1 length(currPopBurst)]/1000, [3 3], 'linewidth', 1, 'color', 'r', 'linestyle', ':')
    
    ylim([-1 8])
    
    if iperiod == 1; yh= ylabel('population firing rate(z)', 'fontsize', 8); set(yh, 'rotation', 0); end

    set(ax2.(currPeriod), 'XColor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)

    
%     % BD replay scores
%     ax3.(currPeriod) = subplot(12, 3, ((9:10)-1)*3 + iperiod);
%     
%     currBayesianScores  = BDseqscore.(currPeriod).data.unitIDshuffle.replayScore.prctilescore;
%     BDscoreStepFunction = nan(floor(diff(behavior.time(iperiod, :)))*1000,1);
%     
%     for ipbe = 1:size(PBEtimes.(currPeriod), 1)
%         BDscoreStepFunction(floor((PBEtimes.(currPeriod)(ipbe, 1)*1000 : PBEtimes.(currPeriod)(ipbe, 2)*1000)-behavior.time(iperiod, 1)*1000)) = currBayesianScores(ipbe);
%     end
%     
% 
%     currBDscoreStepFunction = BDscoreStepFunction(startTimes(iperiod)*1000+1:(startTimes(iperiod)+duration2Plot)*1000);
%     plot([1:length(currBDscoreStepFunction)]/1000, currBDscoreStepFunction, 'linewidth', 1.5, 'color', [35 110 150]/255) %  
%     
%     hold on
%     
%     currBayesianScores  = BDseqscore.(currPeriod).data.unitIDshuffle.weightedCorr.prctilescore;
%     BDscoreStepFunction = nan(floor(diff(behavior.time(iperiod, :)))*1000,1);
%     
%     for ipbe = 1:size(PBEtimes.(currPeriod), 1)
%         BDscoreStepFunction(floor((PBEtimes.(currPeriod)(ipbe, 1)*1000 : PBEtimes.(currPeriod)(ipbe, 2)*1000)-behavior.time(iperiod, 1)*1000)) = currBayesianScores(ipbe);
%     end
%     
%     currBDscoreStepFunction = BDscoreStepFunction(startTimes(iperiod)*1000+1:(startTimes(iperiod)+duration2Plot)*1000);
%     plot([1:length(currBDscoreStepFunction)]/1000, currBDscoreStepFunction, 'linewidth', 1.5, 'color', [21 178 211]/255) %  
%     

%     
%     hold on
%     
%     currBayesianScores  = BDseqscore.(currPeriod).data.wPBEtimeswap.replayScore.prctilescore;
%     BDscoreStepFunction = nan(floor(diff(behavior.time(iperiod, :)))*1000,1);
%     
%     for ipbe = 1:size(PBEtimes.(currPeriod), 1)
%         BDscoreStepFunction(floor((PBEtimes.(currPeriod)(ipbe, 1)*1000 : PBEtimes.(currPeriod)(ipbe, 2)*1000)-behavior.time(iperiod, 1)*1000)) = currBayesianScores(ipbe);
%     end
%     
%     currBDscoreStepFunction = BDscoreStepFunction(startTimes(iperiod)*1000+1:(startTimes(iperiod)+duration2Plot)*1000);
%     plot([1:length(currBDscoreStepFunction)]/1000, currBDscoreStepFunction, 'linewidth', 1.5, 'color', [247 95 28]/255) %  
%     
%     hold on
%  
%     
%     currBayesianScores  = BDseqscore.(currPeriod).data.wPBEtimeswap.weightedCorr.prctilescore;
%     BDscoreStepFunction = nan(floor(diff(behavior.time(iperiod, :)))*1000,1);
%     
%     for ipbe = 1:size(PBEtimes.(currPeriod), 1)
%         BDscoreStepFunction(floor((PBEtimes.(currPeriod)(ipbe, 1)*1000 : PBEtimes.(currPeriod)(ipbe, 2)*1000)-behavior.time(iperiod, 1)*1000)) = currBayesianScores(ipbe);
%     end
%     
%     currBDscoreStepFunction = BDscoreStepFunction(startTimes(iperiod)*1000+1:(startTimes(iperiod)+duration2Plot)*1000);
%     plot([1:length(currBDscoreStepFunction)]/1000, currBDscoreStepFunction, 'linewidth', 1.5, 'color', [255 154 0]/255) %  
%     
%     
%     if iperiod == 3
%         legend('UI-RT', 'UI-WC', 'TS-RT', 'TS-WC', 'box', 'off', 'linewidth', 2)
%     end
%     
%     
%     if iperiod == 1; yh= ylabel('BD score', 'fontsize', 8); set(yh, 'rotation', 0); end
%     set(ax3.(currPeriod), 'box', 'off', 'TickDir', 'out', 'linewidth', 2) % , 'XColor', 'none'

    
    % BD posterior probabilities
    
    ax3.(currPeriod) = subplot(12, 3, ((9:10)-1)*3 + iperiod);
    
    currPosterior = postPr.(currPeriod).integrated;
    
    currBinCenters = (startTimes(iperiod):binDur.fineTS:(startTimes(iperiod) + duration2Plot -binDur.fineTS)) + binDur.fineTS/2;
    
    cl(1) = min(currPosterior(:));
    cl(2) = min(currPosterior(:)) + 0.1*(max(currPosterior(:))-min(currPosterior(:)));
    
    imagesc(currBinCenters-startTimes(iperiod), 1:size(currPosterior, 1), currPosterior, cl); colormap(ax3.(currPeriod), 'hot')
    
end
    

linkaxes([ax0.PRE ax0.RUN ax0.POST], 'y')
linkaxes([ax1.PRE ax1.RUN ax1.POST], 'y')
linkaxes([ax2.PRE ax2.RUN ax2.POST], 'y')
linkaxes([ax3.PRE ax3.RUN ax3.POST], 'y')


% linkaxes([ax0.PRE ax1.PRE ax2.PRE ax3.PRE ax0.RUN ax1.RUN ax2.RUN ax3.RUN ax0.POST ax1.POST ax2.POST ax3.POST], 'x')

linkaxes([ax0.PRE ax1.PRE ax2.PRE ax3.PRE], 'x')
linkaxes([ax0.RUN ax1.RUN ax2.RUN ax3.RUN], 'x')
linkaxes([ax0.POST ax1.POST ax2.POST ax3.POST], 'x')


%% visualizing spike raster, reactivation strengths and Bayesian decoding in population bursts around high reactivation strength points

iperiod = 3;
currPeriod = periodNames{iperiod};

% % points of high reactivation (we have already calculated this)

% population bursts around the high reactivation points


selectICs = 1:5;
hiActEvts = cell(numel(selectICs), 1);

for ic = 1:numel(selectICs)
    
    temp = findHiPCactPeriods(ICpeakReactivationPoints, popBurstTimeCourse, selectICs(ic));
    hiActEvts{ic} = temp.(currPeriod);

end
%%
% measuring replay scores (4 different methods)








% plotting everything together

plotHiPCactPeriods(hiActEvts, activationStrength, binCenters, selectICs, ICweights, spikeStruct, binnedFiring_nonOverlap, popBurstTimeCourse, rippleLFP, binDur, spatialTunings_LR, spatialTunings_RL, behavior, currPeriod)








%% visualizing spike rasters, posterior probability and reactivation strengh of different ICs within individual PBEs

sessionName(sessionName == '_') = '-';

for iperiod = 3
    currentPeriod = periodNames{iperiod};
    
    directory = fullfile(mainDir, 'ICA', currentPeriod);
    if ~exist(directory, 'dir')
        mkdir(directory)
    end

    currBDscores  = BDseqscore.(currentPeriod).data.unitIDshuffle.weightedCorr.prctilescore;
    currICAscores = activationStrPBEs_scores.summedStr.coarseTS.(currentPeriod).actual;
    
    % PBEs with high BD scores
%     PBE2Plot = find(currBDscores > 95);
%     PBE2Plot = PBE2Plot(randperm(numel(PBE2Plot), 40));
      
    [~, sortInd] = sort(currBDscores, 'descend');
    PBE2Plot = sortInd(1:40);

%     [~, sortInd] = sort(currICAscores, 'descend');
%     PBE2Plot = sortInd(1:40);
    
    
    filename  = ['ICsActivationDuringPBEs_' currentPeriod 'highBDscores'];
    %%First 40 PBEs with highest IC reactivation strengths are shown
    textOnTop = sprintf('Session: %s\nPeriod: %s\nFirst 40 PBEs with highest BD replay scores are shown.\n\\color{magenta}UI-WC: BD method based on Unit ID shuffle and weighted correlation\n\\color{blue}ICscore: maximum activation of coarse time scale ICs during the PBE\n', sessionName, currentPeriod);
    
    plotICactivation4PBEs(binnedFiringPBEs.(currentPeriod), activationStrPBEs.coarseTS.(currentPeriod).actual, posteriorProbMatrix.(currentPeriod).data, currBDscores, currICAscores, runTemplate_LR, runTemplate_RL, PBE2Plot, directory, filename, textOnTop)

% 
%     % randomly selected PBEs
%     PBE2Plot = randperm(numel(currBDscores), 40);
%     
%     filename  = ['ICsActivationDuringPBEs_' currentPeriod 'randomPBEs'];
%     textOnTop = sprintf('Session: %s\nPeriod: %s', sessionName, currentPeriod);
%     
%     plotICactivation4PBEs(binnedFiringPBEs.(currentPeriod), activationStrPBEs.fineTS.(currentPeriod).actual, posteriorProbMatrix.(currentPeriod).data, currBDscores, runTemplate_LR, runTemplate_RL, PBE2Plot, directory, filename, textOnTop)

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

colorLimit = [-0.3 0.3];

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
xlim([-1 1])
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
xlim([-1 1])
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
%% functions

%%
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

%%
function smoothedMat = arrayConv(inputMat, avgWindow)

smoothedMat = zeros(size(inputMat));
for ii = 1: size(inputMat, 1)
    smoothedMat(ii, :) = conv(inputMat(ii, :), avgWindow, 'same');
    
    % dealing with the edge effect by replacing the first and last segements of the data (segments equal in length with half length of  avgWindow) with Nans 
    avgWindowlen = length(avgWindow);
    
    smoothedMat(ii, [1:floor(avgWindowlen/2) (end-floor(avgWindowlen/2)):end]) = nan;   
    
end

end

%%
function z = zscore_KM(x, mu, sigma)

z = (x - repmat(mu, [1 size(x, 2)]))./repmat(sigma, [1 size(x, 2)]);

end


%%
function [hiActEvts, hiActBins] = findHiPCactPeriods(ICpeakReactivationPoints, popBurstTimeCourse, selectICs)


periodNames = {'PRE'; 'RUN'; 'POST'};
hiActEvts   = struct('PRE', [], 'RUN', [], 'POST', []);
hiActBins   = struct('PRE', [], 'RUN', [], 'POST', []);

addedFlankPeriod = 500; % in miliseconds, in case the

for iperiod = 1:3
    
    currPeriod          = periodNames{iperiod}; 
    
    currMUATimeCourse   = popBurstTimeCourse.(currPeriod);
    currMUATimePntIdx   = (1:numel(currMUATimeCourse))';
    
    currICpeakReactPnts = cell2mat(ICpeakReactivationPoints(selectICs, iperiod)');
%     curr250msBinCenters = binCenters.coarseTS.(currPeriod);
    
    
    nPeaks         = numel(currICpeakReactPnts);
    
    hiActEvts.(currPeriod) = zeros(nPeaks, 3);
    hiActBins.(currPeriod) = zeros(nPeaks, 4); % columns 3 and 4 corresponding to peak population bursts within 250ms time bin and the whole population burst period
    
    for ipeak = 1:nPeaks
        
        currPeak = currICpeakReactPnts(ipeak);
        
        
        if currMUATimeCourse(currPeak) ~= 0
            
            if currMUATimeCourse(currPeak) > 0
                
                searchStartPnt = max([currPeak-1000 1]);
                searchEndPnt   = min([currPeak+1000 numel(currMUATimeCourse)]);
                
                aroundPeakMUATimeCourse = currMUATimeCourse(searchStartPnt : searchEndPnt);
                arondPeakMUATimePntIdx  = currMUATimePntIdx(searchStartPnt : searchEndPnt);

                preZeroCrossing  = arondPeakMUATimePntIdx(find(aroundPeakMUATimeCourse < 0 & arondPeakMUATimePntIdx < currPeak, 1, 'last'));
                if isempty(preZeroCrossing); preZeroCrossing = currPeak-addedFlankPeriod; end
                
                postZeroCrossing = arondPeakMUATimePntIdx(find(aroundPeakMUATimeCourse < 0 & arondPeakMUATimePntIdx > currPeak, 1, 'first'));
                if isempty(postZeroCrossing); postZeroCrossing = currPeak+addedFlankPeriod; end


                hiActEvts.(currPeriod)(ipeak, 1:2) = [preZeroCrossing postZeroCrossing];
                hiActEvts.(currPeriod)(ipeak, 3)   = max(currMUATimeCourse(preZeroCrossing:postZeroCrossing));
                
                hiActBins.(currPeriod)(ipeak, 1:2) = [currPeak-125 currPeak+125];
                hiActBins.(currPeriod)(ipeak, 3)   = max(currMUATimeCourse(currPeak-125:currPeak+125));
                hiActBins.(currPeriod)(ipeak, 4)   = max(currMUATimeCourse(preZeroCrossing:postZeroCrossing));
                
                
            elseif currMUATimeCourse(currPeak) < 0
                
                
                hiActEvts.(currPeriod)(ipeak, :) = nan;
                
%                 hiActEvts.(currPeriod)(ipeak, 1:2) = [currPeak-addedFlankPeriod currPeak+addedFlankPeriod];
%                 hiActEvts.(currPeriod)(ipeak, 3)   = max(currMUATimeCourse(currPeak-addedFlankPeriod:currPeak+addedFlankPeriod));
            end

        else
            
            hiActEvts.(currPeriod)(ipeak, :) = nan;
        end
        
    end
    
    % finding repeated detections
    temp  = [0; diff(hiActEvts.(currPeriod)(:,1))];
    
    hiActEvts.(currPeriod)(isnan(hiActEvts.(currPeriod)(:,1)) | temp == 0, :) = [];

end


end

%%
function plotHiPCactPeriods(hiActEvts, activationStrength, binCenters, selectICs, ICweights, spikeStruct, binnedFiring, popBurstTimeCourse, rippleLFP, binDur, spatialTunings_LR, spatialTunings_RL, behavior, currPeriod)


if strcmp(currPeriod, 'PRE')
    iperiod = 1;
elseif strcmp(currPeriod, 'RUN')
    iperiod = 2;
elseif strcmp(currPeriod, 'POST')
    iperiod = 3;
end


cm = redblue();
close all
% select a subset of high activation events to plot based on their
% corresponding peak MUAs

nICs = numel(selectICs);
setSize = 100;
selectHiActEvts = cell(nICs, 1);

for ic = 1%:nICs
    
    currIC = selectICs(ic);
    [~, sortInd] = sort(hiActEvts{currIC}(:, 3), 'descend'); % sorting based on max MUA 
    selectHiActEvts{ic} = hiActEvts{currIC}(sortInd(1:setSize), :);
%     selectHiActEvts{ic} = hiActEvts{currIC}(sortInd, :);
end

try
    selectHiActEvts = cell2mat(selectHiActEvts);
catch
    selectHiActEvts = cell2mat(selectHiActEvts');
end


% remove repeated events
[~, sortInd] = sort(selectHiActEvts(:,1), 'ascend');
selectHiActEvts = selectHiActEvts(sortInd, :);

temp  = [0; diff(selectHiActEvts(:,1))]; 
selectHiActEvts(isnan(selectHiActEvts(:,1)) | temp == 0, :) = [];




selectHiActEvts = selectHiActEvts(randperm(50),:);




for ii = 1:size(selectHiActEvts, 1)
    
    figure; 
    set(gcf, 'units', 'centimeters', 'position', [0 0 30 15.3]) % [0 0 53.5 15.3]
    flankingPeriod = 150; % in miliseconds
    
    currStartTime = (selectHiActEvts(ii, 1) - flankingPeriod)*1e-3 + behavior.time(iperiod, 1);
    currEndTime   = (selectHiActEvts(ii, 2) + flankingPeriod)*1e-3 + behavior.time(iperiod, 1);
    
    qclus = 1:3;
    spikes2Plot = find(spikeStruct.t >= currStartTime & spikeStruct.t <= currEndTime & ismember(spikeStruct.qclu, qclus));

    
    
    for ic = 1:nICs
        
        
    % ax1: spike rasters and activation strengths of PCs/ICs around the high PC activation strength periods
    
    ax1 = subplot(12, nICs, ((1:6)-1)*nICs +ic);  
    

    currActivationStrs = activationStrength.coarseTS.(currPeriod).original;
    
    
    binCentersCoarse = binCenters.coarseTS.(currPeriod); 
    
    currCoarseBinIdx = binCentersCoarse >= currStartTime & binCentersCoarse <= currEndTime;
    
    currBinCenters  = binCentersCoarse(currCoarseBinIdx);
    currActivations = currActivationStrs(selectICs, currCoarseBinIdx);
    % currBinCenters = (startTimes(iperiod):currStepSize:(startTimes(iperiod) + duration2Plot -binDur.coarseTS)) + binDur.coarseTS/2;

    
    plot(currBinCenters-currStartTime, currActivations(setdiff(selectICs, selectICs(ic)), :)', 'color', [0 0 0 0.15], 'linewidth',1.5)
    hold on
    plot(currBinCenters-currStartTime, currActivations(selectICs(ic), :)', 'color', 'k', 'linewidth',1.5)
    

    
    for ibin = 1:numel(currBinCenters)
        currPnt = currBinCenters(ibin) - currStartTime;
        
        if mod(ibin, 2) == 1
            
            line([currPnt-0.12 currPnt+0.12], [-20 -20], 'linewidth', 3, 'color', [0 0 0 0.2])
            
        else
            line([currPnt-0.12 currPnt+0.12], [-15 -15], 'linewidth', 3, 'color', [0 0 0 0.2])
        end
        
    end
    
    
    hold on
    
    scatter(spikeStruct.t(spikes2Plot)-currStartTime, spikeStruct.LRtemplateIdx(spikes2Plot)+ 50, 5, ICweights.coarseTS(selectICs(ic), spikeStruct.unit(spikes2Plot)), 'filled')
    colormap(ax1, cm)
    ylim([-20 200])
    caxis([-0.2 0.2])
    
    
    [maxreact, peakPnt] = max(currActivations(selectICs(ic), :));
    peakBinCenter = currBinCenters(peakPnt)-currStartTime;
    
    hold on
    line([flankingPeriod flankingPeriod]*1e-3, [-20 200], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([currEndTime-currStartTime-flankingPeriod*1e-3  currEndTime-currStartTime-flankingPeriod*1e-3], [-20 200], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    
    if maxreact > 20
    hold on
    line([peakBinCenter-0.125 peakBinCenter-0.125], [-20 200], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([peakBinCenter+0.125 peakBinCenter+0.125], [-20 200], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    end
    
    if ic == 1; ylabel('PCs reactivation strength', 'fontsize', 12); end
    set(ax1, 'xcolor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)    

    text(0.1, 200, sprintf('PC %d', selectICs(ic)), 'fontsize', 12)
    
    
    
    % ax2: ripple
    ax2 = subplot(12, nICs, (7-1)*nICs +ic); 
    
    currRipple = rippleLFP(currStartTime*1250+1: currEndTime*1250);
    plot([1:length(currRipple)]/1250, currRipple, 'linewidth', 1.5, 'color', [.5 .5 .5]) %  
    
    ylim([-1000 1000])
    
%     if ic== 1; yh = ylabel('ripple band', 'fontsize', 8); set(yh, 'rotation', 0);  end
    
    set(ax2, 'XColor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)

    
    hold on
    line([flankingPeriod flankingPeriod]*1e-3, [-1000 1000], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([currEndTime-currStartTime-flankingPeriod*1e-3  currEndTime-currStartTime-flankingPeriod*1e-3], [-1000 1000], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    
    if maxreact > 20
    hold on
    line([peakBinCenter-0.125 peakBinCenter-0.125], [-1000 1000], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([peakBinCenter+0.125 peakBinCenter+0.125], [-1000 1000], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    end
    
    

    % ax3: MUA
    ax3 = subplot(12, nICs, (8-1)*nICs +ic); 
    
    
    currPopBurst = popBurstTimeCourse.(currPeriod)((currStartTime- behavior.time(iperiod, 1))*1000+1: (currEndTime-behavior.time(iperiod, 1))*1000);
    plot([1:length(currPopBurst)]/1000, currPopBurst, 'linewidth', 1.5, 'color', [.5 .5 .5]) %  
    hold on
    line([1 length(currPopBurst)]/1000, [3 3], 'linewidth', 1, 'color', 'r', 'linestyle', ':')
    
    ylim([-2 10])
    
%     if ic == 1; yh= ylabel('population firing rate(z)', 'fontsize', 8); set(yh, 'rotation', 0); end

    set(ax3, 'XColor', 'none', 'box', 'off', 'TickDir', 'out', 'linewidth', 2)

    hold on
    line([flankingPeriod flankingPeriod]*1e-3, [-2 10], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([currEndTime-currStartTime-flankingPeriod*1e-3  currEndTime-currStartTime-flankingPeriod*1e-3], [-2 10], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    
    if maxreact > 20
    hold on
    line([peakBinCenter-0.125 peakBinCenter-0.125],  [-2 10], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([peakBinCenter+0.125 peakBinCenter+0.125],  [-2 10], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    end
    
    

    % ax4: BD posterior probability matrices for 20ms time-binned firings
    ax4 = subplot(12, nICs, ((9:10)-1)*nICs + ic);
 
    nPosBins = size(spatialTunings_LR, 2);

    spatialTunings_LR(spatialTunings_LR==0) = 1e-4;
    spatialTunings_RL(spatialTunings_RL==0) = 1e-4;


    
    binCentersfine = binCenters.fineTS.(currPeriod); 
    
    currFineBinIdx = binCentersfine >= currStartTime & binCentersfine <= currEndTime;
        
    currBinCenters =  binCentersfine(currFineBinIdx);
    currFirings = binnedFiring.fineTS.(currPeriod){1,2}(:, currFineBinIdx);

    
    sumFiring = sum(currFirings, 1);
    idx       = find(sumFiring == 0);


    postPr  = struct('LR', [], 'RL', [], 'integrated', []);

    postPr.LR = baysDecoder(currFirings, spatialTunings_LR, binDur.fineTS);
    postPr.RL = baysDecoder(currFirings, spatialTunings_RL, binDur.fineTS);

    postPr.integrated = (postPr.LR + postPr.RL) ./ repmat(sum(postPr.LR + postPr.RL), [nPosBins, 1]);
    postPr.integrated(:, idx) = 0;


    postPr.LR = postPr.LR ./ repmat(sum(postPr.LR, 1), [size(postPr.LR, 1) 1]);
    postPr.RL = postPr.RL ./ repmat(sum(postPr.RL, 1), [size(postPr.RL, 1) 1]);

    postPr.LR(:, idx) = 0;
    postPr.RL(:, idx) = 0;

    
    currPosterior = postPr.integrated;
    
%     currBinCenters = (startTimes(iperiod):binDur.fineTS:(startTimes(iperiod) + duration2Plot -binDur.fineTS)) + binDur.fineTS/2;
    
    cl(1) = min(currPosterior(:));
    cl(2) = min(currPosterior(:)) + 0.1*(max(currPosterior(:))-min(currPosterior(:)));
    
    imagesc(currBinCenters-currStartTime, 1:size(currPosterior, 1), currPosterior, cl); colormap(ax4, 'hot')
    set(gca, 'YDir', 'normal')

    hold on
    line([flankingPeriod flankingPeriod]*1e-3, [0 nPosBins], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([currEndTime-currStartTime-flankingPeriod*1e-3  currEndTime-currStartTime-flankingPeriod*1e-3], [0 nPosBins], 'linestyle', ':', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    
    if maxreact > 20
    hold on
    line([peakBinCenter-0.125 peakBinCenter-0.125],  [0 nPosBins], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    hold on
    line([peakBinCenter+0.125 peakBinCenter+0.125],  [0 nPosBins], 'linestyle', '-', 'linewidth', 2, 'color', [0 0.4 0.8 0.5])
    end
    
    linkaxes([ax1 ax2 ax3 ax4], 'x')

    xlim([0 currEndTime-currStartTime])
    
    end
    
    savepdf(gcf, ['hiActEvents_' num2str(ii)], '-dpdf')
    close all
end



end

