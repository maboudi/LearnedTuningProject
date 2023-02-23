clear;
clc;


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults';

rr = dir(parentDir);

binSizes = [0.015 0.02 0.025 0.05 0.125 0.25];



for sessionNumber = 1:17

sessionName = rr(sessionNumber+2).name


basePath = fullfile(parentDir, sessionName);

s = load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'fileInfo', 'spikes_pyr');
fileInfo = s.fileInfo;
spikes   = s.spikes_pyr;
nUnits   = numel(spikes);

clear s




fileName = fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']);
if isfile(fileName)
    load(fileName, 'clusterQuality')
else 
    clusterQuality = [];
end



behavior = fileInfo.behavior;
linearPos = fileInfo.linearPos;
posLapIdx = linearPos(:, 3);

uniqLapIdx = unique(posLapIdx);
uniqLapIdx = uniqLapIdx(uniqLapIdx > 0);

nLaps = numel(uniqLapIdx);


laps = zeros(nLaps, 3);

for iLap = 1:nLaps
    
    startPosIdx = find(posLapIdx == uniqLapIdx(iLap), 1, 'first');
    endPosIdx   = find(posLapIdx == uniqLapIdx(iLap), 1, 'last');
    
    laps(iLap, 1) = fileInfo.linearPos(startPosIdx, 2);
    laps(iLap, 2) = fileInfo.linearPos(endPosIdx, 2);
    laps(iLap, 3) = uniqLapIdx(iLap);
    
end

laps  = laps(laps(:,2) < behavior.time(2,2), :);
nLaps = size(laps, 1);


if isfield(spikes(1).spatialTuning_smoothed, 'RL') 
    twoPlaceFieldsFlag = true;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.RL);
    directions = {'RL'; 'LR'; 'uni'};
else
    twoPlaceFieldsFlag = false;
    nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
    directions = {'uni'};

end



% config params

posBinSize = 2; 
binDur = 0.02;
runSpeedThresh = 10;
ifshuffle = 0;

thetaPeriods = [];
turningPeriods = [];

% gw = gausswindow(3,9);


%% the cross-validation procedure

% laps = laps(randperm(nLaps), :); % randomizing order of the laps for the cross-validation procedure
% laps(:, 3) = 1: nLaps;

nFolds = min(nLaps, 10);
foldSize = floor(nLaps/nFolds);

binnedActiveRun  = cell(nLaps, 1);
posteriorExclude = cell(nLaps, 1);


for fold = 1:nFolds
  
%    if fold == nFolds
%        decodeLaps = laps((fold-1)*foldSize+1 : nLaps, :); % laps to test (decode based on) the Bayesian encoding models
%    else
%        decodeLaps = laps((fold-1)*foldSize+1 : fold*foldSize, :);
%    end

    decodeLaps = laps;


   encodeLaps = laps; %(~ismember(1:nLaps, decodeLaps), :); 
   
  
   % calculate the spatial tunings
   
   fileInfo.timeUnit = 1;
   
   for iDir = 1:numel(directions) 
        spikes = spatialTuning_1D_V2(spikes, encodeLaps, thetaPeriods, turningPeriods, runSpeedThresh, directions{iDir} , posBinSize, fileInfo); % spatial tuning or Bayesian encoding model
   end
   
   
   nDecodeLaps = size(decodeLaps, 1);
   
   for iLap = 1:nDecodeLaps
       
        currLap = decodeLaps(iLap, 3);
        binnedActiveRun{currLap}  = binRunData_1D(spikes, decodeLaps(iLap, 1:2), runSpeedThresh, turningPeriods, binDur, 0, posBinSize, fileInfo); % need to add theta and turning periods to this ...
        posteriorExclude(currLap) = calculateAssemblyTuning_theta(binnedActiveRun(currLap), spikes, clusterQuality, binDur, ifshuffle);

   end   

   
end


%%

posteriorExclude_cnct = cell2mat(posteriorExclude');
binnedActiveRun_cnct  = cell2mat(binnedActiveRun');
nTimeBins = size(binnedActiveRun_cnct, 2);

learnedTunings = nan(nUnits, nPosBins);
spatialTunings = zeros(nUnits, nPosBins);

for iUnit = 1:nUnits
    
    currUnitPosteriors = posteriorExclude_cnct(:, :, iUnit); 
    curUnitFirings     = binnedActiveRun_cnct(iUnit, :);
    
    weightedSum = nansum(currUnitPosteriors .* repmat(curUnitFirings, [nPosBins 1]), 2)/nTimeBins;
    px = nansum(currUnitPosteriors, 2)/nTimeBins;
    

    spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni; 
    spatialTunings(iUnit, :) = spatialTunings(iUnit, :) ./ max(spatialTunings(iUnit, :));

    learnedTunings(iUnit, :) = weightedSum./px;
%     temp = conv(temp, gw, 'same');
       
end


allCorrs = corr(learnedTunings', spatialTunings');
LTPFcorr = diag(allCorrs);


save(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_activeRun_binDur0.02.mat']), 'learnedTunings')


end


return

%% figures


% non-directional spatial tuning
% spatialTunings = zeros(nUnits, nPosBins);
peakPFLocation = nan(nUnits, 1);

for iUnit = 1:nUnits
%     spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni; 
%     spatialTunings(iUnit, :) = spatialTunings(iUnit, :) ./ max(spatialTunings(iUnit, :));
    [~, peakPFLocation(iUnit)] = max(spatialTunings{sessionNumber}(iUnit, :));
end


learnedTunings{sessionNumber} = learnedTunings{sessionNumber} ./repmat(max(learnedTunings{sessionNumber}, [], 2), 1, nPosBins);


[~, sortIdx] = sort(peakPFLocation, 'ascend');

figure; 

subplot(1,2,1)
imagesc(spatialTunings{sessionNumber}(sortIdx, :))


subplot(1,2,2)
imagesc(learnedTunings{sessionNumber}(sortIdx, :))

colormap('jet')



