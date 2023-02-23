function basicProcedure_GrosmarkOriginal(sessionNumber)


% addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/ReplayPreplayAnalyses'))

% currDir = '/nfs/turbo/umms-kdiba/Kourosh/Grosmark_originalClusters';

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_originalClusters';
cd(currDir)


%% loading session data

% VarList = {'Spikes','Epochs','Position'};

% note that the spikes within the .mat file are only from the sorted units,
% for MUA I need to load the clu and res files again to consider all the
% spikes

dirresults = dir(currDir);

noSessions = 8;
sessionNames = cell(noSessions, 1);


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100

mazeLimits = [-inf inf -inf 13 ;... % Achilles_10252013
              -inf inf -inf inf;... % Achilles_11012013
              -inf inf -140 inf;... % Buddy_06272013
              -inf inf -inf inf;... % Cicero_09012014
              -inf inf -inf inf;... % Cicero_09102014
              -inf inf -inf inf;... % Cicero_09172014
              -inf inf -inf inf;... % Gatsby_08022013
              -20  inf -inf 110];   % Gatsby_08282013

mazeShapes = {'linear', 'circular', 'linear', 'linear', 'circular', 'linear', 'linear', 'circular'};

    
    
sessionName = dirresults(sessionNumber + 2).name;

load(fullfile(currDir, 'NoveltySessInfoMatFiles', [sessionName '_sessInfo.mat']))


spikes = sessInfo.Spikes;
behavior = sessInfo.Epochs;
position = sessInfo.Position;


behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity





%% sessioninfo

fileBase = fullfile(currDir, sessionName, sessionName);

fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
mkdir(mainDir)

% Load the recording configurations from the xml file

Par = LoadXml(fileBase);
% fileinfo.lfpSampleRate = Par.lfpSampleRate;

fileinfo.nCh = Par.nChannels;
fileinfo.Fs = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording
fileinfo.lfpSampleRate = 1;



% RUN high power theta periods

thetaPeriods = importEVTfile(fileBase); % make sure that we know details about procedures and parameters which were used for the detection of these periods 


% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

tpos = position.TimeStamps';


% We need to exclude long periods during which the positional data is
% missing for PBE detection (if we are not using theta)


temp   = [0; diff(isnan(xpos))]; 
starts = find(temp == 1);
ends   = find(temp == -1);

if ends(1) < starts(1)
    starts = [1; starts];
end

if starts(end) > ends(end)
    ends = [ends; length(xpos)];
end

nanPeriods = [tpos(starts) tpos(ends)];
nanDurs    = diff(nanPeriods')'; 
nanPeriods = nanPeriods(nanDurs > 5, :);


fileinfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track boundaries

fileinfo.xyt(find(ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3) | xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1)), 1:2) = NaN;


linearPos = linearizePosition(fileinfo, behavior, mazeShapes{sessionNumber});

fileinfo.xyt2(:, 1) = linearPos;
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


% % refining the boundaries of RUN period
% 
% firstTimePnt = find(~isnan(linearPos), 1, 'first');
% 
% fileinfo.xyt2(:, 1) = linearPos(firstTimePnt: end); 
% fileinfo.xyt2(:, 2) = fileinfo.xyt(firstTimePnt: end, 3);
% 
% 
% % The following line of code was added to remove the Nans introduced
% % because of excluding the part before the animal was put on the track
% 
% behavior.time(2,1) = behavior.time(2,1) + firstTimePnt * diff(fileinfo.xyt2(1:2,2));


% % updated
if ismember(sessionNumber, [2 5 8]) 
    direction = 'uni'; 
    [laps, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
else
    direction = 'bi';
    [lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
    
    if length(lapsStruct.RL) > length(lapsStruct.LR)
       lapsStruct.RL(1,:) = [];
%        behavior.MazeEpoch(1,:) = lapsStruct.LR(1,:);
    end
    
    totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
    laps = zeros(totNumLaps, 2);
    laps(1:2:totNumLaps, :)  = lapsStruct.LR;
    laps(2:2:totNumLaps, :)  = lapsStruct.RL;
end

laps(:, 3) = 1:size(laps, 1); 

% %


fileinfo.xyt2(:, 3) = zeros(size(fileinfo.xyt2(:, 1))); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end

runSpeedThresh = 10; % cm/s

%% Behavioral states

% nrem = 1, Drowsy = 3.5; rem = 2, Intermediate = 1.5, wake = 4


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Can we consider them as NREM???


bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ... % including both the active and quiet wake periods
            3 *   ones(length(behavior.Drowsy), 1); ... % considering Drowsy as a part of NREM, for now 
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering Intermediate as a part of NREM, for now 
            2   *   ones(length(behavior.REM), 1)];



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

spikeStruct = spikeBehaviorAnalysis3(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo);



% We might use multiunit (any detected spike) as well for defining the population burst events
% REMOVE CLUSTER #0 as it might contain artifacts 

numberofShanks = length(dir([fileBase '.res.*']));

MUA.t = [];
for shank = 1:numberofShanks
    
    resContent  = load([fileBase '.res.' num2str(shank)]);
    
    cluContent  = load([fileBase '.clu.' num2str(shank)]);
    clusterIDs  = cluContent(2:end);
    
    spikes2Add = resContent(clusterIDs ~= 0);
    
    MUA.t = [MUA.t; spikes2Add]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');
% The originla sampling frequency is 20000

MUA.t = MUA.t./20000;



%% Spatial tuning of the units


close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)

% decide whether to exclude the turning periods from calculations


%%% 2D spatial tuning

% This section needs to be updated with changes similar to what I did for
% spatialTuning_1D
% 
% if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular
%     spatialTunings2D = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'uni', 1, runSpeedThresh, fileinfo.Fs, subfolder);
% else
%     spatialTunings2D_LR = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'LR', 1, runSpeedThresh, fileinfo.Fs, subfolder);
%     spatialTunings2D_RL = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'RL', 1, runSpeedThresh, fileinfo.Fs, subfolder);
% end


%%% 1D spatial tuning: using linearized position

% the place fields are calculated separately for the left and right
% direction


if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular

    [spatialTunings, unsmoothedSpatialTunings, PF_sorted, runTemplate, spatialInfo, conslapsRatio, diffWithAvg] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings, [], spatialInfo, [], conslapsRatio, [], behavior, [], fileinfo, subfolder);
else

    [spatialTunings_LR, unsmoothedSpatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
    [spatialTunings_RL, unsmoothedSpatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileinfo, subfolder);

    [spatialTunings_biDir, unsmoothedSpatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, ccl);
end

save(fullfile(subfolder, 'biDirectional.mat'), 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir')
% close all




%% Ripple Detection

% After visualizing deteceted events using different metrics (MUA and
% ripple ampplitude) on Neuroscope, I realized that ripple is more reliable


best_channels = BigRippleChannels(fileBase, fileinfo);
fileinfo.RippleChannels = best_channels(1);
%%
threshSD = 3;
[ripplePeriods, rippleLFP, ripplePower] = RippleDetect(fileBase, fileinfo, threshSD); % the results in seconds


%% PBEs


%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the mean


% baseFile = fullfile(unsortedDataDir, sessionName, sessionName);
% 
% thetaPeriods = load([baseFile '.theta.1']);
% thetaPeriods = floor(thetaPeriods./fileinfo.lfpSampleRate);


fileinfo.lfpSampleRate = 1;

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(3,2);


subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, wake=4


velocityFilter = 1;
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);


qclus = [1 2 3];

[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileinfo); 
PBErippleIdx = ifContainRipples(secondaryPBEs, ripplePeriods);

secondaryPBEs(:, 5) = PBErippleIdx;

nPBEs = size(binnedPBEs, 1);


secondaryPBEs = [secondaryPBEs zeros(nPBEs, 4)];

for ii= 1:nPBEs
    
    pbeCenter = secondaryPBEs(ii, 3);
    
    boutInd        = find(bvrTimeList(:,1) < pbeCenter & bvrTimeList(:,2) > pbeCenter, 1, 'first');
    boutBrainState = bvrState(boutInd);
    
    secondaryPBEs(ii, 5 + boutBrainState) = 1;
    
end


baseStruct = struct('data', [], 'p', [], 'ts', [], 'pts', []);


PREbinnedPBEs       = baseStruct;
PREidx              = find(secondaryPBEs(:, 1) > behavior.time(1,1) & secondaryPBEs(:, 2) < behavior.time(1,2));
PREbinnedPBEs.data  = binnedPBEs(PREidx, :);
secondaryPBEs_PRE   = secondaryPBEs(PREidx, :);
PREbinnedPBEs       = genSurrogates(PREbinnedPBEs);


RUNbinnedPBEs       = baseStruct;
RUNidx              = find(secondaryPBEs(:, 1) > behavior.time(2,1) & secondaryPBEs(:, 2) < behavior.time(2,2));
RUNbinnedPBEs.data  = binnedPBEs(RUNidx, :);
secondaryPBEs_RUN   = secondaryPBEs(RUNidx, :);
RUNbinnedPBEs       = genSurrogates(RUNbinnedPBEs);


POSTbinnedPBEs       = baseStruct;
POSTidx              = find(secondaryPBEs(:, 1) > behavior.time(3,1) & secondaryPBEs(:, 2) < behavior.time(3,2));
POSTbinnedPBEs.data  = binnedPBEs(POSTidx, :);
secondaryPBEs_POST   = secondaryPBEs(POSTidx, :);
POSTbinnedPBEs       = genSurrogates(POSTbinnedPBEs);


BayesianReplayDetection_GrosmarkReclu_nov2020


end


function speed = calspeed(position)


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100


% 2D position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

timepnts = position.TimeStamps';


diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difftt = [1; diff(timepnts)];


% calculation of speed in just x direction (we then can use the speed or theta for filtering the events)

velocity = sqrt(diffx.^2 + diffy.^2)./difftt; % 2D velocity 

velocity(isnan(velocity)) = interp1(timepnts(~isnan(velocity)), velocity(~isnan(velocity)), timepnts(isnan(velocity)));

% smoothing the speed

sigma = 20; %%% smoothing the speed, the positon sampling rate is around 40 Hz (0.0256 sec sampling period), so duration of the sigma is about 25.6 ms times sigma (25.6*20 ~ 512 ms)
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
velocity = conv(velocity, smoothwin, 'same'); 


speed.v = velocity;
speed.t = timepnts;


end



function [assemblyTunings, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
%     [assemblyTunings, assemblyTunings_ci, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;
    p_of_x = sum(posteriorProbMatrix, 2)/nTimeBins;
    
    
    assemblyTunings = weightedSummation; 
    
    
%     sigma2 = sum(repmat(unitFirings, [nPosBins 1]) .* ((posteriorProbMatrix ./ repmat(p_of_x, [1 nTimeBins])) - assemblyTunings).^2, 2)/ nTimeBins;
%     sigma = sqrt(sigma2);
%     
%     assemblyTunings_ci = assemblyTunings - sigma/sqrt(nTimeBins) * 1.645;
    
    
end
