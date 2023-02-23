% clear; clc; close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_originalClusters';
cd(currDir)


%% loading session data

% VarList = {'spikes','behavior','position'};

 
% for var = 1 : length(VarList)
%     load([currDir '/CRCNSReclustered-' VarList{var} '.mat'])
% end
 

% note that the spikes within the .mat file are only from the sorted units,
% for MUA I need to load the clu and res files again to consider all the
% spikes

temp = dir(currDir);

noSessions = 8;
sessionNames = cell(noSessions, 1);


% 
% sessionNames = fieldnames(spikes);
% noSessions = numel(sessionNames);

 

cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100

mazeLimits = [-inf inf -inf 13;
   % 0 120 -inf 13 ;... % Achilles_10252013
              -inf inf -inf inf;... % Achilles_11012013
              -inf inf -140 inf;... % Buddy_06272013
              -inf inf -inf inf;... % Cicero_09012014
              -inf inf -inf inf;... % Cicero_09102014
              -inf inf -inf inf;... % Cicero_09172014
              -inf inf -inf inf;... % Gatsby_08022013
              -20  inf -inf 110];   % Gatsby_08282013
          
mazeShapes = {'linear', 'circular', 'linear', 'linear', 'circular', 'linear', 'linear', 'circular'};
 

for currSession = 1%1: noSessions


close all
    
sessionName = temp(currSession + 2).name;

load(fullfile(currDir, 'NoveltySessInfoMatFiles', [sessionName '_sessInfo.mat']))

% sessionName = sessionNames{currSession};

% spikes = eval(['spikes.' sessionName]);
% behavior = eval(['behavior.' sessionName]);
% position = eval(['Position.' sessionName]);

spikes   = sessInfo.Spikes;
behavior = sessInfo.Epochs;
position = sessInfo.Position;

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

% calculate animal's velocity
speed = calspeed(position); 


% Particioning maze period based on animal's velocity, assuming that the animal's
% velocity follows a  bimodal (multimodal) distribution

% runSpeedThresh = multiModalDist(speed.v, 2);

% or we can use a traditinal threshold (less than 5 cm/s for quiet wake and above 10 cm/s for theta)
runSpeedThresh = 10;


%% sessioninfo


fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', [], 'pix2cm', 100); 

% loadDir = fullfile(currDir, sessionName, sessionName);

loadDir = fullfile(currDir, sessionName, sessionName);
mainDir = fullfile(currDir, 'analysesResults', fileinfo.name);
mkdir(mainDir)


% Load the recording configurations from the xml file

Par = LoadXml(loadDir);
% fileinfo.lfpSampleRate = Par.lfpSampleRate; 

fileinfo.lfpSampleRate = 1;

fileinfo.nCh = Par.nChannels;
fileinfo.Fs  = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording


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




% % if the pos includes NaNs use interpolation to replace them 
% 
% nanIdx = find(isnan(xpos));
% 
% xpos = replaceNans(xpos, tpos, nanIdx);
% ypos = replaceNans(ypos, tpos, nanIdx);



fileinfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track


fileinfo.xyt(xpos > mazeLimits(currSession, 2) | xpos < mazeLimits(currSession, 1) | ypos > mazeLimits(currSession, 4) | ypos < mazeLimits(currSession, 3), 1:2) = NaN;


% xpos = fileinfo.xyt(:, 1);
% ypos = fileinfo.xyt(:, 2);

% figure; plot(xpos, ypos, '.', 'markersize', 3)


linearPos = linearizePosition(fileinfo, behavior, mazeShapes{currSession});


fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


laps = calculateLapTimings(fileinfo, speed, mainDir);



fileinfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end


%% Behavioral states

% nrem = 1, Drowsy = 3.5; rem = 2, Intermediate = 1.5, wake = 4


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Can we consider them as NREM???


bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ... % including both the active and quiet wake periods
            1 *   ones(length(behavior.Drowsy), 1); ... % considering Drowsy as a part of NREM, for now 
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering Intermediate as a part of NREM, for now 
            2   *   ones(length(behavior.REM), 1)];


%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

spikeStruct = spikeBehaviorAnalysis3(spikes, speed, laps, qual2consider, fileinfo);

% spikeStruct = spikeBehaviorAnalysis2(spikes, speed, qual2consider, fileinfo);


% We need the multiunit (any detected spike) as well for determining the population burst events

numberofShanks = length(dir([loadDir '.res.*']));

MUA.t = [];
for shank = 1:numberofShanks
    
    currSpikes = load([loadDir '.res.' num2str(shank)]);
    MUA.t = [MUA.t; currSpikes]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');

MUA.t = MUA.t./20000; % The original sampling frequency is 20000



% %% Ripple detection

 
% % Here we intend to detect only ripples that have on average relatively large
% % amplitudes over the all four shanks residing in CA1. Alternatively, we can
% % include also ripples (with usually smaller amplitude) having a high
% amplitude on just one or few shanks.
% 
% 
% 
% best_channels = BigRippleChannels(fileinfo,1);
% fileinfo.bestch = best_channels;
% 
% rippPowThresh = 3; 
% rippleEvents = RippleDetect(fileinfo, rippPowThresh); % in lfpSampleRate (beg peak end normalized amplitude)
% 
% Filename = [mainDir '.spl.evt'];
% MakeEvtFile(rippleEvents(:, 1:3), Filename, {'beg', 'peak', 'end'}, fileinfo.lfpSampleRate, 1)
% 
% 
% eventVelocity = interp1(speed.t, velocity, (rippleEvents(:,2)/fileinfo.lfpSampleRate *1e6) + periods(1,1));
% rippleEvents2 = rippleEvents(eventVelocity < 5, :); % refining the ripples events; excluding the ones that happen during mobile state


%% Spatial tuning of the units

close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)


%%% 2D spatial tuning

% This section needs to be updated with changes similar to what I did for
% spatialTuning_1D

% spatialTunings2D_LR = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, subfolder);
% spatialTunings2D_RL = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, subfolder);



%%% 1D spatial tuning: using linearized position

% the place fields are calculated separately for the left and right
% direction


if ismember(currSession, [2 5 8]) % sessions 2, 5, and 8 are circular
    
    [spatialTunings, PF_sorted, PF_peak, activeUnits_sorted, sparsity] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'uni', 2, runSpeedThresh, fileinfo.Fs, subfolder);

else

    [spatialTunings_LR, PF_sorted_LR, PF_peak_LR, activeUnits_sorted_LR, sparsity_LR] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, subfolder);
    [spatialTunings_RL, PF_sorted_RL, PF_peak_RL, activeUnits_sorted_RL, sparsity_RL] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, subfolder);
end


%% PBEs


%% Determining population burst periods


time_resolution = 0.001; % in second
threshZ = 3; % sdf with 3 std deviation above the mean


% thetaPeriods = load([mainDir '.theta.1']);
% thetaPeriods = floor(thetaPeriods./fileinfo.lfpSampleRate);



%% PRE


fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend   = behavior.time(1,2);

subfolder = fullfile(mainDir, 'PBEs', 'PRE');
mkdir(subfolder)

exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, wake=4


% add to the exclude set all the theta periods (they may overlap with the existing period but doesn't matter)

% exclude = [exclude; thetaPeriods];

[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude      = exclude(sortIdx, :);


[primaryPBEs_pre, sdat_pre] = PBPeriods(MUA, fileinfo, [], time_resolution, threshZ, exclude);


[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, spikeStruct, fileinfo);

savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, fileinfo, subfolder)


plotSurroundMultiunitRate(sdat_pre, secondaryPBEs_pre, fileinfo, subfolder)

[poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE, pooledTSEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE); % generate poission and within-PBE-time-swap surrogate datasets




%% RUN (PBEs + theta sequences)

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend   = behavior.time(2,2);

subfolder = fullfile(mainDir, 'PBEs', 'RUN');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [1 2]),:); % nrem=1, rem=2, quiet=3, wake=4
exclude = [exclude; nanPeriods];

% exclude = [exclude; laps(:,1:2); nanPeriods];

% Add theta periods to exclude 
% exclude = [exclude; thetaPeriods]; 


[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude      = exclude(sortIdx, :);


[primaryPBEs_run, sdat_run] = PBPeriods(MUA, fileinfo, [], time_resolution, 3, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, spikeStruct, fileinfo);

savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, fileinfo, subfolder)

plotSurroundMultiunitRate(sdat_run, secondaryPBEs_run, fileinfo, subfolder)

[poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, pooledTSEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN); % generate poission and within-PBE-time-swap surrogate datasets



%% POST

fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend   = behavior.time(3,2);

subfolder = fullfile(mainDir, 'PBEs', 'POST');
mkdir(subfolder)

exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4


% add theta periods to exclude
% exclude = [exclude; thetaPeriods];
[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude      = exclude(sortIdx, :);


[primaryPBEs_post, sdat_post]= PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude);


[eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, spikeStruct, fileinfo);

savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, fileinfo, subfolder)


plotSurroundMultiunitRate(sdat_post, secondaryPBEs_post, fileinfo, subfolder)

[poissonEventsBinnedfiringPOST, timeSwapEventsBinnedfiringPOST, pooledTSEventsBinnedfiringPOST] = genSurrogates(eventsBinnedfiringPOST); % generate poission and within-PBE-time-swap surrogate datasets



%% RUN (behavioral time scale binning)

runBouts = laps;

%%% binning the spikes within each run bout

binDur = 0.125; 
qclus = [1 2 3]; % Only the pyramidal neurons were included

runBoutsBinnedfiring  = timeBinning(runBouts, spikeStruct, qclus, binDur, fileinfo);
runBoutsBinnedfiring  = removeSideSilentBins(runBoutsBinnedfiring, runBouts(:, 1), binDur, fileinfo);



%%  Bayesian decoding

if ismember(currSession, [2 5 8])

    BayesianDecodingSequenceAnalysis_1D

else

    BayesianDecodingSequenceAnalysis

end



%% Hidden Markov model

HMMsequenceAnalysis


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