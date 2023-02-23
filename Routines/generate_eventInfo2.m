function generate_eventInfo2(sessionNumber)


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

addpath(genpath(fullfile(parentDir, '/ReplayPreplayAnalyses')))

datasetDir = fullfile(parentDir, '/Grosmark_reclustered');
cluFilesDir = '/data/datasets/GrossmarkData_BapunsClusters';


stateDetectionResultsDir = fullfile(parentDir, 'StateDetectionResults');

outputDir = fullfile(parentDir, 'sessions_calculated_PBEinfo3');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end



%% loading session data

VarList = {'spikes','behavior','position'};

% note that the spikes within the .mat file are from just the pyramidal
% units, for MUA I need to load the clu and res files again which contains
% all instable units, multiunits, and also noise cluster which should be
% excluded

% Although no difference regarding the behavior and position data


for var = 1 : length(VarList)
    load(fullfile(datasetDir, ['/CRCNSReclustered-' VarList{var} '.mat']))
end

sessionNames = fieldnames(spikes);



% MAZE SHAPES

% the commented ones are those that were not used in the reculstered version of the dataset

mazeLimits = [-inf inf -inf 13; ... % Achilles_10252013
              -inf inf -inf inf;... % Achilles_11012013
              -inf inf -140 inf;... % Buddy_06272013 
%               -inf inf -inf inf;... % Cicero_09012014
%               -15 inf -inf 115;... % Cicero_09102014
              -inf inf -inf inf;... % Cicero_09172014
%               -inf inf -inf inf;... % Gatsby_08022013
              -20  inf -inf 110];   % Gatsby_08282013 
          
mazeShapes = {
            'linear'; ... % Achilles_10252013
            'circular'; ... % Achilles_11012013
            'linear'; ... % Buddy_06272013
%             'linear'; ... % Cicero_09012014
%             'circular'; ... % Cicero_09102014
            'linear'; ... % Cicero_09172014
%             'linear'; ... % Gatsby_08022013
            'circular' ... % Gatsby_08282013
            };

        
        
% CURRENT SESSION        
        
sessionName = sessionNames{sessionNumber}

baseFile = fullfile(cluFilesDir, sessionName, sessionName);


storageDir = fullfile(outputDir,  sessionName);
if ~exist(storageDir, 'dir'); mkdir(storageDir); end


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['Position.' sessionName]);


behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];


% give a warning if there is a gap between the end of one epoch and start
% of the next epoch

if behavior.time(2,1) ~= behavior.time(1,2) || behavior.time(3,1) ~= behavior.time(2,2) 
    warning('There is gap between epochs')
end



%% BRAIN STATES / THETA-NONTHETA

% We are not going to calculate the periods belonging to each brain state, as the calculated boundaries lacked enough accuracy for the purpose of PBE detection. Instead we use
% [almost instantaneous] estimations of theta/delta ratio and slow wave amplitude to estimate the state in which each population burst event occurred  


% loading theta info
filename = sprintf([sessionName '.thetaStates.mat']);
load(fullfile(stateDetectionResultsDir, sessionName, filename), 'timePnts', 'theratio'); 
brainState.thetaTimePnts = timePnts;
brainState.theratio      = theratio;


% loading slow wave power/slope 
filename = sprintf([sessionName '.SleepScoreMetrics.mat']);
load(fullfile(stateDetectionResultsDir, sessionName, filename), 'SleepScoreMetrics')


brainState.slowWave        = SleepScoreMetrics.broadbandSlowWave;
brainState.sw_emg_timePnts = SleepScoreMetrics.t_clus;

try 
  load(fullfile(stateDetectionResultsDir, sessionName, sprintf([sessionName '.swthresh.mat'])), 'swthresh'); 
    brainState.swthresh = swthresh;
catch
    brainState.swthresh = SleepScoreMetrics.histsandthreshs.swthresh;
end



% loading emg
brainState.emg = SleepScoreMetrics.EMG;

try
    load(fullfile(stateDetectionResultsDir, sessionName, sprintf([sessionName '.emgThresh.mat'])), 'emgThresh');
    brainState.emgThresh = emgThresh;
catch
    brainState.emgThresh = SleepScoreMetrics.histsandthreshs.EMGthresh; 
end



%%% temporary
fileName = sprintf('%s.brainState.mat', sessionName);
save(fullfile(parentDir, 'assemblyTuning_finalResults', sessionName, 'spikes', fileName), 'brainState')

return


%%%

% bad channels to exclude for lpf event detections (e.g., ripples)

filename = sprintf('%s.sessionInfo.mat', sessionName);
load(fullfile(stateDetectionResultsDir, sessionName, filename), 'sessionInfo')

badChannels = sessionInfo.badchannels;

%% SESSION'S INFORMATION


% load spikes and fileInfo
folderName = fullfile(storageDir,  'spikes');
fileName = fullfile(folderName, [sessionName '.spikes.mat']);

if exist(fileName, 'file')

    load(fileName, 'spikes', 'spikes_pyr', 'fileInfo') % , 'clusterQuality'

else
    mkdir(folderName)


    Par = LoadXml(baseFile); % Load the recording configurations from the xml file

    fileInfo.sessionName   = sessionName;
    fileInfo.animal        = sessionName(1:strfind(sessionName, '_')-1);
    fileInfo.xyt           = [];
    fileInfo.linearPos     = []; % linearized position and lap information=
    fileInfo.speed         = [];
    fileInfo.tbegin        = behavior.time(1,1);
    fileInfo.tend          = behavior.time(3,2);
    fileInfo.Fs            = Par.SampleRate; 
    fileInfo.lfpSampleRate = Par.lfpSampleRate;
    fileInfo.timeUnit      = 1; % in sec
    fileInfo.nCh           = Par.nChannels;
    fileInfo.badChannels   = badChannels;
    fileInfo.behavior      = behavior;


    
    % ripples
    if ~isfield(fileInfo, 'RippleChannels')
       best_channels = BigRippleChannels(baseFile, fileInfo);
       fileInfo.RippleChannels = best_channels; % check if there are no bad channels among the best_channels 
        
       save(fileName, 'fileInfo') % the value of fileInfo will be updated later
    end



    % Position data

    xpos = position.TwoDLocation(:,1)* 100; % convert to centimeters
    ypos = position.TwoDLocation(:,2)* 100;
    tpos = position.TimeStamps';

    fileInfo.xyt = [xpos ypos tpos]; 
    position_samplingRate = 1/(tpos(2) - tpos(1));

    speed = calspeed(fileInfo, position_samplingRate); % animal's velocity

    fileInfo.speed = speed;



    %% LINEARIZE POSTION AND CALCULATE LAPS

    fileInfo.xyt(xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1) | ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3), :) = [];


    linearPos = linearizePosition(fileInfo, behavior, mazeShapes{sessionNumber});


    fileInfo.linearPos(:, 1) = linearPos; 
    fileInfo.linearPos(:, 2) = fileInfo.xyt(:, 3);


    if ismember(sessionNumber, [2 5 8]) 
        direction = 'uni'; 
        laps = calculateLapTimings(fileInfo, direction, storageDir); 
    else
        direction = 'bi';
        lapsStruct = calculateLapTimings(fileInfo, direction, storageDir); 

        if length(lapsStruct.RL) > length(lapsStruct.LR)
           lapsStruct.RL(1,:) = [];
        end

        totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
        laps = zeros(totNumLaps, 2);
        laps(1:2:totNumLaps, :)  = lapsStruct.LR;
        laps(2:2:totNumLaps, :)  = lapsStruct.RL;

    end

    laps(:, 3) = 1:size(laps, 1); 
    

    fileInfo.linearPos(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

    for ii = 1: length(laps)
       idx =  fileInfo.linearPos(:, 2) > laps(ii, 1) & fileInfo.linearPos(:, 2) < laps(ii, 2);
       fileInfo.linearPos(idx, 3) = laps(ii, 3);
    end

    runSpeedThresh = 10; % cm/s


    %% calculating position, velocity, etc at each spike (these will be needed for calculation of place fields)


    nUnits = numel(spikes);

    for unit = 1:nUnits

        spikeTimes = spikes(unit).time;
        spikes(unit).x = interp1(fileInfo.xyt(:, 3), fileInfo.xyt(:, 1), spikeTimes);
        spikes(unit).y = interp1(fileInfo.xyt(:, 3), fileInfo.xyt(:, 2), spikeTimes);
        spikes(unit).linearPos = interp1(fileInfo.linearPos(:, 2), fileInfo.linearPos(:, 1), spikeTimes);

        spikes(unit).speed = interp1(speed.t, speed.v, spikeTimes);

        spikes(unit).lap = zeros(numel(spikeTimes), 1);
        for spk = 1:numel(spikeTimes)
            index = find(laps(:, 1) < spikeTimes(spk) & laps(:, 2) > spikeTimes(spk));
            if ~isempty(index)
                spikes(unit).lap(spk) = laps(index, 3);
            end
        end

    end


    
    %% cluster quality

    clusterQuality = calClusterQuality(baseFile);


    %% Spatial tuning of the units
    
    spikes_pyr = spikes([spikes.quality] <= 3); % limiting to only pyramidal units
    fileInfo.okUnits = 1:numel(spikes_pyr); % all pyramidal units regardless of their stability will be considered in Bayesian decoding

  
    close all

    subfolder = fullfile(storageDir, 'placeFields');
    if ~exist(subfolder, 'dir')
        mkdir(subfolder)
    end

    %%% 1D spatial tuning: using linearized position

    posBinSize = 2;

    if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular

        spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);

    else 
        spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'LR', posBinSize, fileInfo, subfolder);
        spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'RL', posBinSize, fileInfo, subfolder);

        spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);
    end

    close all
    
    fileName = fullfile(folderName, [sessionName '.spikes.mat']);
    save(fileName, 'spikes', 'spikes_pyr', 'clusterQuality', 'fileInfo', '-v7.3')
    
    
    
end



%% generate a clu and res file for interneurons, basically all the clusters excluding 0 and 1 and pyramidal units

subfolder = fullfile(storageDir, 'placeFields');
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end


pyrUnitIDs = [spikes.id]';
pyrUnitIDs = reshape(pyrUnitIDs, [2 numel(pyrUnitIDs)/2])';

nShanks    = length(dir([baseFile '.res.*']));

res_t = [];
clu_t = [];

tt = 0;

for shank = 1:nShanks
    
    resContent = load([baseFile '.res.' num2str(shank)]);
    cluContent = load([baseFile '.clu.' num2str(shank)]);
    
    clusterIDs = cluContent(2:end); % the first array is the number of clusters
     
    
    uniqClusters = unique(clusterIDs);
    intClusters  = setdiff(uniqClusters, [0 pyrUnitIDs(2, pyrUnitIDs(1,:) == shank)]); % the noise cluster and cluster #1 (unsorted spikes) were excluded
    
    nIntClusters = numel(intClusters);
    
    
    idx = ismember(clusterIDs, intClusters);
    
    resContent = resContent(idx);
    clusterIDs = clusterIDs(idx);
    
    
    % renumber the clusters
    newClusterIDs = zeros(size(clusterIDs));
    for ii = 1:nIntClusters
        newClusterIDs(clusterIDs == uniqClusters(ii)) = ii+tt;
    end
    
    tt = tt + nIntClusters;
    
    res_t = [res_t; resContent];
    clu_t = [clu_t; newClusterIDs];
    
end
    
[res_t, ind] = sort(res_t);
clu_t = clu_t(ind);  


Saveres([subfolder '/' fileInfo.sessionName 'Int.201'  '.res.3'], res_t);
SaveClu([subfolder '/' fileInfo.sessionName 'Int.201'  '.clu.3'], clu_t);




%% RIPPLE DETECTION
       
folderName = fullfile(storageDir,  'ripples');
fileName   = fullfile(folderName, [fileInfo.sessionName '.rippleEvents.mat']);

if exist(fileName, 'file')
    load(fileName, 'rippleLFP', 'rplDetectParams', 'rippleEvents')
end



overwrite = 0;
if ~exist('rippleLFP', 'var') || overwrite == 1
    
mkdir(folderName);
    
    rplDetectParams.method                 = 'summedChan';
    rplDetectParams.threshold              = 2;

    rplDetectParams.minDur                 = 0.04; % number of bins
    rplDetectParams.maxDur                 = 0.60;
    rplDetectParams.binDur                 = 0.02;
    
    rplDetectParams.thRatioThresh          = 1;
    rplDetectParams.rippleAThresh_low      = 1;
    rplDetectParams.rippleAThresh_high     = 3;
    
    
    
    [ripplePeriods, rawLFP, ripplePower_summed] = RippleDetect(baseFile, fileInfo, rplDetectParams.threshold); % the results in seconds
    

    % downsample the ripple band lfp or power to 1kHz
    
    ts    = 1/fileInfo.lfpSampleRate;

    tend  = length(rawLFP)/fileInfo.lfpSampleRate;
    tPnts = ts:ts:tend;

    t_1ms    = 1e-3:1e-3:tend;

    
    rawLFP2   = zeros(numel(t_1ms), size(rawLFP, 2));
    for ich = 1:size(rawLFP, 2)
        rawLFP2(:, ich)   = interp1(tPnts, rawLFP(:, ich), t_1ms);
    end
    rawLFP = rawLFP2;
    
    
    power_summed = interp1(tPnts, ripplePower_summed, t_1ms);
   
    rippleLFP.rawLFP       = rawLFP;
    rippleLFP.power_summed = power_summed;
    rippleLFP.timePnts     = t_1ms;
      
end


%% DETECTING POPULATION BURST EVENTS (PBEs)

folderName = fullfile(storageDir,  'PopulationBurstEvents');
overwrite = 1;

if exist(folderName, 'dir') && overwrite == 0
    
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo.mat']), 'PBEInfo_Bayesian', 'sdat');
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo.mat']), 'PBEInfo')
    
else
        
    mkdir(folderName)
    
    
    % PYR and MUA for detection of PBEs
    % MUA
    nShanks = length(dir([baseFile '.res.*']));

    MUA.time = [];
    for shank = 1:nShanks

        resContent  = load([baseFile '.res.' num2str(shank)]);

        cluContent  = load([baseFile '.clu.' num2str(shank)]);
        clusterIDs  = cluContent(2:end); % the first entry is the total number of clusters

        spikes2Add = resContent(~ismember(clusterIDs, 0)); % to exclude the noise cluster 

        MUA.time = [MUA.time; spikes2Add]; % the first is the total number of clusters detected for the shank   
    end
    MUA.time = sort(MUA.time, 'ascend');
    MUA.time = MUA.time./fileInfo.Fs;



    % detection parameters
    
    pbeDetectParams.time_resolution        = 0.001;
    pbeDetectParams.threshold              = 2;
    pbeDetectParams.exclude                = [];

    pbeDetectParams.smoothingSigma         = 0.02;
    
    pbeDetectParams.minDur                 = 0.04; % number of bins
    pbeDetectParams.maxDur                 = 0.70;
    
    pbeDetectParams.maxSilence             = 0.05;
    
    pbeDetectParams.maxVelocity            = 10; 
    pbeDetectParams.thRatioThresh          = 1;
    pbeDetectParams.rippleAThresh_low      = 1;
    pbeDetectParams.rippleAThresh_high     = 3;
    
    pbeDetectParams.binDur                 = 0.02;
    pbeDetectParams.minNTimeBins           = 4;
    pbeDetectParams.minFiringUnits         = max(5, 0.1*numel(fileInfo.okUnits)); % okUnits are all pyramidal units regardless of their stability
    

    
    % detecting primary PBEs, just the boundaries of potential PBEs before
    % applying any criteria

    [primaryPBEs, sdat] = PBPeriods(MUA, pbeDetectParams, fileInfo); % already filtering based on the duration in this step
    
    fileName = fullfile(folderName, [fileInfo.sessionName '.PBEs_before_splitting.mat']);
    save(fileName, 'primaryPBEs', 'pbeDetectParams', 'sdat', '-v7.3')
    
    
    % splitting the detected PBEs if there was a long silence in the middle
    PBEs_splitted = splitPBEs(primaryPBEs, spikes_pyr, sdat, pbeDetectParams);
    

    
    % calculate a number of features for the PBEs; brain state, number of
    % acive units, binned firing rate , etc
    
    [PBEInfo, PBEInfo_Bayesian, idx_Bayesian] = genEventStructure(PBEs_splitted, spikes_pyr, sdat, rippleLFP, pbeDetectParams, brainState, fileInfo, folderName, 'pbe');
    
end



%% calculate for each ripple concurrent brain state, binned firing rate, number of active units, etc

folderName = fullfile(storageDir,  'ripples');
fileName   = fullfile(folderName, [fileInfo.sessionName '.rippleEvents.mat']);

if ~exist(fileName, 'file') 
    
    rippleEvents = genEventStructure(ripplePeriods, spikes_pyr, sdat, rippleLFP, rplDetectParams, brainState, fileInfo, folderName, 'ripple');
    
    save(fileName, 'rippleLFP', 'rplDetectParams', 'rippleEvents', '-v7.3')
end



%% adding variables to the EEG file
% 
% nCh = fileInfo.nCh;
% 
% 
% EEGdata = readmulti([baseFile '.eeg'], nCh, 1); % load .eeg trace
% totalT  = length(EEGdata)/fileInfo.lfpSampleRate;
% 
% 
% samplingDur  = 1/fileInfo.lfpSampleRate;
% samplingPnts = samplingDur:samplingDur:totalT;
% 
% sdat_mua_2add = interp1((1:length(sdat))*1e-3, sdat, samplingPnts);
% speed_2add    = interp1(fileInfo.speed.t, fileInfo.speed.v, samplingPnts);
% ripple_2add   = interp1(rippleLFP.timePnts, rippleLFP.power_summed, samplingPnts);
% theratio_2add = interp1(brainState.thetaTimePnts, brainState.theratio, samplingPnts);
% 
% 
% sdat_mua_2add(isnan(sdat_mua_2add)) = 0;
% sdat_mua_2add(sdat_mua_2add < 0) = 0;
% sdat_mua_2add(sdat_mua_2add > 2) = 2;
% sdat_mua_2add = (2*(sdat_mua_2add-min(sdat_mua_2add))/(max(sdat_mua_2add)-min(sdat_mua_2add))-1)*2^10;
% 
% 
% 
% speed_2add(speed_2add < 10) = 0;
% speed_2add = (2*(speed_2add-min(speed_2add))/(max(speed_2add)-min(speed_2add))-1)*2^11;
% 
% 
% ripple_2add(isnan(ripple_2add)) = 0;
% ripple_2add(ripple_2add > 4) = 4;
% ripple_2add(ripple_2add < 0) = 0;
% ripple_2add = (2*(ripple_2add - min(ripple_2add))/(max(ripple_2add)-min(ripple_2add))-1)*2^10;
% 
% 
% theratio_2add(isnan(theratio_2add)) = 0;
% theratio_2add(theratio_2add < 1) = 0;
% theratio_2add = (2*(theratio_2add - min(theratio_2add))/(max(theratio_2add)-min(theratio_2add))-1)*2^11;
% 
% 
% 
% inputFileHandle  = fopen([baseFile,'.eeg']);
% outputFileHandle = fopen([baseFile,'-2.eeg'],'w');
% doneFrame = 0;
% bufferSize = 4096;
% while ~feof(inputFileHandle)
%     data = (fread(inputFileHandle,[nCh,bufferSize],'int16'))';
%     
%     for frame=1:size(data,1)
%         fwrite(outputFileHandle,[data(frame,:), sdat_mua_2add(frame+doneFrame), speed_2add(frame+doneFrame), ripple_2add(frame+doneFrame), theratio_2add(frame+doneFrame)]','int16'); 
%     end
%     doneFrame=doneFrame+frame;
% 
% end
% 
% fclose(inputFileHandle);
% fclose(outputFileHandle);
% 

end


function speed = calspeed(fileInfo, position_samplingRate)


xpos = fileInfo.xyt(:, 1);
ypos = fileInfo.xyt(:, 2);

timepnts = fileInfo.xyt(:, 3);


diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difft = [1; diff(timepnts)];


velocity = sqrt(diffx.^2 + diffy.^2)./difft; 

% interpolate the velocty at time point with missing position data
nanIdx = isnan(velocity);
velocity(nanIdx) = interp1(timepnts(~nanIdx), velocity(~nanIdx), timepnts(nanIdx));


% smoothing the speed
sigma = 0.25/(1/position_samplingRate); 
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
speed.v = conv(velocity, smoothwin, 'same'); 


speed.v(nanIdx) = nan;

% speed.v = velocity;
speed.t = timepnts;


end
