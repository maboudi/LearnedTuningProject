function generate_eventInfo_HiroSessions_longDuration(animalID, sessionNumber)


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

addpath(genpath(fullfile(parentDir, '/ReplayPreplayAnalyses')))

datasetDir  = fullfile(parentDir, 'Hiro_Dataset/Hiro_10hrPost_reclustered/');
cluFilesDir = '/data/datasets/HiroDatasets_cluFiles';


outputDir = fullfile(parentDir, 'sessions_calculated_PBEinfo3');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end



%% loading session data

rats = {'Roy','Ted', 'Kevin'};
% allsessionNumbers = {[1 2 3],[1 2 3], [1]}; 

mazeShape.Roy = {'linear';'linear';'linear'};
mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
mazeShape.Kevin = {'linear'};




%% session info

rat = rats{animalID};

%%% load the data

sessionName = [rat 'Maze' num2str(sessionNumber)];
sessionName2 = [rat '-maze' num2str(sessionNumber)];


baseFile = fullfile(cluFilesDir, sessionName2, sessionName2);


storageDir = fullfile(outputDir,  sessionName2);
if ~exist(storageDir, 'dir'); mkdir(storageDir); end



VarList = {'spikes','behavior','position','speed','basics','ripple'};


for var = 1 : length(VarList)
    load([datasetDir '/wake-' VarList{var} '.mat'])
end


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['position.' sessionName]);
speed    = eval(['speed.' sessionName]);
basics   = eval(['basics.' sessionName]);

% ripple   = eval(['ripple.' sessionName]);
% rippleEvents = ripple.time;


%% SESSION INFORMATION

% load spikes and fileInfo
folderName = fullfile(storageDir,  'spikes');
fileName   = fullfile(folderName, [sessionName2 '.spikes.mat']);

if exist(fileName, 'file')

    load(fileName, 'spikes', 'spikes_pyr', 'fileInfo') % , 'clusterQuality'

else
    mkdir(folderName)

    fileInfo.sessionName   = sessionName2;
    fileInfo.animal        = sessionName2(1:strfind(sessionName2, '-')-1);
    fileInfo.xyt           = [];
    fileInfo.linearPos     = []; % linearized position and lap information=
    fileInfo.speed         = [];
    
    
    if exist([baseFile '.xml'], 'file')
        Par = LoadXml(baseFile); % Load the recording configurations from the xml file

        fileInfo.Fs            = Par.SampleRate; 
        fileInfo.lfpSampleRate = Par.lfpSampleRate;
        fileInfo.nCh           = Par.nChannels;
    else
        fileInfo.Fs            = basics.SampleRate; 
        fileInfo.lfpSampleRate = basics.lfpSampleRate;
        fileInfo.nCh           = basics.nChannels;
    end

    fileInfo.timeUnit      = 1e6; % converted sampling rate

    %     fileInfo.badChannels   = badChannels;

    fileInfo.pix2cm        = 0.3861;


    %%% behavior 

    if strcmp(rat, 'Kevin') % there are two wake periods, we need to concatenate them
       behavior.time(2,2) = behavior.time(3,2);
       behavior.time(3,1) = behavior.time(4,1);
       behavior.time(3,2) = behavior.time(4,2);
       behavior.time(4,:) = [];
    end

    fileInfo.tbegin        = behavior.time(1,1); 
    fileInfo.tend          = behavior.time(3,2);
    
    %%%% should change this section
    
    behavior.time(2,1) = behavior.time(2,1) + 1e6; % slight modification to make sure this time period is only limited to the maze (1e6 = 1 sec)
    behavior.time(2,2) = behavior.time(2,2) - 1e6;
    
    %%%%%%
    
    
    behavior.time = (behavior.time - fileInfo.tbegin)/fileInfo.timeUnit;
    
    %%%% should change this
    fileInfo.behavior.time = behavior.time; 
    %%%%%%


    % brain states

    brainStates = behavior.list;
    brainStates(:, 1:2) = (brainStates(:, 1:2) - fileInfo.tbegin)/fileInfo.timeUnit;
    fileInfo.brainStates = brainStates;

    
    % ripples
    if ~isfield(fileInfo, 'RippleChannels')
        best_channels = BigRippleChannels(baseFile, fileInfo);
        fileInfo.RippleChannels = best_channels; % check if there are no bad channels among the best_channels 
        
        save(fileName, 'fileInfo') % the value of fileInfo will be updated later
    end
    
    
    
    %%% position info and preprocessings

    xpos = position.x * fileInfo.pix2cm; % convert to centimeters
    ypos = position.y * fileInfo.pix2cm;
    tpos = (position.t - fileInfo.tbegin)/fileInfo.timeUnit;

    fileInfo.xyt = [xpos; ypos; tpos]'; 
    % position_samplingRate = 1/(tpos(2) - tpos(1));

    % speed = calspeed(fileInfo, position_samplingRate); % animal's velocity

    speed.t = (speed.t - fileInfo.tbegin)/fileInfo.timeUnit;
    fileInfo.speed = speed;


    %% LINEARIZE POSTION AND CALCULATE LAPS

    mazeShapes    = mazeShape.(rat);
    currMazeShape = mazeShapes{sessionNumber};

    linearPos = linearizePosition(fileInfo, behavior, currMazeShape); % click on the middle of the platforms

    fileInfo.linearPos(:, 1) = linearPos; 
    fileInfo.linearPos(:, 2) = fileInfo.xyt(:, 3);



    % % updated

    direction = 'bi';
    lapsStruct = calculateLapTimings(fileInfo, direction, storageDir); 

    if length(lapsStruct.RL) > length(lapsStruct.LR)
       lapsStruct.RL(1,:) = [];
    end

    totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
    laps = zeros(totNumLaps, 2);
    laps(1:2:totNumLaps, :)  = lapsStruct.LR;
    laps(2:2:totNumLaps, :)  = lapsStruct.RL;

    laps(:, 3) = 1:size(laps, 1); 


    fileInfo.linearPos(:, 3) = zeros(size(fileInfo.linearPos(:, 1))); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

    for ii = 1: length(laps)
       idx =  find(fileInfo.linearPos(:, 2) > laps(ii, 1) & fileInfo.linearPos(:, 2) < laps(ii, 2));
       fileInfo.linearPos(idx, 3) = laps(ii, 3);         
    end


    runSpeedThresh = 10; % cm/s


    %% calculating position, velocity, etc at each spike (these will be needed for calculation of place fields)


    nUnits = numel(spikes);

    for unit = 1:nUnits

        spikes(unit).time = transpose(spikes(unit).time - fileInfo.tbegin)/fileInfo.timeUnit;
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

    clusterQuality = [];
    % try
    %     clusterQuality = calClusterQuality(baseFile);
    % catch
    % end

    %% Spatial tuning of the units

    spikes_pyr = spikes([spikes.quality] <= 3); % limiting to only pyramidal units
    fileInfo.okUnits = 1:numel(spikes_pyr); % all pyramidal units regardless of their stability will be considered in Bayesian decoding



    close all

    subfolder = fullfile(storageDir, 'PlaceFields');

    if ~exist(subfolder, 'dir')
        mkdir(subfolder)
    end


    %%% 1D spatial tuning: using linearized position

    posBinSize = 2;

    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'LR', posBinSize, fileInfo, subfolder);
    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'RL', posBinSize, fileInfo, subfolder);

    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);


    close all
    
    
    fileInfo.tend   = (fileInfo.tend - fileInfo.tbegin)/fileInfo.timeUnit;
    fileInfo.tbegin = 0;

    fileName = fullfile(folderName, [fileInfo.sessionName '.spikes.mat']);
    save(fileName, 'spikes', 'spikes_pyr', 'clusterQuality', 'fileInfo', '-v7.3')

end



%% DETECTING POPULATION BURST EVENTS (PBEs)

folderName = fullfile(storageDir,  'PopulationBurstEvents');
overwrite = 0;

if exist(folderName, 'dir') && overwrite == 0
    
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo.mat']), 'PBEInfo_Bayesian', 'sdat');
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo.mat']), 'PBEInfo')
    
else
        
    mkdir(folderName)
    
    MUA.time = []; for unit = 1:numel(spikes); MUA.time = [MUA.time; spikes(unit).time]; end
    MUA.time = sort(MUA.time, 'ascend');
    
% %     PYR and MUA for detection of PBEs
% %     MUA
%     nShanks = length(dir([baseFile '.res.*']));
% 
%     
%     MUA.time = [];
%     for shank = 1:nShanks
% 
%         resContent  = load([baseFile '.res.' num2str(shank)]);
% 
%         cluContent  = load([baseFile '.clu.' num2str(shank)]);
%         clusterIDs  = cluContent(2:end); % the first entry is the total number of clusters
% 
%         spikes2Add = resContent(~ismember(clusterIDs, 0)); % to exclude the noise cluster 
% 
%         MUA.time = [MUA.time; spikes2Add]; % the first is the total number of clusters detected for the shank   
%     end
%     MUA.time = sort(MUA.time, 'ascend');
%     MUA.time = MUA.time./fileInfo.Fs;

    
    
    % detection parameters
    
    pbeDetectParams.time_resolution        = 0.001;
    pbeDetectParams.threshold              = 2;
    pbeDetectParams.exclude                = [];

    pbeDetectParams.smoothingSigma         = 0.02;
    
    pbeDetectParams.minDur                 = 0.04; % number of bins
    pbeDetectParams.maxDur                 = 0.70;
    
    pbeDetectParams.maxSilence             = 0.05;
    
    pbeDetectParams.maxVelocity            = 10; 
    
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
    
    
    [PBEInfo, PBEInfo_Bayesian, idx_Bayesian] = genEventStructure_Hiro(PBEs_splitted, spikes_pyr, sdat, pbeDetectParams, fileInfo, folderName);
    
end


% %% add new channels corresponding to theta/ripple/speed to the eeg file
% 
% nCh = fileInfo.nCh;
% 
% 
% EEGdata = readmulti([baseFile '.eeg'], nCh, 1); % load .eeg trace
% totalT  = length(EEGdata)/fileInfo.lfpSampleRate;




