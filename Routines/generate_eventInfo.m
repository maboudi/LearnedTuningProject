function generate_eventInfo(sessionNumber)



parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

addpath(genpath(fullfile(parentDir, '/ReplayPreplayAnalyses')))

datasetDir = fullfile(parentDir, '/Grosmark_reclustered');
cluFilesDir = '/data/datasets/GrossmarkData_BapunsClusters';


stateDetectionResultsDir = fullfile(parentDir, 'StateDetectionResults');

outputDir = fullfile(parentDir, 'sessions_calculated_PBEinfo2');
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
        
sessionName = sessionNames{sessionNumber};

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


% loading slow wave power/slope 
filename = sprintf([sessionName '.SleepScoreMetrics.mat']);
load(fullfile(stateDetectionResultsDir, sessionName, filename), 'SleepScoreMetrics')


slowWave    = SleepScoreMetrics.broadbandSlowWave;
swthresh    = SleepScoreMetrics.histsandthreshs.swthresh;
sw_timePnts = SleepScoreMetrics.t_clus;


% loading emg 
emg = SleepScoreMetrics.EMG;
load(fullfile(stateDetectionResultsDir, sessionName, sprintf([sessionName '.emgThresh.mat'])), 'emgThresh');


% bad channels to exclude for lpf event detections (e.g., ripples)

filename = sprintf([sessionName '.sessionInfo.mat']);
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
    fileInfo.xyt2          = []; % linearized position and lap information=
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
    
    
%     spikes = spikes([spikes.quality] <= 3); % limiting to only pyramidal units

%     fileInfo.stablePyr = find([spikes.StablePrePostFiner] == 1); % pyramidal units that are stable
%     fileInfo.okUnits   = 1:numel(spikes);



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


    fileInfo.xyt2(:, 1) = linearPos; 
    fileInfo.xyt2(:, 2) = fileInfo.xyt(:, 3);


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
    

    fileInfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

    for ii = 1: length(laps)
       idx =  fileInfo.xyt2(:, 2) > laps(ii, 1) & fileInfo.xyt2(:, 2) < laps(ii, 2);
       fileInfo.xyt2(idx, 3) = laps(ii, 3);
    end

    runSpeedThresh = 10; % cm/s


    %% calculating position, velocity, etc at each spike (these will be needed for calculation of place fields)


    nUnits = numel(spikes);

    for unit = 1:nUnits

        spikeTimes = spikes(unit).time;
        spikes(unit).x = interp1(fileInfo.xyt(:, 3), fileInfo.xyt(:, 1), spikeTimes);
        spikes(unit).y = interp1(fileInfo.xyt(:, 3), fileInfo.xyt(:, 2), spikeTimes);
        spikes(unit).linearPos = interp1(fileInfo.xyt2(:, 2), fileInfo.xyt2(:, 1), spikeTimes);

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
fileName = fullfile(folderName, [fileInfo.sessionName '.rippleInfo.mat']);

if exist(fileName, 'file')
    load(fileName, 'detectParams', 'rippleLFP_adjust', 'ripplePower_summed_adjust', 'rippleInfo', 'tq')
end


if ~exist('rippleInfo', 'var')
    

    threshSD = 2;
    [ripplePeriods, rippleLFP, ripplePower_summed] = RippleDetect(baseFile, fileInfo, threshSD); % the results in seconds


    % downsample the ripple band lfp or power to 1kHz (1ms sampling period) % maybe we don't
    % need to downsample at all?

    ts    = 1/fileInfo.lfpSampleRate;
    
    tend  = length(rippleLFP)/fileInfo.lfpSampleRate;
    tPnts = ts:ts:fileInfo.tend;

    tq    = 1e-3:1e-3:fileInfo.tend;


    rippleLFP_adjust = zeros(numel(tq), size(rippleLFP, 2));
    for ich = 1:size(rippleLFP, 2)
        rippleLFP_adjust(:, ich) = interp1(tPnts, rippleLFP(:, ich), tq);
    end

    ripplePower_summed_adjust = interp1(tPnts, ripplePower_summed, tq);
    
    
    
    
    %% GENERATE RIPPLE INFO
    
    mkdir(folderName);
    
    nRipples = size(ripplePeriods, 1);

    rippleInfo = struct('sessionName', [], 'rippleno', [], ...
                        'startT', [], 'endT', [], 'peakT', [], 'duration', [], 'peakRippleA', [], 'peakMUA', [], ...
                        'rippleWaveForm', [], 'epoch', []);


    for rpl = 1:nRipples    

        rippleInfo(rpl).sessionName = sessionName;
        rippleInfo(rpl).rippleno    = rpl;

        rippleInfo(rpl).startT   = ripplePeriods(rpl, 1);
        rippleInfo(rpl).endT     = ripplePeriods(rpl, 2);
        rippleInfo(rpl).peakT    = ripplePeriods(rpl, 3);
        rippleInfo(rpl).duration = rippleInfo(rpl).endT - rippleInfo(rpl).startT;

        rippleInfo(rpl).peakRippleA = ripplePeriods(rpl, 4);

%         rippleInfo(rpl).peakMUA  = max(sdat(floor(rippleInfo(rpl).startT*1e3):floor(rippleInfo(rpl).endT*1e3)));


        try
            rippleInfo(rpl).rippleWaveForm = rippleLFP_adjust(floor(rippleInfo(rpl).startT*1e3)-500:floor(rippleInfo(rpl).endT*1e3)+500, 1);
        catch
            rippleInfo(rpl).rippleWaveForm = [];
        end

        if rippleInfo(rpl).startT > behavior.time(1,1) && rippleInfo(rpl).startT < behavior.time(1,2)
           rippleInfo(rpl).epoch = 'pre';
        elseif rippleInfo(rpl).startT > behavior.time(2,1) && rippleInfo(rpl).startT < behavior.time(2,2)
           rippleInfo(rpl).epoch = 'run';
        elseif rippleInfo(rpl).startT > behavior.time(3,1) && rippleInfo(rpl).startT < behavior.time(3,2)
           rippleInfo(rpl).epoch = 'post';
        end

    end


    fileName = fullfile(folderName, [fileInfo.sessionName '.rpl.evt']);
    MakeEvtFile(ripplePeriods(:, 1:3), fileName, {'beg', 'end', 'peak'}, 1, 1)


    detectParams.threshold = threshSD;
    detectParams.method    = 'summedChan';


    fileName = fullfile(folderName, [fileInfo.sessionName '.rippleInfo.mat']);
    save(fileName, 'detectParams', 'rippleLFP_adjust', 'ripplePower_summed_adjust', 'tq', 'rippleInfo','-v7.3')   % , 'ripplePower_adjust'

end


%% DETECTING POPULATION BURST EVENTS (PBEs)

folderName = fullfile(storageDir,  'PopulationBurstEvents');
overwrite = 1;

if exist(folderName, 'dir') && overwrite == 0
    
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo_MUA.mat']), 'PBEInfo_MUA', 'detectParams')
    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo_PYR.mat']), 'PBEInfo_PYR_accepted', 'sdat');

    load(fullfile(folderName, [fileInfo.sessionName '.PBEInfo_PYR_all.mat']), 'PBEInfo_PYR_all')
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

        spikes2Add = resContent(~ismember(clusterIDs, 0)); % to exclude the noise cluster and maybe the mutiunit

        MUA.time = [MUA.time; spikes2Add]; % the first is the total number of clusters detected for the shank   
    end
    MUA.time = sort(MUA.time, 'ascend');
    MUA.time = MUA.time./fileInfo.Fs;

    
    % PYR
    PYR.time = [];

    for unit = 1:numel(spikes_pyr)
        PYR.time = [PYR.time; spikes_pyr(unit).time];
    end
    PYR.time = sort(PYR.time, 'ascend');
    
    


    % detection parameters
    
    clear detectParams
    
%     detectParams.threshold              = 2;
    
    detectParams.threshold.MUA          = 2;
    detectParams.threshold.PYR          = 1;
    
    detectParams.MUAorPyr               = '';
    detectParams.time_resolution        = 0.001;
    detectParams.exclude                = [];
    detectParams.smoothingSigma.MUA     = 0.01;
    detectParams.smoothingSigma.PYR     = 0.01;
    detectParams.minDur                 = 0.04; % number of bins
    detectParams.maxDur                 = 0.60;
    detectParams.minNumTimeBins         = 3;
    detectParams.maxNumTimeBins         = 30;
    detectParams.minFiringUnits         = max(5, 0.1*numel(fileInfo.okUnits)); % okUnits are all pyramidal units regardless of their stability
    detectParams.maxVelocity            = 10; 
    detectParams.thRatioThresh          = 1;
    detectParams.rippleAThresh_low      = 1;
    detectParams.rippleAThresh_high     = 3;


    
    % detecting primary PBEs, just the boundaries of potential PBEs before
    % applying any criteria

    [primaryPBEs, sdat] = PBPeriods_MUA_PYR2(MUA, PYR, detectParams, fileInfo); % already filtering based on the duration in this step

    
    %%
    spikes_pyr_pooled = cell(numel(spikes_pyr), 1);
    for iUnit = 1:numel(spikes_pyr)
        spikes_pyr_pooled{iUnit} = spikes_pyr(iUnit).time';
    end
    spikes_pyr_pooled = cell2mat(spikes_pyr_pooled');
    spikes_pyr_pooled = sort(spikes_pyr_pooled, 'ascend');
    
    
    
    PBEs_splitted = [];
    for ipbe = 1:size(primaryPBEs.MUA, 1)
        
        currSpikes = spikes_pyr_pooled(spikes_pyr_pooled >= primaryPBEs.MUA(ipbe, 1) & spikes_pyr_pooled <= primaryPBEs.MUA(ipbe, 2));
        cutIdx     = find([0 diff(currSpikes)] > 0.03); % indices of the spikes that occurred longer than 30 ms after the previous spike
        
        if isempty(cutIdx)
            % the original PBE is passed
            
            PBEs_splitted = [PBEs_splitted; primaryPBEs.MUA(ipbe, :)];
        else
            % the PBE is split to multiple
            nseg = numel(cutIdx)+1;
            
            split_starts = [primaryPBEs.MUA(ipbe, 1) currSpikes(cutIdx)];
            split_ends   = [currSpikes(cutIdx-1) primaryPBEs.MUA(ipbe, 2)];
            
            for isplit = 1:nseg
                
                [max_MUA, peakIdx] = max(sdat(floor(split_starts(isplit)*1e3) : floor(split_ends(isplit)*1e3)));
                
                duration =  split_ends(isplit)- split_starts(isplit);
                if duration >= detectParams.minDur
                    PBEs_splitted = [PBEs_splitted; [split_starts(isplit) split_ends(isplit) split_starts(isplit)+peakIdx*1e-3 max_MUA]];
                end
                
            end
        end
    end
    
    
    
    
    %% MUA-based detected PBEs
    % This set of detected PBEs is going to be used for assessing how unit participation evolves over sleep


    nPBEs = size(primaryPBEs.MUA, 1);

    PBEInfo = struct('sessionName', [], 'PBEno', [], ..., 
                    'startT', [], 'endT', [], 'peakT', [], ...
                    'peakMUA', [], 'MUAwave', [], ...
                    'duration', [], ...
                    'peakRippleA', [], 'rippleWave', [], ...
                    'thetaRatio', [], 'SWA', [], 'emg', [], ...
                    'epoch', [], 'brainState', []);


    for ipbe = 1:nPBEs

        PBEInfo(ipbe).sessionName = sessionName;
    %     PBEInfo(ipbe).PBEno       = ipbe;

        PBEInfo(ipbe).startT      = primaryPBEs.MUA(ipbe, 1);
        PBEInfo(ipbe).endT        = primaryPBEs.MUA(ipbe,2);
        PBEInfo(ipbe).peakT       = primaryPBEs.MUA(ipbe,3);
        PBEInfo(ipbe).duration    = PBEInfo(ipbe).endT - PBEInfo(ipbe).startT;
        
        PBEInfo(ipbe).peakMUA     = primaryPBEs.MUA(ipbe,4);

        
%         [PBEInfo(ipbe).peakRippleA, maxChIdx] = max(max(ripplePower_summed_adjust(floor(PBEInfo(ipbe).startT*1e3): floor(PBEInfo(ipbe).endT*1e3), :), [], 1));

        try
            PBEInfo(ipbe).rippleWave = rippleLFP_adjust(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500), 1);
        catch
            PBEInfo(ipbe).rippleWave = [];
        end
        PBEInfo(ipbe).peakRippleA = max(ripplePower_summed_adjust(floor(PBEInfo(ipbe).startT*1e3): floor(PBEInfo(ipbe).endT*1e3)));
        
        
        PBEInfo(ipbe).thetaRatio  = interp1(timePnts, theratio, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).SWA         = interp1(sw_timePnts, slowWave, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).emg         = interp1(sw_timePnts, emg, PBEInfo(ipbe).peakT);


        try
            PBEInfo(ipbe).MUAwave = sdat.MUA(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500));
        catch
            PBEInfo(ipbe).MUAwave = [];
        end

        

        PBEInfo(ipbe).velocity    = interp1(fileInfo.speed.t, fileInfo.speed.v, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).velocity(isnan(PBEInfo(ipbe).velocity)) = -1;


        if PBEInfo(ipbe).startT > behavior.time(1,1) && PBEInfo(ipbe).startT < behavior.time(1,2)
           PBEInfo(ipbe).epoch = 'pre';
        elseif PBEInfo(ipbe).startT > behavior.time(2,1) && PBEInfo(ipbe).startT < behavior.time(2,2)
           PBEInfo(ipbe).epoch = 'run';
        elseif PBEInfo(ipbe).startT > behavior.time(3,1) && PBEInfo(ipbe).startT < behavior.time(3,2)
           PBEInfo(ipbe).epoch = 'post';
        end


        if PBEInfo(ipbe).SWA > swthresh && ismember(PBEInfo(ipbe).epoch, {'pre'; 'post'})
            PBEInfo(ipbe).brainState = 'NREM';
        elseif PBEInfo(ipbe).thetaRatio < detectParams.thRatioThresh || (PBEInfo(ipbe).peakRippleA > detectParams.rippleAThresh_high && PBEInfo(ipbe).peakRippleA < 20) % hte ceiling threshold on ripple is to exclude the artifacts
            PBEInfo(ipbe).brainState = 'QW';
        elseif PBEInfo(ipbe).emg < emgThresh
            PBEInfo(ipbe).brainState = 'REM';
        else
            PBEInfo(ipbe).brainState = 'WAKE';
        end
        
        if ipbe == 1
            fprintf('\n dividing the PBEs into 20 ms time bins ..')
        end

        if mod(ipbe, 1000) == 0
            fprintf('.')
        end
        

    end


    PBEInfo_MUA = PBEInfo(ismember({PBEInfo.brainState}, {'NREM'; 'QW'}) & ...
                            [PBEInfo.velocity] < detectParams.maxVelocity & ...
                            [PBEInfo.peakRippleA] > detectParams.rippleAThresh_low);


    for ipbe = 1:numel(PBEInfo_MUA)
        PBEInfo_MUA(ipbe).PBEno = ipbe;
    end

    fileName = fullfile(folderName, [fileInfo.sessionName '.PBEInfo_MUA2.mat']);
    save(fileName, 'PBEInfo_MUA', 'detectParams', '-v7.3')
    
    
    % generate event files of vidualizing on Neuroscope 
    nPBEs = numel(PBEInfo_MUA);

    PBEInfo_MUA_neurscope = zeros(nPBEs, 3);
    PBEInfo_MUA_neurscope(:, 1) = [PBEInfo_MUA.startT];
    PBEInfo_MUA_neurscope(:, 2) = [PBEInfo_MUA.endT];
    PBEInfo_MUA_neurscope(:, 3) = [PBEInfo_MUA.peakT];
    
    
    fileName = fullfile(folderName, [fileInfo.sessionName '.mu2.evt']);
    MakeEvtFile(PBEInfo_MUA_neurscope(:, 1:3), fileName, {'beg', 'end', 'peak'}, 1, 1)
    
    
    
    %% PYR-based detected PBEs
    % This set of detected PBEs is going to be used mostly for Bayesian decoding. We save both vesrion before and after appplying some filters based on duration, participation, etc 

    clear PBEInfo

    nPBEs = size(primaryPBEs.PYR, 1);
    PBEInfo = struct('sessionName', [], 'PBEno', [], ...
                        'startT', [], 'endT', [], 'peakT', [], ...
                        'startT_adj', [], 'endT_adj', [], 'peakT_adj', [], ...
                        'peakMUA', [], 'MUAwave', [], ...
                        'nFiringUnits', [], 'duration', [], 'duration_adj', [], ...
                        'peakRippleA', [], 'rippleWave', [], ...
                        'thetaRatio', [], 'SWA', [], 'emg', [], ...
                        'epoch', [], 'brainState', [], ...
                        'fr_1msbin', [], 'fr_20msbin', []);



    binDur = 0.02; % in sec

    for ipbe = 1:nPBEs

        PBE_primary = primaryPBEs.PYR(ipbe, :);

        [binnedPBE, PBE_adj, nFiringUnits, PBElength] = finalBinningResult(PBE_primary, spikes_pyr, binDur, fileInfo); 


        PBEInfo(ipbe).sessionName = sessionName;
    %     PBEInfo(ipbe).PBEno       = ipbe;

        PBEInfo(ipbe).startT      = PBE_primary(1);
        PBEInfo(ipbe).endT        = PBE_primary(2);
        PBEInfo(ipbe).peakT       = PBE_primary(3);

        PBEInfo(ipbe).startT_adj  = PBE_adj(1);
        PBEInfo(ipbe).endT_adj    = PBE_adj(2);
        PBEInfo(ipbe).peakT_adj   = PBE_adj(3);


        PBEInfo(ipbe).duration     = PBE_primary(2) - PBE_primary(1);
        PBEInfo(ipbe).duration_adj = PBElength; 

        PBEInfo(ipbe).nFiringUnits = nFiringUnits;

        
%         [PBEInfo(ipbe).peakRippleA, maxChIdx] = max(max(ripplePower_adjust(floor(PBEInfo(ipbe).startT*1e3): floor(PBEInfo(ipbe).endT*1e3), :), [], 1));

        try
            PBEInfo(ipbe).rippleWave = rippleLFP_adjust(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500), 1);
        catch
            PBEInfo(ipbe).rippleWave = [];
        end
        
        PBEInfo(ipbe).peakRippleA = max(ripplePower_summed_adjust(floor(PBE_primary(1)*1e3): floor(PBE_primary(2)*1e3)));

        

        PBEInfo(ipbe).peakMUA = PBE_adj(4);
        try
            PBEInfo(ipbe).MUAwave = sdat.MUA(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500));
        catch
            PBEInfo(ipbe).MUAwave = [];
        end


        PBEInfo(ipbe).thetaRatio  = interp1(timePnts, theratio, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).SWA         = interp1(sw_timePnts, slowWave, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).emg         = interp1(sw_timePnts, emg, PBEInfo(ipbe).peakT);


        PBEInfo(ipbe).linearPos   = interp1(fileInfo.xyt2(:, 2), fileInfo.xyt2(:, 1), PBEInfo(ipbe).peakT);

        PBEInfo(ipbe).velocity    = interp1(fileInfo.speed.t, fileInfo.speed.v, PBEInfo(ipbe).peakT);
        PBEInfo(ipbe).velocity(isnan(PBEInfo(ipbe).velocity)) = -1;


        PBEInfo(ipbe).fr_1msbin   = binnedPBE{1};
        PBEInfo(ipbe).fr_20msbin  = binnedPBE{2};


        if PBEInfo(ipbe).startT > behavior.time(1,1) && PBEInfo(ipbe).startT < behavior.time(1,2)
           PBEInfo(ipbe).epoch = 'pre';
        elseif PBEInfo(ipbe).startT > behavior.time(2,1) && PBEInfo(ipbe).startT < behavior.time(2,2)
           PBEInfo(ipbe).epoch = 'run';
        elseif PBEInfo(ipbe).startT > behavior.time(3,1) && PBEInfo(ipbe).startT < behavior.time(3,2)
           PBEInfo(ipbe).epoch = 'post';
        end


        % estimating brain states of each PBE
        if PBEInfo(ipbe).SWA > swthresh && ismember(PBEInfo(ipbe).epoch, {'pre'; 'post'})
            PBEInfo(ipbe).brainState = 'NREM';
        elseif PBEInfo(ipbe).thetaRatio < detectParams.thRatioThresh || (PBEInfo(ipbe).peakRippleA > detectParams.rippleAThresh_high && PBEInfo(ipbe).peakRippleA < 20)
            PBEInfo(ipbe).brainState = 'QW';
        elseif PBEInfo(ipbe).emg < emgThresh
            PBEInfo(ipbe).brainState = 'REM';
        else
            PBEInfo(ipbe).brainState = 'WAKE';
        end


        if ipbe == 1
            fprintf('\n dividing the PBEs into 20 ms time bins ..')
        end

        if mod(ipbe, 1000) == 0
            fprintf('.')
        end

    end
    
%%

    PBEInfo_PYR_all = PBEInfo(ismember({PBEInfo.brainState}, {'NREM'; 'QW'}) & ...
                            [PBEInfo.velocity] < detectParams.maxVelocity & ...
                            [PBEInfo.peakRippleA] > detectParams.rippleAThresh_low);
                        
                        
    for ipbe = 1:numel(PBEInfo_PYR_all)
        PBEInfo_PYR_all(ipbe).PBEno = ipbe;
    end
    
    
    % A smaller subset of PBEs with stricter criteria, on which we are going to
    % apply Bayesian decoding
    acceptedIdx = [PBEInfo_PYR_all.duration_adj] >= detectParams.minNumTimeBins & ...
                    [PBEInfo_PYR_all.duration_adj] <= detectParams.maxNumTimeBins & ...
                    [PBEInfo_PYR_all.nFiringUnits] >= detectParams.minFiringUnits;

    PBEInfo_PYR_accepted = PBEInfo_PYR_all(acceptedIdx);



    % save a .mat file contatining PBEs that passed the criteria
    fileName = fullfile(folderName, [fileInfo.sessionName '.PBEInfo_PYR2.mat']);
    save(fileName, 'PBEInfo_PYR_accepted', 'detectParams', 'sdat', 'tq', '-v7.3')
    
    
     % store all the PBEs before doing further filterings based on duration
    % or unit participation
    fileName = fullfile(folderName, [fileInfo.sessionName '.PBEInfo_PYR_all2.mat']);
    save(fileName, 'PBEInfo_PYR_all', 'acceptedIdx', '-v7.3')

    
    
    
    % generate event files of vidualizing on Neuroscope (all the PBEs not only the accepted ones for Bayesian decoding)
    nPBEs = numel(PBEInfo_PYR_all);

    PBEInfo_PYR_all_neurscope = zeros(nPBEs, 3);
    PBEInfo_PYR_all_neurscope(:, 1) = [PBEInfo_PYR_all.startT];
    PBEInfo_PYR_all_neurscope(:, 2) = [PBEInfo_PYR_all.endT];
    PBEInfo_PYR_all_neurscope(:, 3) = [PBEInfo_PYR_all.peakT];


    % make an .evt file for the PBEs
    fileName = fullfile(folderName, [fileInfo.sessionName '.py2.evt']);
    MakeEvtFile(PBEInfo_PYR_all_neurscope(:, 1:3), fileName, {'beg', 'end', 'peak'}, 1, 1)
   
    

end



%% adding MUA at the time of the ripple to the ripple Info

if ~isfield(rippleInfo, 'peakMUA')

    for rpl = 1:nRipples 

        try
            rippleInfo(rpl).peakMUA  = max(sdat(floor(rippleInfo(rpl).startT*1e3):floor(rippleInfo(rpl).endT*1e3)));
        catch
            rippleInfo(rpl).peakMUA  = [];
        end

    end

    folderName = fullfile(storageDir,  'ripples');
    fileName = fullfile(folderName, [fileInfo.sessionName '.rippleInfo.mat']);
    save(fileName, 'rippleInfo','-v7.3')
end




%% adding variables to the EEG file

nCh = fileInfo.nCh;


EEGdata = readmulti([baseFile '.eeg'], nCh, 1); % load .eeg trace
totalT  = length(EEGdata)/fileInfo.lfpSampleRate;


samplingDur  = 1/fileInfo.lfpSampleRate;
samplingPnts = samplingDur:samplingDur:totalT;

sdat_mua_2add = interp1((1:length(sdat.MUA))*1e-3, sdat.MUA, samplingPnts);
sdat_pyr_2add = interp1((1:length(sdat.PYR))*1e-3, sdat.PYR, samplingPnts);
speed_2add    = interp1(fileInfo.speed.t, fileInfo.speed.v, samplingPnts);
ripple_2add   = interp1(tq, ripplePower_summed_adjust, samplingPnts);
theratio_2add = interp1(timePnts, theratio, samplingPnts);


sdat_mua_2add(isnan(sdat_mua_2add)) = 0;
sdat_mua_2add(sdat_mua_2add < 0) = 0;
sdat_mua_2add(sdat_mua_2add > 2) = 2;
sdat_mua_2add = (2*(sdat_mua_2add-min(sdat_mua_2add))/(max(sdat_mua_2add)-min(sdat_mua_2add))-1)*2^10;


sdat_pyr_2add(isnan(sdat_pyr_2add)) = 0;
sdat_pyr_2add(sdat_pyr_2add < 0) = 0;
sdat_pyr_2add(sdat_pyr_2add > 2) = 2;
sdat_pyr_2add = (2*(sdat_pyr_2add-min(sdat_pyr_2add))/(max(sdat_pyr_2add)-min(sdat_pyr_2add))-1)*2^10;


speed_2add(isnan(speed_2add)) = 0;
speed_2add(speed_2add < 10) = 0;
speed_2add = (2*(speed_2add-min(speed_2add))/(max(speed_2add)-min(speed_2add))-1)*2^11;


ripple_2add(isnan(ripple_2add)) = 0;
ripple_2add(ripple_2add > 3) = 3;
ripple_2add(ripple_2add < 0) = 0;
ripple_2add = (2*(ripple_2add - min(ripple_2add))/(max(ripple_2add)-min(ripple_2add))-1)*2^10;


theratio_2add(isnan(theratio_2add)) = 0;
theratio_2add(theratio_2add < 1) = 0;
theratio_2add = (2*(theratio_2add - min(theratio_2add))/(max(theratio_2add)-min(theratio_2add))-1)*2^11;



inputFileHandle  = fopen([baseFile,'.eeg']);
outputFileHandle = fopen([baseFile,'-2.eeg'],'w');
doneFrame = 0;
bufferSize = 4096;
while ~feof(inputFileHandle)
    data = [fread(inputFileHandle,[nCh,bufferSize],'int16')]';
    
    for frame=1:size(data,1)
        fwrite(outputFileHandle,[data(frame,:), sdat_mua_2add(frame+doneFrame), sdat_pyr_2add(frame+doneFrame), speed_2add(frame+doneFrame), ripple_2add(frame+doneFrame), theratio_2add(frame+doneFrame)]','int16');
    end
    doneFrame=doneFrame+frame;

end

fclose(inputFileHandle);
fclose(outputFileHandle);



%% Bayesian decoding


% BayesianReplayDetection_GrosmarkReclu_nov2020




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
