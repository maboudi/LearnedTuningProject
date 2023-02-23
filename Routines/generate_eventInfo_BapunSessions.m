function generate_eventInfo_BapunSessions(sessionNumber)


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

addpath(genpath(fullfile(parentDir, '/ReplayPreplayAnalyses')))

datasetDir = fullfile(parentDir, '/Bapun_NSD_datasets');


stateDetectionResultsDir = fullfile(parentDir, 'StateDetectionResults');

outputDir = fullfile(parentDir, 'sessions_calculated_PBEinfo2');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

sessionNames = {'RatN_Day2_2019-10-11_03-58-54'; 'RatS-Day2-2020-11-27_10-22-29'};



% MAZE SHAPES
mazeLimits = [-inf inf -inf 128; ... % RatN
              -inf inf -inf inf]; % Rat S
mazeShapes = {'L-shape'; 'L-shape'};


% CURRENT SESSION        
sessionName = sessionNames{sessionNumber};




%% loading session data

baseFile = fullfile(datasetDir, sessionName, sessionName);

varList = {'spike'; 'epochs'; 'positions'};
for var = 1 : length(varList)
    load([baseFile '.' varList{var} '.mat'])
end


storageDir = fullfile(outputDir,  sessionName);
if ~exist(storageDir, 'dir'); mkdir(storageDir); end


try
    behavior.time = double([pre; maze; post]);
catch
    behavior.time = double([pre; maze1; post]);
end


% give a warning if there is a gap between the end of one epoch and start
% of the next epoch

if behavior.time(2,1) ~= behavior.time(1,2) || behavior.time(3,1) ~= behavior.time(2,2) 
    warning('There are gaps between the epochs')
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
sw_timePnts = SleepScoreMetrics.t_clus;

try 
  load(fullfile(stateDetectionResultsDir, sessionName, sprintf([sessionName '.swthresh.mat'])), 'swthresh'); 
catch
    swthresh = SleepScoreMetrics.histsandthreshs.swthresh;
end


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
fileName   = fullfile(folderName, [sessionName '.spikes.mat']);

if exist(fileName, 'file')

    load(fileName, 'spikes', 'spikes_pyr', 'fileInfo') % , 'clusterQuality'

else
     mkdir(folderName); 

    Par = LoadXml(baseFile); % Load the recording configurations from the xml file

    fileInfo.sessionName   = sessionName;
    fileInfo.animal        = sessionName(1:strfind(sessionName, '_')-1);
    fileInfo.xyt           = [];
    fileInfo.xyt2          = []; % linearized position and lap information
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
    
    if sessionNumber == 1 % rat N
        position_samplingRate = 120;
    elseif sessionNumber == 2 % rat S
        position_samplingRate = 60;
    end

    y = positions(1,:)';
    x = positions(2,:)';
    time = (0:1/position_samplingRate:(length(positions)-1)/position_samplingRate)';



    fileInfo.xyt = [x y time]; 
      

    speed = calspeed(fileInfo, position_samplingRate); % animal's velocity
    
    fileInfo.speed = speed;
    


    %% LINEARIZE POSTION AND CALCULATE LAPS

    fileInfo.xyt(x > mazeLimits(sessionNumber, 2) | x < mazeLimits(sessionNumber, 1) | y > mazeLimits(sessionNumber, 4) | y < mazeLimits(sessionNumber, 3), :) = [];


    linearPos = linearizePosition2(fileInfo, behavior, mazeShapes{sessionNumber});


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

    spikes = spike; % just to make the variable names the same across different datasets
    nUnits = numel(spikes);


    stabilityPeriod_pre  = [max(behavior.time(1,2)-3*60*60, behavior.time(1,1)) behavior.time(1,2)];
    stabilityPeriod_maze =  behavior.time(2,:);
    stabilityPeriod_post = [behavior.time(3,1) behavior.time(3,1)+3*60*60];

    stb = [stabilityPeriod_pre; stabilityPeriod_maze; stabilityPeriod_post];


    for unit = 1:nUnits

        spikes(unit).time = spikes(unit).t';
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

        spikes(unit).id = [spikes(unit).shank spikes(unit).id];
        spikes(unit).quality = spikes(unit).q;


        %
        spikes(unit).fr_by_epoch = zeros(1,3); 
        for ip = 1:size(stb, 1)
           spikes(unit).fr_by_epoch(ip) = numel(find(spikes(unit).time > stb(ip, 1) & spikes(unit).time <= stb(ip, 2)))/diff(stb(ip, :));
        end

        spikes(unit).fr_pre_post_ratio = spikes(unit).fr_by_epoch(3)/ spikes(unit).fr_by_epoch(1); % the difference between the pre and post epoch next the the maze is enough


        if spikes(unit).fr_pre_post_ratio > 1
            spikes(unit).fr_pre_post_bias = 1; % smaller as a percent of the greater firing rate
            spikes(unit).fr_pre_post_ratio = 1/spikes(unit).fr_pre_post_ratio;
        else
            spikes(unit).fr_pre_post_bias = -1;
        end


        if spikes(unit).fr_pre_post_ratio > 0.3 && spikes(unit).fr_by_epoch(2) > 0.02 && strcmp(spikes(unit).group, 'good')
            spikes(unit).fr_isStable = 1;
        else
            spikes(unit).fr_isStable = 0;
        end

    end

    spikes = rmfield(spikes, {'t'; 'shank'; 'q'; 'sh'; 'Var1'});

    
    %% cluster quality

    % wait to recieve the features files from Bapun 

    % clusterQuality = calClusterQuality(baseFile);



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

    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], fileInfo.speed, runSpeedThresh, 'LR', posBinSize, fileInfo, subfolder);
    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], fileInfo.speed, runSpeedThresh, 'RL', posBinSize, fileInfo, subfolder);

    spikes_pyr = spatialTuning_1D_June21(spikes_pyr, behavior, [], [], fileInfo.speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);


    close all
    
    
    % save spikes and file info
    fileName = fullfile(folderName, [sessionName '.spikes.mat']);
    save(fileName, 'spikes', 'spikes_pyr', 'fileInfo', '-v7.3') % , 'clusterQuality'
    

end



%% generate a clu and res file for interneurons and mua for visualizing population bursts

subfolder = fullfile(storageDir, 'placeFields');
if ~exist(subfolder, 'dir') 
    mkdir(subfolder)
end

int_and_MUA_Idx = find(ismember([spikes.quality], [6 8]));

res_t = [];
clu_t = [];

for unit = 1:numel(int_and_MUA_Idx)
    
    spikeTimes = floor(spikes(int_and_MUA_Idx(unit)).time*fileInfo.Fs);
    cluIdx     = unit*ones(numel(spikeTimes), 1);
    
    res_t = [res_t; spikeTimes];
    clu_t = [clu_t; cluIdx];
    
end

[res_t, ind] = sort(res_t);
clu_t = clu_t(ind)-1;  

Saveres(fullfile(subfolder, [fileInfo.sessionName 'mua.201'  '.res.3']), res_t);
SaveClu(fullfile(subfolder, [fileInfo.sessionName 'mua.201'  '.clu.3']), clu_t);



%% RIPPLE DETECTION

folderName = fullfile(storageDir,  'ripples');
fileName = fullfile(folderName, [fileInfo.sessionName '.rippleInfo.mat']);

if exist(fileName, 'file')
    load(fileName, 'detectParams', 'rippleLFP_adjust', 'ripplePower_summed_adjust', 'rippleInfo', 'tq')
end


if ~exist('rippleInfo', 'var')
    
    % if want to detect ripples invidually for each channel
%     threshSD = 5;
%     [ripplePeriods, rippleLFP, ripplePower, ripplePower_summed] = RippleDetect2(baseFile, fileInfo, threshSD);
    
    threshSD = 2;
    [ripplePeriods, rippleLFP, ripplePower_summed] = RippleDetect(baseFile, fileInfo, threshSD); % the results in seconds
    

    % downsample the ripple band lfp or power to 1kHz (1ms sampling period) % maybe we don't
    % need to downsample at all?
    
    ts    = 1/fileInfo.lfpSampleRate;

    tend  = length(rippleLFP)/fileInfo.lfpSampleRate;
    tPnts = ts:ts:tend;

    tq    = 1e-3:1e-3:tend;

    
    rippleLFP_adjust   = zeros(numel(tq), size(rippleLFP, 2));
%     ripplePower_adjust = zeros(numel(tq), size(ripplePower, 2));
    
    for ich = 1:size(rippleLFP, 2)
        rippleLFP_adjust(:, ich)   = interp1(tPnts, rippleLFP(:, ich), tq);
%         ripplePower_adjust(:, ich) = interp1(tPnts, ripplePower(:, ich), tq);
    end

    ripplePower_summed_adjust = interp1(tPnts, ripplePower_summed, tq);
   
    
    
    %% GENERATE RIPPLE INFO
    
    mkdir(folderName);


    nRipples = size(ripplePeriods, 1);

    rippleInfo = struct('sessionName', [], 'rippleno', [], ...
                        'startT', [], 'endT', [], 'peakT', [], 'duration', [], 'peakRippleA', [], 'peakChIdx', [], 'peakMUA', [], ...
                        'rippleWaveForm', [], 'epoch', []);


    for rpl = 1:nRipples    

        rippleInfo(rpl).sessionName = sessionName;
        rippleInfo(rpl).rippleno    = rpl;

        rippleInfo(rpl).startT   = ripplePeriods(rpl, 1);
        rippleInfo(rpl).endT     = ripplePeriods(rpl, 2);
        rippleInfo(rpl).peakT    = ripplePeriods(rpl, 3);
        rippleInfo(rpl).duration = rippleInfo(rpl).endT - rippleInfo(rpl).startT;

        rippleInfo(rpl).peakRippleA = ripplePeriods(rpl, 4);
        
%         rippleInfo(rpl).peakChIdx   = ripplePeriods(rpl, 5);
        

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

    rippleInfo([rippleInfo.startT] > behavior.time(3,2)) = [];


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
    PYR.time = [];
    MUA.time = [];

    for unit = 1:numel(spikes)

        if spikes(unit).quality <= 3
            PYR.time = [PYR.time; spikes(unit).time];
        end

        MUA.time = [MUA.time; spikes(unit).time];
    end

    PYR.time = sort(PYR.time, 'ascend');
    MUA.time = sort(MUA.time, 'ascend');



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
        PBEInfo(ipbe).endT        = primaryPBEs.MUA(ipbe, 2);
        PBEInfo(ipbe).peakT       = primaryPBEs.MUA(ipbe, 3);
        PBEInfo(ipbe).duration    = PBEInfo(ipbe).endT - PBEInfo(ipbe).startT;
        
        PBEInfo(ipbe).peakMUA     = primaryPBEs.MUA(ipbe, 4);

        
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



speed_2add(speed_2add < 10) = 0;
speed_2add = (2*(speed_2add-min(speed_2add))/(max(speed_2add)-min(speed_2add))-1)*2^11;


ripple_2add(isnan(ripple_2add)) = 0;
ripple_2add(ripple_2add > 4) = 4;
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
velocity(isnan(velocity)) = interp1(timepnts(~isnan(velocity)), velocity(~isnan(velocity)), timepnts(isnan(velocity)));


% smoothing the speed
sigma = 0.25/(1/position_samplingRate); 
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
speed.v = conv(velocity, smoothwin, 'same'); 

% speed.v = velocity;
speed.t = timepnts;


end
