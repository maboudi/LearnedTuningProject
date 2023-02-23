function PBEdetection_MUAvsonlyPyr(sessionNumber)



parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

addpath(genpath(fullfile(parentDir, '/ReplayPreplayAnalyses')))

datasetDir = fullfile(parentDir, '/Grosmark_reclustered');
dataPath_OriginalClusters = fullfile(parentDir, 'Grosmark_originalClusters');
stateDetectionResultsDir = fullfile(parentDir, 'StateDetectionResults');

outputDir = fullfile(parentDir, 'MUAvsonlyPyr');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end




%% LOADING SESSION DATA

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
emg         = SleepScoreMetrics.EMG;
load(fullfile(stateDetectionResultsDir, sessionName, sprintf([sessionName '.emgThresh.mat'])), 'emgThresh');


% bad channels to exclude for lpf event detections (e.g., ripples)

filename = sprintf([sessionName '.sessionInfo.mat']);
load(fullfile(stateDetectionResultsDir, sessionName, filename), 'sessionInfo')

badChannels = sessionInfo.badchannels;




%% SESSION'S INFORMATION

baseFile = fullfile(dataPath_OriginalClusters, sessionName, sessionName);
Par = LoadXml(baseFile); % Load the recording configurations from the xml file

fileInfo.sessionName   = sessionName;
fileInfo.animal        = sessionName(1:strfind(sessionName, '_')-1);
fileInfo.xyt           = [];
fileInfo.xyt2          = []; % linearized position and lap information
fileInfo.tbegin        = behavior.time(1,1);
fileInfo.tend          = behavior.time(3,2);
fileInfo.Fs            = Par.SampleRate; 
fileInfo.lfpSampleRate = Par.lfpSampleRate;
fileInfo.timeUnit      = 1; % in sec
fileInfo.nCh           = Par.nChannels;
fileInfo.badChannels   = badChannels;
fileInfo.behavior      = behavior;



spikes = spikes([spikes.quality] <= 3); % limiting to only pyramidal units

fileInfo.stablePyr = find([spikes.StablePrePostFiner] == 1); % pyramidal units that are stable
fileInfo.okUnits   = 1:numel(spikes);

% Position data

xpos = position.TwoDLocation(:,1)* 100; % convert to centimeters
ypos = position.TwoDLocation(:,2)* 100;
tpos = position.TimeStamps';

fileInfo.xyt = [xpos ypos tpos]; 


speed = calspeed(fileInfo); % animal's velocity





%% LINEARIZE POSTION AND CALCULATE LAPS

fileInfo.xyt(xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1) | ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3), 1:2) = NaN;


linearPos = linearizePosition(fileInfo, behavior, mazeShapes{sessionNumber});


fileInfo.xyt2(:, 1) = linearPos; 
fileInfo.xyt2(:, 2) = fileInfo.xyt(:, 3);


if ismember(sessionNumber, [2 5 8]) 
    direction = 'uni'; 
    laps = calculateLapTimings(fileInfo, speed, direction, storageDir); 
else
    direction = 'bi';
    lapsStruct = calculateLapTimings(fileInfo, speed, direction, storageDir); 
    
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



%% CALCULATING POSITION, VELOCITY, ETC AT EACH SPIKE (THESE WILL BE NEEDED FOR CALCULATION OF PLACE FIELDS)

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



%% MULTIUINT SPIKES FOR DETECTION OF POPULATION BURSTS


nShanks = length(dir([baseFile '.res.*']));


MUA.time = [];
pyrInt.time = [];
for shank = 1:nShanks
    
    resContent  = load([baseFile '.res.' num2str(shank)]);
    
    cluContent  = load([baseFile '.clu.' num2str(shank)]);
    clusterIDs  = cluContent(2:end); % the first entry is the total number of clusters
    
    
    % all MUA spikes: pyr + int + unclustered
    spikes2Add = resContent(clusterIDs ~=0); % to exclude the noise cluster
    MUA.time = [MUA.time; spikes2Add]; % the first is the total number of clusters detected for the shank   
    
    % pyr + int
    spikes2Add = resContent(~ismember(clusterIDs, [0 1])); % to exclude the noise cluster
    pyrInt.time = [pyrInt.time; spikes2Add]; % the first is the total number of clusters detected for the shank
    
end

MUA.time = sort(MUA.time, 'ascend');
MUA.time = MUA.time./fileInfo.Fs;


pyrInt.time = sort(pyrInt.time, 'ascend');
pyrInt.time = pyrInt.time ./fileInfo.Fs;


%% SPATIAL TUNNIG OF THE UNITS


close all

subfolder = fullfile(storageDir, 'placeFields');
mkdir(subfolder)


%%% 1D spatial tuning: using linearized position

posBinSize = 2;

if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular

    spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);

else 
    spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'LR', posBinSize, fileInfo, subfolder);
    spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'RL', posBinSize, fileInfo, subfolder);
    
    spikes = spatialTuning_1D_June21(spikes, behavior, [], [], speed, runSpeedThresh, 'uni', posBinSize, fileInfo, subfolder);
end

close all



%% RIPPLE DETECTION
 

best_channels = BigRippleChannels(baseFile, fileInfo);
fileInfo.RippleChannels = best_channels; % check if there are no bad channels among the best_channels 


threshSD = 2;
[ripplePeriods, rippleLFP, ripplePower] = RippleDetect(baseFile, fileInfo, threshSD); % the results in seconds


% downsample the ripple band lfp or power to 1kHz (1ms sampling period) % maybe we don't
% need to downsample at all?

ts    = 1/Par.lfpSampleRate;
tPnts = ts:ts:fileInfo.tend;

tq    = 1e-3:1e-3:fileInfo.tend;


rippleLFP_adjust = zeros(numel(tq), size(rippleLFP, 2));
for ich = 1:size(rippleLFP, 2)
    rippleLFP_adjust(:, ich) = interp1(tPnts, rippleLFP(:, ich), tq);
end

ripplePower_adjust = interp1(tPnts, ripplePower, tq);
  



%% GENERATE RIPPLE INFO

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
    
%     rippleInfo(rpl).peakMUA  = max(sdat(floor(rippleInfo(rpl).startT*1e3):floor(rippleInfo(rpl).endT*1e3)));
    
    
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
              


%% DETERMINING POPULATION BURST EVENTS (PBEs)

fileInfo.tbegin = behavior.time(1,1); 
fileInfo.tend   = behavior.time(3,2);


detectParams.time_resolution    = 0.001;
detectParams.exclude            = [];
detectParams.smoothingSigma     = 0.02;
detectParams.smoothinghalfwidth = 0.06;
detectParams.minDur             = 3*0.02; % number of bins
detectParams.maxDur             = 30*0.02;
detectParams.minFiringUnits     = max(5, 0.1*numel(fileInfo.stablePyr));
detectParams.maxVelocity        = 10; 
detectParams.thRatioThresh      = 1;
detectParams.rippleAThresh      = 3;



thresholds = [1 2 3];
detectionUnits = {'MUA'; 'pyrInt';'onlyPyr'};
rippleAThreshs = [2 5]; 


rippleOverlap_ratio_of_ripples = zeros(numel(thresholds), numel(detectionUnits), numel(rippleAThreshs));
rippleOverlap_ratio_of_PBEs    = zeros(numel(thresholds), numel(detectionUnits), numel(rippleAThreshs));


for ithresh = 1:numel(thresholds)
    ithresh
    for idu = 1:numel(detectionUnits)


        detectParams.threshold = thresholds(ithresh);
        detectParams.MUAorPyr  = detectionUnits{idu};

        
        if strcmp(detectParams.MUAorPyr, 'MUA') 
            [primaryPBEs, sdat] = PBPeriods(MUA, detectParams, fileInfo);
        elseif strcmp(detectParams.MUAorPyr, 'pyrInt')
            [primaryPBEs, sdat] = PBPeriods(pyrInt, detectParams, fileInfo);
        elseif strcmp(detectParams.MUAorPyr, 'onlyPyr')
            [primaryPBEs, sdat] = PBPeriods(spikes, detectParams, fileInfo);
        end


        nPBEs = size(primaryPBEs, 1);
        PBEInfo = struct('sessionName', [], 'PBEno', [], ...
                            'startT', [], 'endT', [], 'peakT', [], ...
                            'startT_adj', [], 'endT_adj', [], 'peakT_adj', [], ...
                            'peakMUA', [], 'MUAwave', [], ...
                            'nFiringUnits', [], 'duration', [], ...
                            'peakRippleA', [], 'rippleWave', [], ...
                            'thetaRatio', [], 'SWA', [], 'emg', [], ...
                            'epoch', [], 'brainState', [], ...
                            'fr_1msbin', [], 'fr_20msbin', []);            


                        
        binDur = 0.02; % in sec

        PBEs_all = [];

        for ipbe = 1:nPBEs

            PBE_primary = primaryPBEs(ipbe, :);

            [binnedPBE, PBE_adj, nFiringUnits, PBElength] = finalBinningResult(PBE_primary, spikes, binDur, fileInfo); 

        %     if PBElength > 0
        %         PBEs_all = [PBEs_all; PBE_adj];
        %     else
        %         PBEs_all = [PBEs_all; nan(size(PBE_adj))];
        %     end


        % %
            PBEs_all = [PBEs_all; PBE_primary];
        % %

            PBEInfo(ipbe).sessionName = sessionName;
            PBEInfo(ipbe).PBEno       = ipbe;

            PBEInfo(ipbe).startT   = PBE_primary(1);
            PBEInfo(ipbe).endT     = PBE_primary(2);
            PBEInfo(ipbe).peakT    = PBE_primary(3);

            PBEInfo(ipbe).startT_adj = PBE_adj(1);
            PBEInfo(ipbe).endT_adj   = PBE_adj(2);
            PBEInfo(ipbe).peakT_adj  = PBE_adj(3);



            PBEInfo(ipbe).peakMUA = PBE_primary(4);

            try
                PBEInfo(ipbe).MUAwave = sdat(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500));
            catch
                PBEInfo(ipbe).MUAwave = [];
            end

            try
                PBEInfo(ipbe).rippleWave = rippleLFP_adjust(floor((PBEInfo(ipbe).startT*1e3)-500) : (floor(PBEInfo(ipbe).endT*1e3)+500), 1);
            catch
                PBEInfo(ipbe).rippleWave = [];
            end


        %     PBEInfo(ipbe).duration     = PBElength;

        % %
            PBEInfo(ipbe).duration     = PBE_primary(2) - PBE_primary(1);
        % %   

            PBEInfo(ipbe).nFiringUnits = nFiringUnits;


        %     PBEInfo(ipbe).peakRippleA = interp1(tq, ripplePower_adjust, PBEInfo(ipbe).peakT); 
            PBEInfo(ipbe).peakRippleA = max(ripplePower_adjust(floor(PBEInfo(ipbe).startT*1e3): floor(PBEInfo(ipbe).endT*1e3)));


            PBEInfo(ipbe).thetaRatio  = interp1(timePnts, theratio, PBEInfo(ipbe).peakT);
            PBEInfo(ipbe).SWA         = interp1(sw_timePnts, slowWave, PBEInfo(ipbe).peakT); % what exactly is the SWA?
            PBEInfo(ipbe).emg         = interp1(sw_timePnts, emg, PBEInfo(ipbe).peakT);

            PBEInfo(ipbe).velocity    = interp1(speed.t, speed.v, PBEInfo(ipbe).peakT);
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
            if PBEInfo(ipbe).SWA > swthresh
                PBEInfo(ipbe).brainState = 'NREM';
            elseif PBEInfo(ipbe).thetaRatio < detectParams.thRatioThresh || PBEInfo(ipbe).peakRippleA > detectParams.rippleAThresh
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


        
        % store the PBEs to visualize using Neuroscope
        
        PBEs_accepted = PBEs_all(ismember({PBEInfo.brainState}, {'NREM'; 'QW'}) & ...
                                [PBEInfo.duration] >= detectParams.minDur & [PBEInfo.duration] <= detectParams.maxDur & ...
                                [PBEInfo.velocity] < detectParams.maxVelocity, :);
        %                         [PBEInfo.nFiringUnits] >= detectParams.minFiringUnits & ...


        % find the overlap between ripple events and detected PBEs using any of the
        % methods

        nPBE = size(PBEs_accepted, 1);

        for irt = 1:numel(rippleAThreshs)

            currentRippleThresh = rippleAThreshs(irt);

            highAmpRipples = rippleInfo([rippleInfo.peakRippleA] > currentRippleThresh);

            nRipples = numel(highAmpRipples);


            rippleIncluded = zeros(nRipples, 1);
            for ir = 1:nRipples 

                temp = find(PBEs_accepted(:,1) < highAmpRipples(ir).peakT & PBEs_accepted(:,2) > highAmpRipples(ir).peakT, 1, 'first');

                if ~isempty(temp)
                    rippleIncluded(ir) = 1;
                end

            end

            rippleOverlap_ratio_of_ripples(ithresh, idu, irt) = length(find(rippleIncluded))/nRipples;



            highAmpRipples_peakT = [highAmpRipples.peakT];

            PBEoverlapped = zeros(nPBE, 1);
            for ip = 1:nPBE 

                temp = find(highAmpRipples_peakT > PBEs_accepted(ip, 1) & highAmpRipples_peakT < PBEs_accepted(ip, 2), 1, 'first');

                if ~isempty(temp)
                    PBEoverlapped(ip) = 1;
                end

            end

            rippleOverlap_ratio_of_PBEs(ithresh, idu, irt)    = length(find(PBEoverlapped))/size(PBEs_accepted, 1);

        end


    end
end



folderName = fullfile(storageDir,  'PopulationBurstEvents');
if ~exist(folderName, 'dir'); mkdir(folderName); end


fileName = [sessionName '_rippleOverlap_fraction_of_ripples.mat'];
save(fullfile(folderName, fileName), 'rippleOverlap_ratio_of_ripples', 'thresholds', 'detectionUnits', 'rippleAThreshs')

fileName = [sessionName '_rippleOverlap_fraction_of_PBEs.mat'];
save(fullfile(folderName, fileName), 'rippleOverlap_ratio_of_PBEs', 'thresholds', 'detectionUnits', 'rippleAThreshs')



%%
figure; 


% ripple overlap as a ratio of total number of ripples
% all ripples: ripple power > 2
subplot(2,2, 1)

plot(rippleOverlap_ratio_of_ripples(:, 1, 1), 'r', 'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 2, 1), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 3, 1), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('all ripples(>2z)')
xlim([0 4])
xticks([1 2 3])

legend('MUA', 'pyr+int', 'only pyr', 'location', 'best')


% high amplitude ripples : ripple power > 5
subplot(2,2, 3)

plot(rippleOverlap_ratio_of_ripples(:, 1, 2), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 2, 2), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_ripples(:, 3, 2), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of ripples'; 'overlapping w PBEs'})
title('stronger ripples(>5z)')
xlim([0 4])
xticks([1 2 3])


% ripple overlap as a ratio of total number of PBEs
% ripple power  > 2
subplot(2,2, 2)

plot(rippleOverlap_ratio_of_PBEs(:, 1, 1), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 2, 1), 'g', 'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 3, 1), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)


xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('all ripples(>2z)')
xlim([0 4])
xticks([1 2 3])


%ripple power > 5
subplot(2,2, 4)

plot(rippleOverlap_ratio_of_PBEs(:, 1, 2), 'r',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 2, 2), 'g',  'marker', 's', 'linewidth', 2, 'markersize', 5)
hold on
plot(rippleOverlap_ratio_of_PBEs(:, 3, 2), 'b',  'marker', 's', 'linewidth', 2, 'markersize', 5)

xlabel('PBE thresholds(z)')
ylabel({'fraction of PBEs'; 'overlapping w ripples'})
title('stronger ripples(>5z)')
xlim([0 4])
xticks([1 2 3])


suptitle(sessionName)

folderName = fullfile(storageDir,  'PopulationBurstEvents');
if ~exist(folderName, 'dir'); mkdir(folderName); end


fileName = [sessionName '_rippleOverlap_MUAvsPyr'];
savepdf(gcf, fullfile(folderName, fileName), '-dpdf')




end


function speed = calspeed(fileInfo)


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
sigma = 5; %%% smoothing the speed, the positon sampling rate is around 40 Hz (0.0256 sec sampling period), so duration of the sigma is about 25.6 ms times sigma (25.6*15 ~ 384 ms)
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
speed.v = conv(velocity, smoothwin, 'same'); 

% speed.v = velocity;
speed.t = timepnts;

end
