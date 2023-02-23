function basicProcedure_GrosmarkReclu_May21(sessionNumber)


% addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/ReplayPreplayAnalyses'))
% 
% currDir = '/nfs/turbo/umms-kdiba/Kourosh/Grosmark_reclustered';
% cd(currDir)
% 
% originalDataDir = '/nfs/turbo/umms-kdiba/Kourosh/Grosmark_originalClusters';


addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/ReplayPreplayAnalyses'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_BapunReclustered';
cd(currDir)

originalDataDir = '/data/GrossmarkOriginalData';



%% loading session data

VarList = {'spikes','behavior','position'};


% note that the spikes within the .mat file are from just the pyramidal
% units, for MUA I need to load the clu and res files again


% Although no difference regarding the behavior and position data

for var = 1 : length(VarList)
    load(fullfile(currDir, ['/CRCNSReclustered-' VarList{var} '.mat']))
end



sessionNames = fieldnames(spikes);
nSessions = numel(sessionNames);


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100


mazeLimits = [-inf inf -inf 13; ...% 0 120 -inf 13 ;... % Achilles_10252013
              -inf inf -inf inf;... % Achilles_11012013
              -inf inf -140 inf;... % Buddy_06272013 % bad place fields
%               -inf inf -inf inf;... % Cicero_09012014
%               -15 inf -inf 115;... % Cicero_09102014
              -inf inf -inf inf;... % Cicero_09172014
%               -inf inf -inf inf;... % Gatsby_08022013
              -20  inf -inf 110];   % Gatsby_08282013
          
mazeShapes = {'linear', 'circular', 'linear', 'linear', 'circular'};




sessionName = sessionNames{sessionNumber};


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['Position.' sessionName]);

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];


%% Behavioral states


filename = sprintf('%s.thetaStates.mat', sessionName);
load(fullfile(originalDataDir, sessionName, filename), 'timePnts', 'theratio'); % t should be corrected


filename = sprintf('%s.SleepScoreMetrics.mat', sessionName);
load(fullfile(originalDataDir, sessionName, filename), 'SleepScoreMetrics')


slowWave    = SleepScoreMetrics.broadbandSlowWave;
swthresh    = SleepScoreMetrics.histsandthreshs.swthresh;
sw_timePnts = SleepScoreMetrics.t_clus;



%% sessioninfo

unsortedData = fullfile(originalDataDir, sessionName, sessionName);

fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
mkdir(mainDir)



Par = LoadXml(unsortedData); % Load the recording configurations from the xml file

fileinfo.nCh = Par.nChannels;


fileinfo.Fs  = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording
fileinfo.lfpSampleRate = 1;




% RUN high power theta periods

thetaPeriods = importEVTfile(unsortedData); % make sure that we know details about procedures and parameters which were used for the detection of these periods 
% thetaPeriods = [];



% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

tpos = position.TimeStamps';


fileinfo.xyt = [xpos ypos tpos]; 


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


speed = calspeed(position); % animal's velocity

%% Based on the lookup table exclude the positions outside the track


fileinfo.xyt(xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1) | ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3), 1:2) = NaN;


% xpos = fileinfo.xyt(:, 1);
% ypos = fileinfo.xyt(:, 2);

% figure; plot(xpos, ypos, '.', 'markersize', 3)


linearPos = linearizePosition(fileinfo, behavior, mazeShapes{sessionNumber});

fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


% laps = calculateLapTimings(fileinfo, speed, mainDir);

if ismember(sessionNumber, [2 5]) 
    direction = 'uni'; 
    [laps, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
else
    direction = 'bi';
    [lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
    
    if length(lapsStruct.RL) > length(lapsStruct.LR)
       lapsStruct.RL(1,:) = [];
%        behavior.MazeEpoch(1) = lapsStruct.LR(1,1);
    end
    
    totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
    laps = zeros(totNumLaps, 2);
    laps(1:2:totNumLaps, :)  = lapsStruct.LR;
    laps(2:2:totNumLaps, :)  = lapsStruct.RL;
end

laps(:, 3) = 1:size(laps, 1); 


fileinfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end

runSpeedThresh = 10; % cm/s



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

[spikeStruct,okUnits] = spikeBehaviorAnalysis2(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo);

temp = [spikes.id];
shanks = temp(2*okUnits - 1);
clus   = temp(2*okUnits);



%% Spatial tuning of the units

close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)

%%% 1D spatial tuning: using linearized position


% the place fields are calculated separately for the left and right
% direction


if ismember(sessionNumber, [2 5]) % sessions 2, 5, and 8 are circular

    [spatialTunings, PF_sorted, runTemplate, spatialInfo, conslapsRatio, diffWithAvg] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings, [], spatialInfo, [], conslapsRatio, [], behavior, [], fileinfo, subfolder);
else

    [spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
    [spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
    
    [spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);

%     firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileinfo, subfolder);
end


close all



%% Ripple Detection

% After visualizing deteceted events using different metrics (MUA and
% ripple ampplitude) on Neuroscope, I realized that ripple is more reliable
 
best_channels = BigRippleChannels(unsortedData, fileinfo);
fileinfo.RippleChannels = best_channels(1:12);


threshSD = 2;
[ripplePeriods, rippleLFP, ripplePower] = RippleDetect(unsortedData, fileinfo, threshSD); % the results in seconds


ts = 1/Par.lfpSampleRate;
tPnts = ts:ts:fileinfo.tend;

tq = 1e-3:1e-3:fileinfo.tend;

ripplePower_adjust = interp1(tPnts, ripplePower, tq);




%% PBEs

%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the mean



fileinfo.lfpSampleRate = 1;

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(3,2);




exclude =[];
% exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4


[primaryPBEs, sdat] = PBPeriods(spikeStruct, time_resolution, threshZ, exclude, fileinfo);



qclus = [1 2 3];
[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileinfo); 

PBErippleIdx = ifContainRipples(secondaryPBEs, ripplePeriods);
secondaryPBEs(:, 5) = PBErippleIdx;


nPBEs = size(binnedPBEs, 1);

PBE_ripplePower = zeros(nPBEs, 1);
for ii = 1:nPBEs
    currPBE = floor(secondaryPBEs(ii, 1:2)*1e3);
    PBE_ripplePower(ii) = max(ripplePower_adjust(currPBE(1): currPBE(2)));
end



% the brain state at which each PBE happened
nPBEs = size(binnedPBEs, 1);

secondaryPBEs = [secondaryPBEs zeros(nPBEs, 2)];


theratio_at_pbe = zeros(nPBEs, 1);
sw_at_pbe       = zeros(nPBEs, 1);

for ii = 1:nPBEs
    
   pbe = secondaryPBEs(ii, 1:3);
   
   theratio_at_pbe(ii) = interp1(timePnts, theratio, pbe(3));
   sw_at_pbe(ii)       = interp1(sw_timePnts, slowWave, pbe(3)); 
    
end

figure; scatter(theratio_at_pbe, PBE_ripplePower, 5, 'k', 'filled'); xlabel('theta ratio at PBE'); ylabel('ripple power at PBE');

secondaryPBEs2 = secondaryPBEs(theratio_at_pbe < 1 | sw_at_pbe  > swthresh | PBE_ripplePower > 2.5, :);
% secondaryPBEs2 = secondaryPBEs( PBE_ripplePower < 2.5 & theratio_at_pbe < 0.5, :);


folderName = fullfile(mainDir,  'PopulationBurstEvents');
mkdir(folderName)

savePBEs(primaryPBEs, secondaryPBEs2, fileinfo, folderName)



%%  Bayesian decoding

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