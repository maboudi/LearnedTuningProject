function basicProcedure_BapunNSD(sessionNumber)



addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/ReplayPreplayAnalyses'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Bapun_NSD_datasets';
cd(currDir)


sessionNames = {'RatN_Day2_2019-10-11_03-58-54'; 'RatS-Day2-2020-11-27_10-22-29'};


sessionName = sessionNames{sessionNumber};

sessionDir = fullfile(currDir, sessionName);

%% loading session data

varList = {'spike'; 'stability'; 'epochs'; 'position'; 'ripplesKM'};


for var = 1 : length(varList)
    load(fullfile(sessionDir, 'mat_files', [sessionName '.' varList{var} '.mat']))
end


% cov2cmfac = 1; 

mazeLimits = [-inf inf -inf 128; ... % RatN
              -inf inf -inf inf]; % Rat S

mazeShapes = {'L-shape'; 'L-shape'};

try
    behavior.time = [pre; maze; post];
catch
    behavior.time = [pre; maze1; post];
end

% calculate speed
position.TwoDLocation(:,1) = x;
position.TwoDLocation(:,2) = y;
position.TimeStamps        = time;


% filename = sprintf('%s.states_usingHMMforTheta.mat', sessionName);
% load(fullfile(sessionDir, filename), 'intervals')
% 
% WAKEstate  = intervals.WAKEstate;
% QWAKEstate = intervals.QWAKEstate;
% NREMstate  = intervals.NREMstate;
% REMstate   = intervals.REMstate;
% 
% 
% bvrTimeList = [WAKEstate; QWAKEstate; NREMstate; REMstate];
% bvrState = [4   *   ones(length(WAKEstate), 1); ...
%             3 *   ones(length(QWAKEstate), 1); ... 
%             1   *   ones(length(NREMstate), 1); ... 
%             2   *   ones(length(REMstate), 1)];
%         
% [~, ind] = sort(bvrTimeList(:, 1), 'ascend');
% bvrTimeList = bvrTimeList(ind, :);
% bvrTimeList(:,1) = bvrTimeList(:,1) - 1;
% 
% bvrState = bvrState(ind);


speed = calspeed(position); % animal's velocity


% % % maze theta 
% 
% filename = sprintf('%s.MAZE.theta.mat', sessionName);
% load(fullfile(sessionDir, filename), 'inTh', 't', 'thetaPeriods'); % t should be corrected 
% maze_theta_idx = inTh;
% maze_theta_t   = t + double(behavior.time(2,1));
% maze_theta_periods = thetaPeriods/1250;


% % theta
filename = sprintf('%s.thetaStates.mat', sessionName);
load(fullfile(sessionDir, filename), 'inTh', 'timePnts', 'thetaPeriods', 'theratio'); % t should be corrected


filename = sprintf('%s.SleepScoreMetrics.mat', sessionName);
load(fullfile(sessionDir, 'mat_files', filename), 'SleepScoreMetrics')


slowWave = SleepScoreMetrics.broadbandSlowWave;
swthresh    = SleepScoreMetrics.histsandthreshs.swthresh;
sw_timePnts = SleepScoreMetrics.t_clus;


%% sessioninfo


fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
mkdir(mainDir)


fileBase = fullfile(sessionDir, sessionName);

Par = LoadXml(fileBase); % Load the recording configurations from the xml file

fileinfo.nCh = Par.nChannels;


fileinfo.Fs  = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording
fileinfo.lfpSampleRate = 1;


% RUN high power theta periods
thetaPeriods = [];



% Position data

xpos = x';
ypos = y';
tpos = time';



fileinfo.xyt = [xpos ypos tpos]; 

% figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track


fileinfo.xyt(xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1) | ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3), 1:2) = NaN;


linearPos = linearizePosition2(fileinfo, behavior, mazeShapes{sessionNumber});


fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


[lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, 'bi', mainDir); 

turningPeriods(diff(turningPeriods, [], 2) < 1, :) = [];


if length(lapsStruct.RL) > length(lapsStruct.LR)
   lapsStruct.RL(1,:) = [];
end


totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
laps = zeros(totNumLaps, 2);
laps(1:2:totNumLaps, :)  = lapsStruct.LR;
laps(2:2:totNumLaps, :)  = lapsStruct.RL;


laps(:, 3) = 1:size(laps, 1); 


fileinfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end


runSpeedThresh = 10; % cm/s


%% combining nanPeriods with turning periods


linearpos = fileinfo.xyt2(:,1);

temp   = [0; diff(isnan(linearpos))]; 
starts = find(temp == 1);
ends   = find(temp == -1);

if ends(1) < starts(1)
    starts = [1; starts];
end

if starts(end) > ends(end)
    ends = [ends; length(xpos)];
end

nanPeriods = [tpos(starts) tpos(ends)];

nanPeriods(diff(nanPeriods, [], 2) < 1, :) = [];



nan_or_turningPeriods = [turningPeriods; nanPeriods];

[~, sortIdx] = sort(nan_or_turningPeriods(:, 1));

nan_or_turningPeriods = nan_or_turningPeriods(sortIdx, :);




primary = nan_or_turningPeriods;
secondary = [];

tmp_period = nan_or_turningPeriods(1, :);
    
for ip=2:size(primary,1)
    if (primary(ip,1)-tmp_period(2)) < 0
        tmp_period = [tmp_period(1), max(tmp_period(2), primary(ip,2))]; % merge two ripples
%     elseif abs(primary(ii,1)-tmp_period(1)) < min_isw_period
%         tmp_period = [min([tmp_period(1) primary(ii,1)]),max([tmp_period(2) primary(ii,2)])]; % merge two ripples
    else
        secondary = [secondary;tmp_period];
        tmp_period = primary(ip,:);
    end
end
    
nan_or_turningPeriods = secondary; 



% % plot the turning periods and nan periods separately

hold on
for ip = 1:size(turningPeriods, 1)
    
    line([turningPeriods(ip, :)], [160 160], 'linewidth', 3, 'color', 'm')
    
end

hold on

for ip = 1:size(nanPeriods, 1)
    
    line([nanPeriods(ip, :)], [163 163], 'linewidth', 3, 'color', 'c')
    
end

hold on

for ip = 1:size(nan_or_turningPeriods, 1)
    
    line([nan_or_turningPeriods(ip, :)], [166 166], 'linewidth', 3, 'color', 'g')
    
end



%% Behavioral states

% 
% bvrTimeList(:,1) = [states.startT]';
% bvrTimeList(:,2) = [states.endT]';
% 
% bvrState = [states.state]'; % 1:nrem, 2:rem, 3:quiet, 4:active



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'pyramidal';


for ii = 1:numel(spike)
    spike(ii).isStable =  prod(isStable(ii, :));
end

[spikeStruct,okUnits] = spikeBehaviorAnalysis_BG(spike, speed, laps, thetaPeriods, qual2consider, fileinfo);


shanks = [spike.shank]';
ids    = [spike.id]';
stableIdx = [spike.isStable]';

% temp   = [spike.id];
shanks = shanks(okUnits);
clus   = ids(okUnits);
stableIdx = stableIdx(okUnits);


%% Spatial tuning of the units

close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)


[spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], nan_or_turningPeriods, speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
[spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], nan_or_turningPeriods, speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
    
[spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], nan_or_turningPeriods, speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);

% firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileinfo, subfolder);

close all




%% Ripple Detection

% After visualizing deteceted events using different metrics (MUA and
% ripple amplitude) on Neuroscope, I realized that ripple is more reliable

% best_channels = BigRippleChannels(fileBase, fileinfo);
% 
% fileinfo.RippleChannels = best_channels;
% 
% 
% 
% threshSD = 2; % if we are summing the ripple envelope over channels for replay detection, a threshold of 2 s.d. works better
% [ripplePeriods, rippleLFP, ripplePower] = RippleDetect(fileBase, fileinfo, threshSD); % the results in seconds
% 
% 
% folderName = fullfile(mainDir, fileinfo.name,  'ripple');
% mkdir(folderName)
% 
% Filename = fullfile(folderName, [fileinfo.name '.prm.evt']);
% MakeEvtFile(ripplePeriods(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)


%%

% fileinfo.tbegin = double(behavior.time(1,1)); 
% fileinfo.tend = length(rippleLFP)/1250;
% 
% 
% sampleDur = 1/1250;
% samplingTimes = fileinfo.tbegin:sampleDur:fileinfo.tend;
% samplingTimes(1) = [];
% 
% 
% newSamplingTimes = fileinfo.tbegin:1e-3:fileinfo.tend;
% newSamplingTimes(1) = [];
% 
% 
% rippleLFP   = interp1(samplingTimes, rippleLFP, newSamplingTimes);
% rippleAmp   = interp1(samplingTimes, ripplePower, newSamplingTimes);
% 
% detectionParams = struct('lowband', 120, 'highband', 250, 'maxThresh', 2, 'minThresh', 1, 'minPeriod', 15, 'maxPeriod', 1000, 'isw', 0);
% 
% save([fileinfo.name '.ripplesKM.mat'], 'rippleLFP', 'rippleAmp', 'ripplePeriods', 'best_channels', 'detectionParams', '-v7.3')
% 


%% PBEs

%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the mean



fileinfo.lfpSampleRate = 1;

fileinfo.tbegin = double(behavior.time(1,1)); 
fileinfo.tend = double(behavior.time(3,2));



% exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4
exclude = [];

[primaryPBEs, sdat] = PBPeriods(spikeStruct, time_resolution, threshZ, exclude, fileinfo);


% bin the PBEs and removing flanking silence periods
[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, [1 2 3], fileinfo); 

% check overlap with ripples for each PBE
PBErippleIdx = ifContainRipples(secondaryPBEs, ripplePeriods);


% add the ripple overlap idx as a column to the PBE array 
secondaryPBEs(:, 5) = PBErippleIdx;




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



% secondaryPBEs = [secondaryPBEs zeros(nPBEs, 4)];
% 
% for ii= 1:nPBEs
%     
%     pbeCenter = secondaryPBEs(ii, 3);
%     
%     boutInd        = find(bvrTimeList(:,1) < pbeCenter & bvrTimeList(:,2) > pbeCenter, 1, 'first');
%     boutBrainState = bvrState(boutInd);
%     
%     secondaryPBEs(ii, 5 + boutBrainState) = 1;
%     
% end


baseStruct = struct('data', [], 'p', [], 'ts', [], 'pts', []);

PREbinnedPBEs       = baseStruct;
PREidx              = find(secondaryPBEs(:, 3) > behavior.time(1,1) & secondaryPBEs(:, 3) < behavior.time(1,2));
PREbinnedPBEs.data  = binnedPBEs(PREidx, :);
secondaryPBEs_PRE   = secondaryPBEs(PREidx, :);
% PREbinnedPBEs       = genSurrogates(PREbinnedPBEs);

RUNbinnedPBEs       = baseStruct;
RUNidx              = find(secondaryPBEs(:, 3) > behavior.time(2,1) & secondaryPBEs(:, 3) < behavior.time(2,2));
RUNbinnedPBEs.data  = binnedPBEs(RUNidx, :);
secondaryPBEs_RUN   = secondaryPBEs(RUNidx, :);
% RUNbinnedPBEs       = genSurrogates(RUNbinnedPBEs);

POSTbinnedPBEs       = baseStruct;
POSTidx              = find(secondaryPBEs(:, 3) > behavior.time(3,1) & secondaryPBEs(:, 3) < behavior.time(3,2));
POSTbinnedPBEs.data  = binnedPBEs(POSTidx, :);
secondaryPBEs_POST   = secondaryPBEs(POSTidx, :);
% POSTbinnedPBEs       = genSurrogates(POSTbinnedPBEs);


% [xpos_at_PBE, ypos_at_PBE, velocity_at_PBE] = calposvelatPBE(secondaryPBEs_RUN, fileinfo);
%%
mazetheta_at_PBE = inInterval(maze_theta_periods, secondaryPBEs_RUN);
%%
% mazetheta_at_PBE = interp1(maze_theta_t, maze_theta_idx, secondaryPBEs_RUN(:, 3));
% mazetheta_at_PBE = round(mazetheta_at_PBE);

% % storing the PBEs that ocurred only during the QW or NREM brain states

% secondaryPBEs_PRE  = secondaryPBEs_PRE(secondaryPBEs_PRE(:, 6) == 1 | secondaryPBEs_PRE(:, 8) == 1, :);
% secondaryPBEs_RUN  = secondaryPBEs_RUN(mazetheta_at_PBE == 0 , :);
% secondaryPBEs_POST = secondaryPBEs_POST(secondaryPBEs_POST(:, 6) == 1 | secondaryPBEs_POST(:, 8) == 1, :);
% 
% 
% secondaryPBEs2 = [secondaryPBEs_PRE; secondaryPBEs_RUN; secondaryPBEs_POST];
secondaryPBEs2 = secondaryPBEs(theratio_at_pbe < 1 | sw_at_pbe  > swthresh, :);

folderName = fullfile(mainDir,  'PopulationBurstEvents');
mkdir(folderName)

savePBEs(primaryPBEs, secondaryPBEs2, fileinfo, folderName)

%  
% PBErippleIdx = ifContainRipples(primaryPBEs, ripplePeriods);
% primaryPBEs = primaryPBEs(PBErippleIdx == 1, :);
% primaryPBEs(:, 5) = PBErippleIdx;
% 
% 
% primaryPBEs = [primaryPBEs zeros(size(primaryPBEs, 1), 4)];
% 
% for ii= 1:size(primaryPBEs, 1)
%     
%     pbeCenter = primaryPBEs(ii, 3);
%     
%     boutInd        = find(bvrTimeList(:,1) < pbeCenter & bvrTimeList(:,2) > pbeCenter, 1, 'first');
%     boutBrainState = bvrState(boutInd);
%     
%     primaryPBEs(ii, 5 + boutBrainState) = 1;
%     
% end
% 
% PREidx = find(primaryPBEs(:, 1) > double(behavior.time(1,1)) & primaryPBEs(:, 2) < double(behavior.time(1,2)));
% primaryPBEs_PRE   = primaryPBEs(PREidx, :);
% primaryPBEs_PRE = primaryPBEs_PRE(primaryPBEs_PRE(:, 6) == 1 | primaryPBEs_PRE(:, 8) == 1, :);
% 
% RUNidx = find(primaryPBEs(:, 1) > double(behavior.time(2,1)) & primaryPBEs(:, 2) < double(behavior.time(2,2)));
% primaryPBEs_RUN   = primaryPBEs(RUNidx, :);
% 
% POSTidx = find(primaryPBEs(:, 1) > double(behavior.time(3,1)) & primaryPBEs(:, 2) < double(behavior.time(3,2)));
% primaryPBEs_POST   = primaryPBEs(POSTidx, :);
% primaryPBEs_POST = primaryPBEs_POST(primaryPBEs_POST(:, 6) == 1 | primaryPBEs_POST(:, 8) == 1, :);
% 
% 
% qclus = [1 2 3];
% 
% [PREbinnedPBEs, secondaryPBEs_PRE]   = finalBinningResult(primaryPBEs_PRE, spikeStruct, qclus, fileinfo); 
% 
% [RUNbinnedPBEs, secondaryPBEs_RUN]   = finalBinningResult(primaryPBEs_RUN, spikeStruct, qclus, fileinfo);
% 
% [POSTbinnedPBEs, secondaryPBEs_POST] = finalBinningResult(primaryPBEs_POST, spikeStruct, qclus, fileinfo); 
% 
% 
% 
% primaryPBEs = [primaryPBEs_PRE; primaryPBEs_RUN; primaryPBEs_POST];
% secondaryPBEs = [secondaryPBEs_PRE; secondaryPBEs_RUN; secondaryPBEs_POST];
% 
% folderName = fullfile(mainDir, fileinfo.name,  'PopulationBurstEvents');
% mkdir(folderName)
% 
% savePBEs(primaryPBEs, secondaryPBEs, fileinfo, folderName)
% 
% 
% 
% %%  Bayesian decoding
% 
% BayesianReplayDetection_GrosmarkReclu_nov2020


end



%% functions


function speed = calspeed(position)


cov2cmfac = 1; % the position are in m in the dataset, hence positions must be multiplied with 100


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

sigma = 15; %%% smoothing the speed, the positon sampling rate is around 40 Hz (0.0256 sec sampling period), so duration of the sigma is about 25.6 ms times sigma (25.6*20 ~ 512 ms)
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
velocity = conv(velocity, smoothwin, 'same'); 


speed.v = velocity;
speed.t = timepnts;


end

function [xpos_at_PBE, ypos_at_PBE, velocity_at_PBE] = calposvelatPBE(PBEs, fileinfo)


tbegin = fileinfo.tbegin;
tend   = fileinfo.tend;

xyt = fileinfo.xyt;



withinRange = find(xyt(:,3) >= tbegin & xyt(:,3) <= tend); % selecting only the position samples occured during neural recording
xyt = xyt(withinRange, :);

timepnts = (xyt(:,3) - tbegin)/fileinfo.Fs; % u seconds to seconds


xpos = xyt(:,1); % Converting the position data (from pixel number to centimeter)
ypos = xyt(:,2);


linearpos = fileinfo.xyt2(:, 2);
linearpos = linearpos(withinRange);


%%%% SPEED %%%%%%%%%%%%

diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difftt = [1; diff(timepnts)];

speed = sqrt(diffx.^2 + diffy.^2)./difftt; 

% smoothing the speed

sigma = 15; %%% smoothing the speed, the positon sampling rate is around 30 Hz, so duration of the sigma is about 500 ms
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
speed = conv(speed, smoothwin, 'same');


% % track position and speed at the PBEs 
nPBEs = size(PBEs, 1);

xpos_at_PBE = zeros(nPBEs, 1);
ypos_at_PBE = zeros(nPBEs, 1);

velocity_at_PBE = zeros(nPBEs, 1);

for ipbe = 1:nPBEs
    
    currPBE = PBEs(ipbe, 3);
    
    velocity_at_PBE(ipbe) = interp1(timepnts, speed, currPBE);
    xpos_at_PBE(ipbe) = interp1(timepnts, xpos, currPBE);
    ypos_at_PBE(ipbe) = interp1(timepnts, ypos, currPBE);
    
    if mod(ipbe, 50) == 0
        fprintf('.')
    end
    
end


end

function mazetheta_at_PBE = inInterval(intervals, pbes)

npbes = size(pbes, 1);
mazetheta_at_PBE = zeros(npbes, 1);

for ipbe = 1:npbes
    
    intIdx = find((pbes(ipbe, 1) < intervals(:, 1) & pbes(ipbe, 2) > intervals(:, 1)) | ...
                  (pbes(ipbe, 1) < intervals(:, 2) & pbes(ipbe, 2) > intervals(:, 2)) | ...
                  (pbes(ipbe, 1) > intervals(:, 1) & pbes(ipbe, 2) < intervals(:, 2)), 1, 'first');
    
    if ~isempty(intIdx)
        mazetheta_at_PBE(ipbe) = 1;
    end
    
end

end
