clear;
clc;


parentDir = '/home/kouroshmaboudi/Documents/NCMLproject';

rr = dir(fullfile(parentDir, 'assemblyTuning_finalResults'));


for sessionNumber = 11

sessionName = rr(sessionNumber+2).name

basePath = fullfile(parentDir, 'assemblyTuning_finalResults', sessionName);



% place fields
s = load(fullfile(basePath, 'spikes', [sessionName '.spikes.mat']), 'fileInfo', 'spikes_pyr');
fileInfo = s.fileInfo;
spikes   = s.spikes_pyr;
position = fileInfo.linearPos;

%
epochs   = fileInfo.behavior.time;
maze     = epochs(2,:);



nPosBins = numel(spikes(1).spatialTuning_smoothed.uni);
nUnits   = numel(spikes); 


% non-directional spatial tuning

spatialTunings = nan(nUnits, nPosBins);
peakPFfiring   = nan(nUnits, 1);
peakLoc        = nan(nUnits, 1);

for iUnit = 1:nUnits
    spatialTunings(iUnit, :) = spikes(iUnit).spatialTuning_smoothed.uni; 
    [~, peakLoc(iUnit)] = max(spatialTunings(iUnit, :));
    
    peakPFfiring(iUnit) = spikes(iUnit).peakFR.uni;

end
[~, PFSortIdx] = sort(peakLoc, 'ascend');



% cluster quality
% load(fullfile(basePath, 'spikes', [sessionName '.clusterQuality.mat']), 'clusterQuality')


% re-maze time periods
load(fullfile(parentDir, 'Bapun_NSD_datasets', sessionName, [sessionName '.epochs.mat']))

if exist('re_maze', 'var')
    remaze = double(re_maze);
elseif exist('maze2', 'var')
    remaze = double(maze2);
else
    error('There is no remaze in this session!')
end


% if exist('maze', 'var')
%     remaze = double(maze);
% elseif exist('maze1', 'var')
%     remaze = double(maze1);
% else
%     error('There is no remaze in this session!')
% end


% theta periods (we need specifically those ocurred in the remaze period)


%(1)
% load(fullfile(parentDir, 'StateDetectionResults', sessionName, [sessionName '.thetaStates_including_remaze.mat']), 'thetaPeriods')


% % (2)
% % load(fullfile(basePath, 'spikes', [sessionName '.brainState.mat']), 'brainState')
% % allPeriods = brainState.periods;
% % thetaPeriods = allPeriods(allPeriods(:, 3)== 3, 1:2);
% 
% 
% 
% (3)

% loading theta info
% filename = sprintf([sessionName '.thetaStates_including_remaze_stricterDetection.mat']);
% load(fullfile(parentDir, 'StateDetectionResults', sessionName, filename), 'timePnts', 'theratio', 'thetaPeriods_all'); 
% brainState.thetaTimePnts = timePnts;
% brainState.theratio      = theratio;


%%%%
fileName = fullfile(parentDir, 'StateDetectionResults', sessionName, [sessionName '.brainStateDetection_HMMtheta_EMG_SWS_SchmidtTrigger.mat']);

if isfile(fileName)

    load(fileName, 'brainState'); 
    bouts = brainState.bouts;
    theta_periods = bouts(bouts(:,3) == 3, 1:2);

else

    brainState = fileInfo.brainStates;
    theta_periods = brainState(brainState(:,3) == 4, 1:2); % in Hiro's datatsets, Active wake is 4

end
%%%%%


% 
% 
% % loading slow wave power/slope 
% filename = sprintf([sessionName '.SleepScoreMetrics.mat']);
% load(fullfile(parentDir, 'StateDetectionResults', sessionName, filename), 'SleepScoreMetrics')
% 
% 
% brainState.slowWave        = SleepScoreMetrics.broadbandSlowWave;
% brainState.sw_emg_timePnts = SleepScoreMetrics.t_clus;
% 
% try 
%   load(fullfile(parentDir, 'StateDetectionResults', sessionName, sprintf([sessionName '.swthresh.mat'])), 'swthresh'); 
%     brainState.swthresh = swthresh;
% catch
%     brainState.swthresh = SleepScoreMetrics.histsandthreshs.swthresh;
% end
% 
% 
% 
% % loading emg 
% 
% brainState.emg = SleepScoreMetrics.EMG;
% try
%     load(fullfile(parentDir, 'StateDetectionResults', sessionName, sprintf([sessionName '.emgThresh.mat'])), 'emgThresh');
%     brainState.emgThresh = emgThresh;
% catch
%     brainState.emgThresh = SleepScoreMetrics.histsandthreshs.EMGthresh; 
% end
% 
% 
% 
% % loading MUA
% 
% [sdat, sdat_t] = calculateMUA(spikes);
% 
% 
% % Bouts with different brain states
% 
% % thetaRatio2 = interp1(brainState.thetaTimePnts, brainState.theratio, brainState.sw_emg_timePnts);
% 
% stateIdx = zeros(size(brainState.sw_emg_timePnts));
% 
% 
% timePoints = sdat_t;
% theRatio   = interp1(brainState.thetaTimePnts, brainState.theratio, timePoints);
% slowWave   = interp1(brainState.sw_emg_timePnts, brainState.slowWave, timePoints);
% emg        = interp1(brainState.sw_emg_timePnts, brainState.emg, timePoints);
% linearPos  = interp1(position(:, 2), position(:, 1), timePoints);
% % sdat     = interp1(sdat_t, sdat, timePoints);
% 
% swThresh  = brainState.swthresh;
% 
% 
% % test something ...
% 
% % theRatio_run = nan(size(theRatio));
% % 
% % idx = timePoints > remaze(1) & timePoints < remaze(2) & theRatio > 1 & emg < brainState.emgThresh+0.2; 
% % theRatio_run(idx) = linearPos(idx);
% % 
% % idx = timePoints > remaze(1) & timePoints < remaze(2);
% % linearPos2 = nan(size(linearPos));
% % linearPos2(idx) = linearPos(idx);
% % linearPos = linearPos2;
% % 
% % figure; scatter(linearPos, timePoints, 2, 'filled')
% % 
% % hold on
% % scatter(theRatio_run, timePoints, 2, 'r', 'filled')
% 
% 
% 
% %
% 
% 
% 
% if sessionNumber == 8
%     theThresh = 1.5;
%     emgThresh = brainState.emgThresh+0.1;
% elseif sessionNumber == 9
%     theThresh = 1;
%     emgThresh = brainState.emgThresh;
% end
% 
% 
% stateIdx(slowWave' >= swThresh) = 1; %NREM
% stateIdx(slowWave' < swThresh & theRatio' < theThresh) = 4; % QW
% stateIdx(slowWave' < swThresh & theRatio' >= theThresh & emg' < emgThresh) = 2; %REM
% stateIdx(slowWave' < swThresh & theRatio' >= theThresh & emg' > emgThresh) = 3; % WAKE
% 
% 
% brainState.periods = cell(4,1);
% 
% for iState = 1:4
%    
%     currInd = stateIdx == iState;
%     
%     crossup   = find([0; diff(currInd)] == 1);
%     crossdown = find([0; diff(currInd)] == -1);
%     
%     if crossdown(1) < crossup(1)
%         crossdown(1) = [];
%     end
%     
%     if crossup(end) > crossdown(end)
%         crossup(end) = [];
%     end
%         
%     brainState.periods{iState} = [timePoints([crossup crossdown])  iState*ones(numel(crossup), 1)];    
% end
% 
% 
% brainState.periods = cell2mat(brainState.periods);
% [~, sortIdx] = sort(brainState.periods(:, 1), 'ascend');
% brainState.periods = brainState.periods(sortIdx, :);
% brainState.names = {'NREM'; 'REM'; 'WAKE'; 'QW'};
% 
% allPeriods = brainState.periods;
% thetaPeriods = allPeriods(allPeriods(:, 3)== 3, 1:2);


% (4)

% fileName = fullfile(basePath, 'PopulationBurstEvents', [sessionName '.mua.evt']);
% 
% allEvents = LoadEvt(fileName, 1e6, {'beg'; 'end'});
% allEvents = allEvents./1e6;
% 
% nEvents = numel(allEvents)/2; 
% 
% allEvents = reshape(allEvents, [2, nEvents]);
% allEvents = allEvents';
% 
% remaze_theta = allEvents(allEvents(:, 1) > remaze(1), :);

%%%




% 
% min_interBout = 0.5;
% temp = remaze_theta(1,:);
% finalremaze_theta = [];
% 
% for ii = 2:size(remaze_theta,1)
%     if (remaze_theta(ii,1)-temp(2)) < min_interBout
%         temp = [temp(1) remaze_theta(ii, 2)]; % merge two PBEs
%     else
%         finalremaze_theta = [finalremaze_theta; temp];
%         temp = remaze_theta(ii,:);
%     end
% end
% 
% finalremaze_theta = [finalremaze_theta; temp]; 

% remaze_theta = finalremaze_theta;

idx = find(theta_periods(:, 1) > maze(1) & theta_periods(:, 1) < maze(2));
remaze_theta  = theta_periods(idx, :);

duration = diff(remaze_theta');
% remaze_theta(duration > 20, :) = [];


% time-binned firing rate of units in the re-maze periods 

binDur = 0.02;
binOverlapRatio = 0;
gw = gausswindow(3,9);

temp = timeBinning_withOverlap(remaze_theta, spikes, binDur, binOverlapRatio, fileInfo); % the theta periods
binned_remaze = temp(:, 2); % we need only data with coarse time bins


ifshuffle = 0;
clusterQuality = []; %%%%%%%%%%%%%%%%
remaze_posteriors = calculateAssemblyTuning_theta(binned_remaze, spikes, clusterQuality, binDur, ifshuffle);


remaze_posteriors_cnct = cell2mat(remaze_posteriors');
binned_remaze_cnct  = cell2mat(binned_remaze');
nTimeBins = size(binned_remaze_cnct, 2);

learnedTunings = nan(nUnits, nPosBins);

for iUnit = 1:nUnits

    
    otherUnitsFiring  = binned_remaze_cnct(setdiff(1:nUnits, iUnit), :);
    otherUnitsFiring  = otherUnitsFiring > 0;
    nOtherUnitsFiring = sum(otherUnitsFiring, 1);
    
%     idx = find(nOtherUnitsFiring > 0);
    
    currUnitPosteriors  = remaze_posteriors_cnct(:, :, iUnit);
    currUnitFirings     = binned_remaze_cnct(iUnit, :);
    
    
    weightedSum = nansum(currUnitPosteriors .* repmat(currUnitFirings, [nPosBins 1]), 2)/nTimeBins;
    px = nansum(currUnitPosteriors, 2)/nTimeBins;
    
    learnedTunings(iUnit, :) = weightedSum./px;
    
    temp = conv(weightedSum./px, gw, 'same');
    learnedTunings(iUnit, :) = temp;
       
end


aa = learnedTunings;
% aa = spatialTunings;
% aa = aa ./repmat(max(aa, [], 2), [1 size(aa, 2)]);
aa = (aa - repmat(mean(aa, 2), [1 size(aa, 2)])) ./ repmat(std(aa, [], 2), [1 size(aa, 2)]);


aa = aa(PFSortIdx, :);

figure; imagesc(aa); colormap('jet') 



allCorrs = corr(learnedTunings', spatialTunings');
learnedTuningPFcorr = diag(allCorrs);



save(fullfile(basePath, 'assemblyTunings', [sessionName '.assemblyTunings_remaze_0.02binDur.mat']), 'learnedTunings', 'learnedTuningPFcorr')


end


function [sdat, timePoints] = calculateMUA(spikes) 

spikeTimes = [];
for unit = 1:numel(spikes)
    spikeTimes = [spikeTimes; spikes(unit).time];
end

spikeTimes = sort(spikeTimes, 'ascend');

binDur = 0.001;

nBins = floor(spikeTimes(end) / binDur);
binEdges = 0:binDur:spikeTimes(end);

timePoints = (1:nBins)*binDur;

spikeBinned = histc(spikeTimes, binEdges);
spikeBinned(end) = []; %% removing the spikes that matches the last edge (we need the counts between the edges)

spikeDensity = spikeBinned/binDur;

sigma     = 0.01 * 1e3; %% in milisecond since we are using 1 milisecond time bins
halfwidth = 3 * sigma ;

smoothwin    = gausswindow(sigma, halfwidth);
spikeDensity = conv(spikeDensity, smoothwin, 'same');

rateMean = mean(spikeDensity); 
rateSD = std(spikeDensity); 

sdat = (spikeDensity - rateMean)/rateSD;


end

