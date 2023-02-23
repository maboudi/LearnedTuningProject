

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/ReplayPreplayAnalyses'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_BapunReclustered';
cd(currDir)

unsortedDataDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_originalClusters';



%% loading session data

VarList = {'spikes','behavior','position'};

% note that the spikes within the .mat file are from just the pyramidal
% units, for MUA I need to load the clu and res files again


% Although no difference regarding the behavior and position data

for var = 1 : length(VarList)
    load([currDir '/CRCNSReclustered-' VarList{var} '.mat'])
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



for sessionNumber = 1:nSessions


sessionName = sessionNames{sessionNumber};


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['Position.' sessionName]);

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity


%% sessioninfo

unsortedData = fullfile(unsortedDataDir, sessionName, sessionName);

fileInfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileInfo.name);
mkdir(mainDir)



Par = LoadXml(unsortedData); % Load the recording configurations from the xml file

fileInfo.nCh = Par.nChannels;


fileInfo.Fs  = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording
fileInfo.lfpSampleRate = 1;




% RUN high power theta periods

thetaPeriods = importEVTfile(unsortedData); % make sure that we know details about procedures and parameters which were used for the detection of these periods 
% thetaPeriods = [];



% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

% 1D pos
% positon1D = position.OneDLocation* cov2cmfac;

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



fileInfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track


fileInfo.xyt(xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1) | ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3), 1:2) = NaN;


% xpos = fileInfo.xyt(:, 1);
% ypos = fileInfo.xyt(:, 2);

% figure; plot(xpos, ypos, '.', 'markersize', 3)


linearPos = linearizePosition(fileInfo, behavior, mazeShapes{sessionNumber});

fileInfo.xyt2(:, 1) = linearPos; 

% fileInfo.xyt2(:, 1) = positon1D;
fileInfo.xyt2(:, 2) = fileInfo.xyt(:, 3);


% laps = calculateLapTimings(fileInfo, speed, mainDir);

if ismember(sessionNumber, [2 5]) 
    direction = 'uni'; 
    [laps, turningPeriods] = calculateLapTimings(fileInfo, speed, direction, mainDir); 
else
    direction = 'bi';
    [lapsStruct, turningPeriods] = calculateLapTimings(fileInfo, speed, direction, mainDir); 
    
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


fileInfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileInfo.xyt2(:, 2) > laps(ii, 1) & fileInfo.xyt2(:, 2) < laps(ii, 2));
   fileInfo.xyt2(idx, 3) = laps(ii, 3);
          
end

runSpeedThresh = 10; % cm/s


%% Behavioral states


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Should we consider them as
% QW or NREM?


bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ...
            3 *   ones(length(behavior.Drowsy), 1); ... % considering Drowsy equivalent to QW
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering Intermediate as NREM
            2   *   ones(length(behavior.REM), 1)];


%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

[spikeStruct,okUnits] = spikeBehaviorAnalysis2(spikes, speed, laps, thetaPeriods, qual2consider, fileInfo);

temp = [spikes.id];
shanks = temp(2*okUnits - 1);
clus   = temp(2*okUnits);

% %

% save('allVariables.mat', 'okUnits', 'shanks', '-append')

% %


% We need the multiunit (any detected spike) as well for defining the population burst events

numberofShanks = length(dir([unsortedData '.res.*']));


MUA.t = [];
for shank = 1:numberofShanks
    
    resContent  = load([unsortedData '.res.' num2str(shank)]);
    
    cluContent  = load([unsortedData '.clu.' num2str(shank)]);
    clusterIDs  = cluContent(2:end);
    
    spikes2Add = resContent(clusterIDs ~= 0);
    
    MUA.t = [MUA.t; spikes2Add]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');
% The original smapling frequency is 20000

MUA.t = MUA.t./20000;



%% Spatial tuning of the units

% close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)

%%% 1D spatial tuning: using linearized position

% the place fields are calculated separately for the left and right
% direction


if ismember(sessionNumber, [2 5]) % sessions 2, 5, and 8 are circular

    [spatialTunings, PF_sorted, runTemplate, spatialInfo, conslapsRatio, diffWithAvg] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileInfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileInfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings, [], spatialInfo, [], conslapsRatio, [], behavior, [], fileInfo, subfolder);
else

    [spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileInfo, behavior, [], [], speed, 'LR', 2, runSpeedThresh, [], fileInfo.Fs, subfolder);
    [spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileInfo, behavior, [], [], speed, 'RL', 2, runSpeedThresh, [], fileInfo.Fs, subfolder);
    
    [spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileInfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileInfo.Fs, subfolder);

%     firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileInfo, subfolder);
end


close all




%% Ripple Detection

% After visualizing deteceted events using different metrics (MUA and
% ripple ampplitude) on Neuroscope, I realized that ripple is more reliable
 
best_channels = BigRippleChannels(unsortedData, fileInfo);
fileInfo.RippleChannels = best_channels(1);


threshSD = 3;
[ripplePeriods, rippleLFP, ripplePower, ripplePower_orig] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds

% %
threshSD = 2;
fileInfo.RippleChannels = best_channels;
[ripplePeriods_1, rippleLFP_1, ripplePower_1] = RippleDetect2(unsortedData, fileInfo,threshSD);

% %

sampleDur = 1/1250;
samplingTimes = fileInfo.tbegin:sampleDur:fileInfo.tend;
samplingTimes(1) = [];

newSamplingTimes = fileInfo.tbegin:1e-3:fileInfo.tend;
newSamplingTimes(1) = [];


rippleLFP   = interp1(samplingTimes, rippleLFP, newSamplingTimes);
ripplePower = interp1(samplingTimes, ripplePower, newSamplingTimes);
ripplePower_orig = interp1(samplingTimes, ripplePower_orig, newSamplingTimes);


ripplePower_2 = zeros(numel(newSamplingTimes), numel(best_channels));
for ii = 1:numel(best_channels)
    
    ripplePower_2(:, ii) = interp1(samplingTimes, ripplePower_1(:, ii), newSamplingTimes);
    
    
end


%% a little experiment with ripples on different channel

threshSD = 3;





fileInfo.RippleChannels = best_channels(1);

[ripplePeriods, rippleLFP, ripplePower, ripplePower_orig] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds



fileInfo.RippleChannels = best_channels(2);

[ripplePeriods2, rippleLFP2, ripplePower2, ripplePower_orig2] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds


fileInfo.RippleChannels = best_channels(3);

[ripplePeriods3, rippleLFP3, ripplePower3, ripplePower_orig3] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds


fileInfo.RippleChannels = best_channels(4);

[ripplePeriods4, rippleLFP4, ripplePower4, ripplePower_orig4] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds


fileInfo.RippleChannels = best_channels(5);

[ripplePeriods5, rippleLFP5, ripplePower5, ripplePower_orig5] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds




fileInfo.RippleChannels = best_channels(1:5);

[ripplePeriodsT, rippleLFPT, ripplePowerT, ripplePower_origT] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds


fileInfo.RippleChannels = best_channels;

[ripplePeriodsT2, rippleLFPT2, ripplePowerT2, ripplePower_origT2] = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds

rippleLFPT2   = interp1(samplingTimes, rippleLFPT2, newSamplingTimes);
ripplePowerT2 = interp1(samplingTimes, ripplePowerT2, newSamplingTimes);
ripplePower_origT2 = interp1(samplingTimes, ripplePower_origT2, newSamplingTimes);


rippleLFP2   = interp1(samplingTimes, rippleLFP2, newSamplingTimes);
ripplePower2 = interp1(samplingTimes, ripplePower2, newSamplingTimes);
ripplePower_orig2 = interp1(samplingTimes, ripplePower_orig2, newSamplingTimes);




figure; 

ax(1) = subplot(7,1,1);
plot(samplingTimes(1:125000), ripplePower(1:125000), 'k')
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')

ax(2) = subplot(7,1,2);
plot(samplingTimes(1:125000), ripplePower2(1:125000), 'k')
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')

ax(3) = subplot(7,1,3);
plot(samplingTimes(1:125000), ripplePower3(1:125000), 'k')
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')

ax(4) = subplot(7,1,4);
plot(samplingTimes(1:125000), ripplePower4(1:125000), 'k')
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')

ax(5) = subplot(7,1,5);
plot(samplingTimes(1:125000), ripplePower5(1:125000), 'k')
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')


ax(6) = subplot(7,1,6);
plot(samplingTimes(1:125000), ripplePowerT(1:125000), 'r', 'linewidth', 2)
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')


ax(7) = subplot(7,1,7);
plot(samplingTimes(1:125000), ripplePowerT2(1:125000), 'g', 'linewidth', 2)
yticks(1:7)
set(gca, 'YGrid', 'on', 'XGrid', 'off')


linkaxes(ax, 'xy')
% xlim([0 100])





%% PBEs

%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileInfo.name '_congruence'], '-dpng')
threshZ = 2; % sdf with 3 std deviation above the mean


% baseFile = fullfile(unsortedDataDir, sessionName, sessionName);
% 
% thetaPeriods = load([baseFile '.theta.1']);
% thetaPeriods = floor(thetaPeriods./fileInfo.lfpSampleRate);


fileInfo.lfpSampleRate = 1;

fileInfo.tbegin = behavior.time(1,1); 
fileInfo.tend = behavior.time(3,2);


subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4


velocityFilter = 1;
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileInfo, [], time_resolution, threshZ, exclude, velocityFilter);


qclus = [1 2 3];

[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileInfo); 
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




%% generate PBEinfo

nPBEs = size(binnedPBEs, 1);

for ipbe = 1:nPBEs
 
    currPBEinfo(ipbe).sessionName = sessionName;
    currPBEinfo(ipbe).ipbe       = ipbe;
    
    currPBEinfo(ipbe).startT  = secondaryPBEs(ipbe, 1);
    currPBEinfo(ipbe).endT    = secondaryPBEs(ipbe, 2);
    currPBEinfo(ipbe).peakT   = secondaryPBEs(ipbe, 3);
    currPBEinfo(ipbe).peakMUA = secondaryPBEs(ipbe, 4);
    
    currPBEinfo(ipbe).rippleOverlap = secondaryPBEs(ipbe, 5);
    
    currPBEinfo(ipbe).NREM    = secondaryPBEs(ipbe, 6);
    currPBEinfo(ipbe).REM     = secondaryPBEs(ipbe, 7);
    currPBEinfo(ipbe).QW      = secondaryPBEs(ipbe, 8);
    currPBEinfo(ipbe).AW      = secondaryPBEs(ipbe, 9);
    
    
    if ismember(ipbe, PREidx); currPBEinfo(ipbe).PRE = 1;  currPBEinfo(ipbe).RUN = 0; currPBEinfo(ipbe).POST = 0; end
    if ismember(ipbe, RUNidx); currPBEinfo(ipbe).PRE = 0;  currPBEinfo(ipbe).RUN = 1; currPBEinfo(ipbe).POST = 0; end
    if ismember(ipbe, POSTidx); currPBEinfo(ipbe).PRE = 0;  currPBEinfo(ipbe).RUN = 0; currPBEinfo(ipbe).POST = 1; end
    
    
    currPBEinfo(ipbe).firing_1msbin  = binnedPBEs{ipbe, 1};
    currPBEinfo(ipbe).firing_20msbin = binnedPBEs{ipbe, 2};
    
    currPBEinfo(ipbe).nFiringUnits   = numel(find(sum(binnedPBEs{ipbe, 2}, 2)));
    currPBEinfo(ipbe).duration       = size(binnedPBEs{ipbe, 1}, 2);
    
end



sessionInfoKM.(sessionName).behavior = behavior;

sessionInfoKM.(sessionName).spikes   = spikeStruct;
sessionInfoKM.(sessionName).fileInfo = fileInfo;
sessionInfoKM.(sessionName).shankIDs = shanks;
sessionInfoKM.(sessionName).cluIDs   = clus;

sessionInfoKM.(sessionName).ripplePeriods = ripplePeriods;
sessionInfoKM.(sessionName).rippleLFP     = rippleLFP;
sessionInfoKM.(sessionName).ripplePower   = ripplePower;


sessionInfoKM.(sessionName).PBEs     = currPBEinfo;

sessionInfoKM.(sessionName).position = linearPos;


if ismember(sessionNumber, [2 5 8])

    sessionInfoKM.(sessionName).spatialTunings = spatialTunings;
    sessionInfoKM.(sessionName).runTemplate    = runTemplate;
    sessionInfoKM.(sessionName).conslapsRatio  = conslapsRatio;
    sessionInfoKM.(sessionName).diffWithAvg    = diffWithAvg;

else
    sessionInfoKM.(sessionName).spatialTunings.RL = spatialTunings_RL;
    sessionInfoKM.(sessionName).runTemplate.RL    = runTemplate_RL;
    sessionInfoKM.(sessionName).conslapsRatio.RL  = conslapsRatio_RL;
    sessionInfoKM.(sessionName).diffWithAvg.RL    = diffWithAvg_RL;
    
    sessionInfoKM.(sessionName).spatialTunings.LR = spatialTunings_LR;
    sessionInfoKM.(sessionName).runTemplate.LR    = runTemplate_LR;
    sessionInfoKM.(sessionName).conslapsRatio.LR  = conslapsRatio_LR;
    sessionInfoKM.(sessionName).diffWithAvg.LR    = diffWithAvg_LR;
    
    sessionInfoKM.(sessionName).spatialTunings.biDir = spatialTunings_biDir;
    sessionInfoKM.(sessionName).runTemplate.biDir    = runTemplate_biDir;
    sessionInfoKM.(sessionName).conslapsRatio.biDir  = conslapsRatio_biDir;
    sessionInfoKM.(sessionName).diffWithAvg.biDir    = diffWithAvg_biDir;
    
end


%% sdat at peak ripple amplitudes and ripple amplitude at peak sdat


ripplePeakTimes = floor(ripplePeriods_1(2:end, 3)*1000);
rippleSurr      = zeros(numel(ripplePeakTimes), 200);

for jj = 1:100%numel(ripplePeakTimes)
    rippleSurr(jj, :) = ripplePeakTimes(jj)-99:ripplePeakTimes(jj)+100;
end

rippleSurr = rippleSurr(1:100, :);

rippleA_around_peakMUA = ripplePower(rippleSurr);





peakMUAtimes = floor([currPBEinfo.peakT]'*1000);
peakMUA      = [currPBEinfo.peakMUA]';
nFiringUnits = [currPBEinfo.nFiringUnits]';


% % calculate ripple amplitude during PBEs

% we can take the peak MUA within each PBE

rippleA_at_peakMUA = max(ripplePower_2(peakMUAtimes, :), [], 2);

% or calculate the maximum ripple amplitude within the PBE boundaries

nPBEs = length(currPBEinfo);
rippleA_at_peakMUA = zeros(nPBEs, 1);

for ii = 1: nPBEs
   
    currMat = ripplePower_2(floor(currPBEinfo(ii).startT*1000): floor(currPBEinfo(ii).endT*1000));
    rippleA_at_peakMUA(ii) = max(currMat(:));
    
end



figure; 
hold on
scatter(rippleA_at_peakMUA, peakMUA, 5, 'filled')

xlabel('ripple amp (z)', 'fontsize', 12)
ylabel('MUA (z)', 'fontsize', 12)
title('MUA - rippleA relationship at MUA > 3 S.D.', 'fontsize', 14)

line([3 3], [2 max(peakMUA)], 'color','r', 'linewidth', 2)
line([min(rippleA_at_peakMUA) max(rippleA_at_peakMUA)], [3 3], 'color', 'r', 'linewidth', 2)


figure; 
hold on
scatter(rippleA_at_peakMUA, nFiringUnits, 5, 'filled')

xlabel('ripple amp (z)', 'fontsize', 12)
ylabel('no. firing units', 'fontsize', 12)
title('no. firing units - rippleA relationship at MUA > 3 S.D.', 'fontsize', 14)

line([3 3], [min(nFiringUnits) max(nFiringUnits)], 'color','r', 'linewidth', 2)
line([min(rippleA_at_peakMUA) max(rippleA_at_peakMUA)], [3 3], 'color', 'r', 'linewidth', 2)


[rho, pi] = corr(rippleA_at_peakMUA, nFiringUnits);
[rho, pi] = corr(peakMUA, nFiringUnits);


figure; scatter(peakMUA, nFiringUnits, 5, 'filled')


figure; hist(rippleA_at_peakMUA(peakMUA<3 & nFiringUnits < 5))
%

peakRippTimes = floor(ripplePeriods_1(:, 3)*1000);
peakRippA     = ripplePeriods_1(:, 4);


% % to calculate MUA at peak ripple times

MUA_at_peakRipple = sdat(peakRippTimes);

% an alternative is to consider the entire ripple event
nRipples = numel(peakRippTimes);
MUA_at_peakRipple = zeros(nRipples, 1);

for ii = 1:nRipples
    
    MUA_at_peakRipple(ii)= max(sdat(floor(ripplePeriods_1(ii, 1)*1e3):floor(ripplePeriods_1(ii, 2)*1e3)));
    
    
end


figure;
hold on
scatter(peakRippA, MUA_at_peakRipple, 5, 'filled')

xlabel('ripple amp (z)', 'fontsize', 12)
ylabel('MUA (z)', 'fontsize', 12)
title('MUA - rippleA relationship at rippleA > 3 S.D.', 'fontsize', 14)

line([3 3], [min(MUA_at_peakRipple) max(MUA_at_peakRipple)], 'color','r', 'linewidth', 2)
line([2 max(peakRippA)], [3 3], 'color', 'r', 'linewidth', 2)



[rho, pval] = corr(peakRippA, MUA_at_peakRipple); 






% 
% 
% bins = linspace(-1, 20, 100);
% 
% binCenters = bins(1:end-1) + diff(bins(1:2))/2;
% 
% h_peakMUA = hist(rippleA_at_peakMUA, binCenters);
% h_rndMUA  = hist(rippleA_at_random, binCenters);
% 
% figure; 
% hold on
% 
% bar(binCenters, h_rndMUA, 'EdgeColor', 'none', 'FaceColor', [.7 .7 .7], 'FaceAlpha', 0.5)
% bar(binCenters, h_peakMUA, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.2)
% 





% clear currPBEinfo

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