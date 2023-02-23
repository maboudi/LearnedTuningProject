

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



for sessionNumber = 4%1:nSessions


sessionName = sessionNames{sessionNumber};


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['Position.' sessionName]);

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity


%% sessioninfo

unsortedData = fullfile(unsortedDataDir, sessionName, sessionName);

fileInfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 

fileInfo.tbegin = behavior.time(1,1); 
fileInfo.tend = behavior.time(3,2);




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

threshSD = 2;
ripplePeriods = RippleDetect(unsortedData, fileInfo, threshSD); % the results in seconds


fileInfo.pix2cm = 1;


% in RippleDetect2, ripples are detected separately on each channel and at
% the end were merged, peak ripple amplitude was calculated again within
% the merged event (maximum across all channels)


threshSD = 2;
velocityFilter = 1;

fileInfo.RippleChannels = best_channels(4:end); % the two last are from two extra set of channels outside of the hippocampus. The actual number of shanks is 12 (6 each hemisphere)
[ripplePeriods, rippleLFP, ripplePower] = RippleDetect2(unsortedData, fileInfo, threshSD, 0, velocityFilter);


% Calculated ripple periods are in seconds
% The rippleLFP and ripplePower are in original sampling rate. They might be easier to work with if their sampling periods are changed to 1 millisecond. 

Par = LoadPar([unsortedData '.xml']); 
lfpSampleRate = Par.lfpSampleRate; % sampling rate of eeg file

samplePer = 1/lfpSampleRate;
samplingTimes = fileInfo.tbegin:samplePer:fileInfo.tend;
samplingTimes(1) = [];

newSamplingTimes = fileInfo.tbegin:1e-3:fileInfo.tend;
newSamplingTimes(1) = [];


rippleLFP   = interp1(samplingTimes, rippleLFP, newSamplingTimes);
ripplePower = interp1(samplingTimes, ripplePower, newSamplingTimes);


%%

figure; 
hold on

for ii = 1:numel(fileInfo.RippleChannels)   
    plot(14-ii+1+rippleLFP(1:1e5, ii)/max(rippleLFP(1:1e5, ii)), 'color', [0.7 0.7 0.7], 'linewidth', 1)
    plot(14-ii+1+ripplePower(1:1e5, ii)/max(ripplePower(1:1e5, ii)), 'k', 'linewidth', 2)
end

rr = floor(ripplePeriods(1:960,1:2)'*1e3);
rr = [rr; nan(1,length(rr))];
rr = rr(:);

tt = 14+1+0.05*ones(960, 1);
tt = repmat(tt, 3,1);
tt = tt(:);
% 
% 
% 
% rr2 = floor(ripplePeriods2(1:60,1:2)'*1e3);
% rr2 = [rr2; nan(1,length(rr2))];
% rr2 = rr2(:);
% 
% tt2 = 12+1+0.1*ones(60, 1);
% tt2 = repmat(tt2, 3,1);
% tt2 = tt2(:);


hold on
plot(rr, tt, 'r', 'linewidth', 2)
% plot(rr2, tt2, 'b', 'linewidth', 2)






%% PBEs

%% defining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileInfo.name '_congruence'], '-dpng')
threshZ = 2; % number of S.D.s above the mean



% Theta periods for exclusion

% baseFile = fullfile(unsortedDataDir, sessionName, sessionName);
% 
% thetaPeriods = load([baseFile '.theta.1']);
% thetaPeriods = floor(thetaPeriods./fileInfo.lfpSampleRate);


fileInfo.lfpSampleRate = 1;


subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4



% first round of detection without any filtering due to participation of
% pyramidal units

velocityFilter = 1;
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileInfo, [], time_resolution, threshZ, exclude, velocityFilter);


% PBEs filtered based on the firing of pyramidal units and divided into
% smaller time bins during which the number of spikes of each unit is
% counted

qclus = [1 2 3]; % pyramidal units
[binnedPBEs, PBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileInfo); 



% Overlap with ripples

% check if the PBEs overlap with any of the detected ripples

PBErippleIdx = ifContainRipples(PBEs, ripplePeriods);
PBEs(:, 5)   = PBErippleIdx;
nPBEs = size(PBEs, 1);



%%
currPBEinfo = struct('sessionName', sessionName, 'ipbe', 0, 'startT', 0, 'endT', 0, 'peakT', 0, ...
    'peakMUA', 0, 'rippleOverlap', 0, 'peakRippleAmp', 0, 'PRE', 0, 'RUN', 0, 'POST', 0, 'NREM', 0, 'REM', 0, 'QW', 0, 'AW', 0, ...
    'nFiringUnits', 0, 'duration', 0, 'spkCount_1msbin', [], 'spkCount_20msbin', []);

currPBEinfo = repmat(currPBEinfo, 1, nPBEs);

for ipbe = 1:nPBEs
 
    currPBEinfo(ipbe).ipbe    = ipbe;
    
    % PBE timing
    currPBEinfo(ipbe).startT  = PBEs(ipbe, 1);
    currPBEinfo(ipbe).endT    = PBEs(ipbe, 2);
    currPBEinfo(ipbe).peakT   = PBEs(ipbe, 3);
    
    
    % peak MUA amplitude
    currPBEinfo(ipbe).peakMUA = PBEs(ipbe, 4);
    

    % overlap with detected ripples
    currPBEinfo(ipbe).rippleOverlap = PBEs(ipbe, 5);
    
    
    % max ripple amplitude within a PBE
    
    % an alternative is to calculate max ripple amplitude during the PBEs (maximum across channels and time)

    temp = ripplePower(floor(currPBEinfo(ipbe).startT*1e3): floor(currPBEinfo(ipbe).endT*1e3), :);
    currPBEinfo(ipbe).peakRippleAmp = max(temp(:));  
    
 
    % behavioral epoch of PBE
    
    currPBEinfo(ipbe).PRE  = 0;
    currPBEinfo(ipbe).RUN  = 0;
    currPBEinfo(ipbe).POST = 0;
    
    if PBEs(:, 1) > behavior.time(1,1) & PBEs(:, 2) < behavior.time(1,2); currPBEinfo(ipbe).PRE = 1; end
    if PBEs(:, 1) > behavior.time(2,1) & PBEs(:, 2) < behavior.time(2,2); currPBEinfo(ipbe).RUN = 1; end
    if PBEs(:, 1) > behavior.time(3,1) & PBEs(:, 2) < behavior.time(3,2); currPBEinfo(ipbe).POST = 1; end
    
    
    % brain state during which a PBE occurred
    boutInd        = find(bvrTimeList(:,1) < currPBEinfo(ipbe).peakT & bvrTimeList(:,2) > currPBEinfo(ipbe).peakT, 1, 'first');
    boutBrainState = bvrState(boutInd);
    
    
    currPBEinfo(ipbe).NREM = 0;
    currPBEinfo(ipbe).REM  = 0;
    currPBEinfo(ipbe).QW   = 0;
    currPBEinfo(ipbe).AW   = 0;
    
    if boutBrainState == 1; currPBEinfo(ipbe).NREM = 1; end
    if boutBrainState == 2; currPBEinfo(ipbe).REM = 1; end
    if boutBrainState == 3; currPBEinfo(ipbe).QW = 1; end
    if boutBrainState == 4; currPBEinfo(ipbe).AW = 1; end
    

    % binned spikes counts during each PBE
    currPBEinfo(ipbe).spkCount_1msbin  = binnedPBEs{ipbe, 1};
    currPBEinfo(ipbe).spkCount_20msbin = binnedPBEs{ipbe, 2};
    
 
    % number of participating units during eah PBE and number of 
    currPBEinfo(ipbe).nFiringUnits   = numel(find(sum(binnedPBEs{ipbe, 2}, 2)));
    currPBEinfo(ipbe).duration       = size(binnedPBEs{ipbe, 1}, 2);
    
end



%% Ripples

nRipples = size(ripplePeriods, 1);
currRippleInfo = struct('sessionName', sessionName, 'iRipp', 0, 'startT', 0, 'endT', 0, 'peakT', 0, 'duration', 0, 'peakRippleAmp', 0, 'peakMUA', 0);
currRippleInfo  = repmat(currRippleInfo , 1, nRipples);

for iRipp = 1:nRipples
    
    currRippleInfo(iRipp).iRipp  = iRipp;
    currRippleInfo(iRipp).startT = ripplePeriods(iRipp, 1);
    currRippleInfo(iRipp).endT   = ripplePeriods(iRipp, 2);
    
    currRippleInfo(iRipp).duration   = (ripplePeriods(iRipp, 2) - ripplePeriods(iRipp, 1))*1000;
    
    currRippleInfo(iRipp).peakT  = ripplePeriods(iRipp, 3);
    
    currRippleInfo(iRipp).peakRippleAmp = ripplePeriods(iRipp, 4);
    
    currRippleInfo(iRipp).peakMUA = max(sdat(floor(ripplePeriods(iRipp, 1)*1e3):floor(ripplePeriods(iRipp, 2)*1e3)));
end



%% generate PBEinfo

sessionInfoKM.(sessionName).behavior = behavior;

sessionInfoKM.(sessionName).spikes   = spikeStruct;
sessionInfoKM.(sessionName).fileInfo = fileInfo;
sessionInfoKM.(sessionName).shankIDs = shanks;
sessionInfoKM.(sessionName).cluIDs   = clus;


sessionInfoKM.(sessionName).ripples     = currRippleInfo;
sessionInfoKM.(sessionName).rippleLFP   = rippleLFP;
sessionInfoKM.(sessionName).ripplePower = ripplePower;


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


%% low and high participation ripples across time 

rippleAmp_at_ripple = [currRippleInfo.peakRippleAmp];
MUA_at_ripple       = [currRippleInfo.peakMUA];
ripplePeak = [currRippleInfo.peakT];


low_MUA_ripples  = find(rippleAmp_at_ripple > 2 & MUA_at_ripple < 3);
low_MUA_ripples_peakT = ripplePeak(low_MUA_ripples);

high_MUA_ripples = find(rippleAmp_at_ripple > 2 & MUA_at_ripple > 3);
high_MUA_ripples_peakT = ripplePeak(high_MUA_ripples);



bins = fileInfo.tbegin:5*60:fileInfo.tend;
binCenters = bins(1:end-1) + diff(bins(1:2))/2;

h_low = histc(low_MUA_ripples_peakT, bins); h_low(end) = [];
h_hi  = histc(high_MUA_ripples_peakT, bins); h_hi(end) = [];

h_low = h_low/(5*60);
h_hi  = h_hi/(5*60);


meanPRE_low  = mean(h_low(binCenters < behavior.MazeEpoch(1)));
meanPOST_low = mean(h_low(binCenters > behavior.MazeEpoch(2)));

meanPRE_high  = mean(h_hi(binCenters < behavior.MazeEpoch(1)));
meanPOST_high = mean(h_hi(binCenters > behavior.MazeEpoch(2)));



figure; 
ax(1) = subplot(2,1,1);
set(gca, 'fontsize', 10)

plot(binCenters/3600, h_low, 'linewidth', 2, 'color', 'k')
line([behavior.MazeEpoch(1) behavior.MazeEpoch(1)]/3600, [0 max(h_low)], 'color', 'r', 'linewidth', 2)
line([behavior.MazeEpoch(2) behavior.MazeEpoch(2)]/3600, [0 max(h_low)], 'color', 'r', 'linewidth', 2)

line([fileInfo.tbegin behavior.MazeEpoch(1)]/3600, [meanPRE_low meanPRE_low], 'color', 'g', 'linewidth', 2)
line([behavior.MazeEpoch(2) fileInfo.tend]/3600, [meanPOST_low meanPOST_low], 'color', 'g', 'linewidth', 2)

ylabel('rate(Hz)', 'fontsize', 12)
title('ripple events (>2z) with low MUA (<3z)')



ax(2) = subplot(2,1,2);
set(gca, 'fontsize', 10)


plot(binCenters/3600, h_hi, 'linewidth', 2, 'color', 'k')
line([behavior.MazeEpoch(1) behavior.MazeEpoch(1)]/3600, [0 max(h_hi)], 'color', 'r', 'linewidth', 2)
line([behavior.MazeEpoch(2) behavior.MazeEpoch(2)]/3600, [0 max(h_hi)], 'color', 'r', 'linewidth', 2)

line([fileInfo.tbegin behavior.MazeEpoch(1)]/3600, [meanPRE_high meanPRE_high], 'color', 'g', 'linewidth', 2)
line([behavior.MazeEpoch(2) fileInfo.tend]/3600, [meanPOST_high meanPOST_high], 'color', 'g', 'linewidth', 2)

linkaxes(ax, 'x')

xlabel('time(hr)', 'fontsize', 12)
ylabel('rate(Hz)', 'fontsize', 12)
title('ripple events (>2z) with high MUA (>3z)')


%% PBEs with low and high ripple power across sleep time

rippleAmp_at_PBE    = [currPBEinfo.peakRippleAmp];
MUA_at_PBE          = [currPBEinfo.peakMUA];
PBEpeakT            = [currPBEinfo.peakT];


low_ripple_PBEs       = find(rippleAmp_at_PBE < 3);
low_ripple_PBEs_peakT = PBEpeakT(low_ripple_PBEs);

hi_ripple_PBEs        = find(rippleAmp_at_PBE > 3);
hi_ripple_PBEs_peakT  = PBEpeakT(hi_ripple_PBEs);


bins = fileInfo.tbegin:5*60:fileInfo.tend;
binCenters = bins(1:end-1) + diff(bins(1:2))/2;

h_low = histc(low_ripple_PBEs_peakT, bins); h_low(end) = [];
h_hi  = histc(hi_ripple_PBEs_peakT, bins); h_hi(end) = [];

h_low = h_low/(5*60);
h_hi  = h_hi/(5*60);


meanPRE_low  = mean(h_low(binCenters < behavior.MazeEpoch(1)));
meanPOST_low = mean(h_low(binCenters > behavior.MazeEpoch(2)));

meanPRE_high  = mean(h_hi(binCenters < behavior.MazeEpoch(1)));
meanPOST_high = mean(h_hi(binCenters > behavior.MazeEpoch(2)));



figure; 
ax(1) = subplot(2,1,1);
set(gca, 'fontsize', 10)

plot(binCenters/3600, h_low, 'linewidth', 2, 'color', 'k')
line([behavior.MazeEpoch(1) behavior.MazeEpoch(1)]/3600, [0 max(h_low)], 'color', 'r', 'linewidth', 2)
line([behavior.MazeEpoch(2) behavior.MazeEpoch(2)]/3600, [0 max(h_low)], 'color', 'r', 'linewidth', 2)

line([fileInfo.tbegin behavior.MazeEpoch(1)]/3600, [meanPRE_low meanPRE_low], 'color', 'g', 'linewidth', 2)
line([behavior.MazeEpoch(2) fileInfo.tend]/3600, [meanPOST_low meanPOST_low], 'color', 'g', 'linewidth', 2)

ylabel('rate(Hz)', 'fontsize', 12)
title('PBEs (>2z) with low ripple amplitude (<3z)')


ax(2) = subplot(2,1,2);
set(gca, 'fontsize', 10)

plot(binCenters/3600, h_hi, 'linewidth', 2, 'color', 'k')
line([behavior.MazeEpoch(1) behavior.MazeEpoch(1)]/3600, [0 max(h_hi)], 'color', 'r', 'linewidth', 2)
line([behavior.MazeEpoch(2) behavior.MazeEpoch(2)]/3600, [0 max(h_hi)], 'color', 'r', 'linewidth', 2)

line([fileInfo.tbegin behavior.MazeEpoch(1)]/3600, [meanPRE_high meanPRE_high], 'color', 'g', 'linewidth', 2)
line([behavior.MazeEpoch(2) fileInfo.tend]/3600, [meanPOST_high meanPOST_high], 'color', 'g', 'linewidth', 2)

linkaxes(ax, 'x')

xlabel('time(hr)', 'fontsize', 12)
ylabel('rate(Hz)', 'fontsize', 12)
title('PBEs (>2z) with high ripple amplitude (>3z)')


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