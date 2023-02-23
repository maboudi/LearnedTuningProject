function basicProcedure(animalID, sessionNumber)


clc; 
close all

addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/ReplayPreplayAnalyses'))
currDir = '/nfs/turbo/umms-kdiba/Kourosh/Hiro-Dataset';


% currDir = '/home/kouroshmaboudi/Documents/HMM_project/Hiro_Dataset/pooled_includingKevin';


cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3],[1 2 3], [1]}; 

mazeShape.Roy = {'linear';'linear';'linear'};
mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
mazeShape.Kevin = {'linear'};



    
rat = rats{animalID}
sessionNumber

%% session info


%%% load the data

sessionName = [rat 'Maze' num2str(sessionNumber)];
sessionName2 = [rat '-maze' num2str(sessionNumber)];

VarList = {'spikes','behavior','position','speed','basics','ripple'};


for var = 1 : length(VarList)
    load([currDir '/wake-' VarList{var} '.mat'])
end


spikes   = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['position.' sessionName]);
speed    = eval(['speed.' sessionName]);
basics   = eval(['basics.' sessionName]);

ripple       = eval(['ripple.' sessionName]);
rippleEvents = ripple.time;


%%% initiate a unified structure for the current session data and make a
%%% master folder for storing the results from different analyses

fileinfo = struct('name', sessionName2, 'animal', rat, 'xyt', [], 'xyt2', [], ...
    'tbegin', [], 'tend', [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'pix2cm', 0.3861); 

mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
mkdir(mainDir)



%%% behavior 

if strcmp(rat, 'Kevin') % there are two wake periods, we need to concatenate them
   behavior.time(2,2) = behavior.time(3,2);
   behavior.time(3,1) = behavior.time(4,1);
   behavior.time(3,2) = behavior.time(4,2);
   behavior.time(4,:) = [];
end


behavior.time(2,1) = behavior.time(2,1) + 1e6; % slight modification to make sure this time period is only limited to the maze (1e6 = 1 sec)
behavior.time(2,2) = behavior.time(2,2) - 1e6;


bvrTimeList = behavior.list(:,[1 2]); % start and end of each behavioral state
bvrState    = behavior.list(:,3); % the label of corresponding states



%%% recording configurations 

% Fs = basics.SampleRate;
fileinfo.Fs = 1e6; % Since in Hiro's data the spike timestamps are already in microsecond we don't need to use the original recording sampling rate. 
            % The behavioral timepoints are also in microsecond therefore easier to handle them.

% lfpSampleRate = basics.lfpSampleRate;

fileinfo.lfpSampleRate = 1e6; % agian we use the time unit after conversion

fileinfo.nCh = basics.nChannels;



%%% position info and preprocessings

fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 


% linearized positions (we need this specially for analyzing the L-shape, U-shape, and circular mazes)
% the information will be loaded into the xyt2 field of fileinfo structure


animalMazeShapes = eval(sprintf('mazeShape.%s', rat));
currMazeShape    = animalMazeShapes{sessionNumber};

linearPos = linearizePosition(fileinfo, behavior, currMazeShape); % click on the middle of the platforms

fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);



% % updated

direction = 'bi';
[lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 

if length(lapsStruct.RL) > length(lapsStruct.LR)
   lapsStruct.RL(1,:) = [];
   behavior.MazeEpoch(1,:) = lapsStruct.LR(1,:);
end

totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
laps = zeros(totNumLaps, 2);
laps(1:2:totNumLaps, :)  = lapsStruct.LR;
laps(2:2:totNumLaps, :)  = lapsStruct.RL;

laps(:, 3) = 1:size(laps, 1); 


fileinfo.xyt2(:, 3) = zeros(size(fileinfo.xyt2(:, 1))); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)
   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);         
end


% calcualting the speed threshold 
% it's usually is 10 cm/s but for some animals it seems to be higher if we
% look at the distribution of speed across the whole maze duration

% runSpeedThresh = multiModalDist(speed.v(speed.t > behavior.time(2,1) & speed.t < behavior.time(2,2)), 2);
runSpeedThresh = 10; % cm/s


%% formating the spike info
% The final format is similar to what Kamran had for his 2006 datasets


unitTypes = 'all';

[spikeStruct, okUnits] = spikeBehaviorAnalysis(spikes, laps, rippleEvents, speed, unitTypes, fileinfo);

temp = [spikes.id];
shanks = temp(2*okUnits - 1);
clus   = temp(2*okUnits);


% save('allVariables.mat', 'okUnits', 'shanks', '-append')


%% Place Fields

close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)


%%% 1D spatial tuning: using linearized position

% the place fields are calculated separately for the left and right
% direction

[spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
[spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileinfo, subfolder);


[spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);


% save(fullfile(subfolder, 'biDirectional.mat'), 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir')


%% PBEs during different behavioral periods


% Detection configurations

time_resolution = 0.001; % in second
threshZ         = 3; % sdf with 3 std deviation above the mean


fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend   = behavior.time(3,2);


subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, qwake=3, wake=4


velocityFilter = 1;
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);


qclus = [1 2 3];

[binnedPBEs, secondaryPBEs] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileinfo); 
PBErippleIdx = ifContainRipples(secondaryPBEs, rippleEvents);

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
RUNidx              = find(secondaryPBEs(:, 1) > behavior.time(2,1) & secondaryPBEs(:, 2) < behavior.time(2,2) & secondaryPBEs(:,5) == 1);
RUNbinnedPBEs.data  = binnedPBEs(RUNidx, :);
secondaryPBEs_RUN   = secondaryPBEs(RUNidx, :);
RUNbinnedPBEs       = genSurrogates(RUNbinnedPBEs);


POSTbinnedPBEs       = baseStruct;
POSTidx              = find(secondaryPBEs(:, 1) > behavior.time(3,1) & secondaryPBEs(:, 2) < behavior.time(3,2));
POSTbinnedPBEs.data  = binnedPBEs(POSTidx, :);
secondaryPBEs_POST   = secondaryPBEs(POSTidx, :);
POSTbinnedPBEs       = genSurrogates(POSTbinnedPBEs);


% save(fullfile(mainDir, 'toAddVariables.mat'), 'fileinfo', 'behavior', 'secondaryPBEs')

%%  Bayesian decoding

% BayesianReplayDetection_nov2020
BayesianReplayDetection_GrosmarkReclu_nov2020


end













