clear; 
clc; 
close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Hiro_Dataset';
cd(currDir)


%% loading session data


rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[3],[], []}; 

mazeShape.Roy = {'linear';'linear';'linear'};
mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
mazeShape.Kevin = {'linear'};



for rr = 1:length(rats) % loop over all the rats and sessions

    
rat = rats{rr}

sessionNumbers = allsessionNumbers{rr};


for sessionNumber = sessionNumbers
    
 close all
 

sessionNumber

%% session info


%%% load the data

sessionName = [rat 'Maze' num2str(sessionNumber)];
sessionName2 = [rat '-maze' num2str(sessionNumber)];


VarList = {'spikes','behavior','position','speed','basics','ripple'};


for var = 1 : length(VarList)
    load([currDir '/pooled_includingKevin/wake-' VarList{var} '.mat'])
end

spikes = eval(['spikes.' sessionName]);
behavior = eval(['behavior.' sessionName]);
position = eval(['position.' sessionName]);
speed = eval(['speed.' sessionName]);
basics = eval(['basics.' sessionName]);

ripple = eval(['ripple.' sessionName]);
rippleEvents = ripple.time;

%%% initiate a unified structure for the current session data and creating a
%%% master folder for storing the results from different analyses

fileinfo = struct('name', sessionName2, 'animal', rat, 'xyt', [], 'xyt2', [], ...
    'tbegin', [], 'tend', [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'pix2cm', 0.3861); 

mainDir = fullfile(currDir, 'analysesResults2', fileinfo.name);
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
bvrState = behavior.list(:,3); % the label of corresponding states



%%% recording configurations 

% Fs = basics.SampleRate;
fileinfo.Fs = 1e6; % Since in Hiro's data the spike timestamps are already in microsecond we don't need to use the original sampling rate. 
            %The behavioral timepoints are also in microsecond therefore easier to handle them.

% lfpSampleRate = basics.lfpSampleRate;

fileinfo.lfpSampleRate = 1e6; % agian we use the time unit after conversion

fileinfo.nCh = basics.nChannels;



%%% position info and preprocessings

fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 


% linearized positions (we need this spacially for analyzing the L-shape, U-shape, and circular mazes)
% the information will be loaded into the xyt2 field of fileinfo structure


animalMazeShapes = eval(sprintf('mazeShape.%s', rat));
currMazeShape = animalMazeShapes{sessionNumber};

linearPos = linearizePosition(fileinfo, behavior, currMazeShape);

fileinfo.xyt2(:, 1) = linearPos; 
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


% calculate timing of the laps (animal travels from one end to the other end of the track)

laps = calculateLapTimings(fileinfo, speed, mainDir);


% in most of the sessions the first lap is RL but rarely LR in which case
% we can omit the first lap

posAtStrEndofFirstLap = interp1(fileinfo.xyt2(:,2), fileinfo.xyt2(:,1), laps(1,:));

if posAtStrEndofFirstLap(1) < posAtStrEndofFirstLap(2)
    laps(1,:) = [];
end


laps = [laps (1:length(laps))'];

% % exclude the laps with long durations
% lapDur = diff(laps(:,1:2)')';
% laps(lapDur > 2*median(lapDur), :) = [];


fileinfo.xyt2(:, 3) = zeros(size(linearPos)); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end


% For ROY sessions, we need to make sure confining the rat on the platform
% doesn't make any difference to the firing characteristics of
% neurons/place cells. Otherwise, we need to split the maze part to two
% One way is Bayesian decoding cross-validation ... 


% %%% see if we want to partition the data to two halves (here based on the time gaps between the laps)
% 
% gapbtwLaps = laps(2:end, 1) - laps(1:end-1, 2);
% 
% figure; 
% plot(2:length(laps), gapbtwLaps)

% breakLap = input('\n Enter the lap number for truncating the data\n(If no partitioning is needed just press enter) \n');  


if exist('breakLap', 'var')
    firstORsecond = input('\nselect a part of data to process (one or two)\n');
else
    firstORsecond = [];
end


% calcualting the speed threshold 
% it's usually is 10 cm/s but for some animals it seems to be higher if we
% look at the distribution of speed across the whole maze duration

% runSpeedThresh = multiModalDist(speed.v(speed.t > behavior.time(2,1) & speed.t < behavior.time(2,2)), 2);
runSpeedThresh = 10; % cm/s


%% formating the spike info
% The final format is similar to what Kamran had for his 2006 datasets


unitTypes = 'all';
spikeStruct = spikeBehaviorAnalysis(spikes, laps, rippleEvents, speed, unitTypes, fileinfo);



%% Place Fields

close all

directory = fullfile(mainDir, 'PlaceFields');
mkdir(directory)



%%% 2D spatial tuning


% spatialTunings2D_LR = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, directory);
% spatialTunings2D_RL = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, directory);


%%% 1D spatial tuning (using linearized positions)

% the place fields are calculated separately for the left and right
% direction

[spatialTunings_LR, PF_sorted_LR, PF_peak_LR, activeUnits_sorted_LR, sparsity_LR] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, directory);
[spatialTunings_RL, PF_sorted_RL, PF_peak_RL, activeUnits_sorted_RL, sparsity_RL] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, directory);



% % plotting lap-by-lap firing of units: analyzing consistence of place
% % fields  across the session

firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, behavior, fileinfo, directory);


%% Determining PBEs during different behavioral periods


%%% Detemining population burst periods

% For detection of the population bursts, Kamran suggested using all the clusters
% even those inclusing MUA and low amplitude spikes. Therefore, it's better to load the spike timestamps from
% res files. I don't have access to the res files for now, so will
% continure for now with the sorted spikes


% Detection configurations

time_resolution = 0.001; % in second
threshZ         = 3; % sdf with 3 std deviation above the mean


%% PRE PBEs

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend   = behavior.time(1,2);

directory = fullfile(mainDir, 'PBEs', 'PRE');
mkdir(directory)


filename = fullfile(directory, [fileinfo.name '-PBEs']);

if isfile(filename)
   load(filename)
else

% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4
[primaryPBEs_pre, sdat_pre] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, 0);


PREbinnedPBEs = struct('data', [], 'p', [], 'ts', [], 'pts', []);
[PREbinnedPBEs.data, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, spikeStruct, fileinfo);
savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, fileinfo, directory)

end


PREbinnedPBEs = genSurrogates(PREbinnedPBEs); % generate poission and within-PBE-time-swap surrogate datasets

plotSurroundMultiunitRate(sdat_pre, secondaryPBEs_pre, fileinfo, directory)



%% RUN PBEs

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);

directory = fullfile(mainDir, 'PBEs', 'RUN');
mkdir(directory)


filename = fullfile(directory, [fileinfo.name '-PBEs']);

if isfile(filename)
   load(filename)
else
    
% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4
[primaryPBEs_run, sdat_run] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, 0);


RUNbinnedPBEs = struct('data', [], 'p', [], 'ts', [], 'pts', []);
[RUNbinnedPBEs.data, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, spikeStruct, fileinfo);
savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, fileinfo, directory)

end

RUNbinnedPBEs = genSurrogates(RUNbinnedPBEs); % generate poission and within-PBE-time-swap surrogate datasets

plotSurroundMultiunitRate(sdat_run, secondaryPBEs_run, fileinfo, directory)



%% POST PBEs


fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);

directory = fullfile(mainDir, 'PBEs', 'POST');
mkdir(directory)


filename = fullfile(directory, [fileinfo.name '-PBEs']);

if isfile(filename)
   load(filename)
else

% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4
[primaryPBEs_post, sdat_post] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, 0);


POSTbinnedPBEs = struct('data', [], 'p', [], 'ts', [], 'pts', []);
[POSTbinnedPBEs.data, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, spikeStruct, fileinfo);
savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, fileinfo, directory)

end

POSTbinnedPBEs = genSurrogates(POSTbinnedPBEs); % generate poission and within-PBE-time-swap surrogate datasets

plotSurroundMultiunitRate(sdat_post, secondaryPBEs_post, fileinfo, directory)


%% RUN (behavioral time scale binning)

runBouts = laps;

%%% binning the spikes within each run bout

binDur = 0.02; 
qclus = [1 2 3]; % Only the pyramidal neurons were included


aRUNbinnedBouts  = struct('data', [], 'p', [], 'ts', [], 'pts', []);

aRUNbinnedBouts.data  = timeBinning(runBouts, spikeStruct, qclus, binDur, fileinfo);
aRUNbinnedBouts.data  = removeSideSilentBins(aRUNbinnedBouts.data, runBouts(:, 1), binDur, fileinfo);

aRUNbinnedBouts = genSurrogates(aRUNbinnedBouts);



%%  Bayesian decoding

BayesianDecodingSequenceAnalysis_polar



%% Hidden Markov model

directory = fullfile(mainDir, 'HMM');
mkdir(directory)

activeRUNmodelquality

% HMMsequenceAnalysis


end

end










