clear; 
clc; 
close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Hiro_Dataset';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3],[1 2 3], 1}; 


mazeShape.Roy = {'linear'; 'linear';'linear'};
mazeShape.Ted = {'linear'; 'L-shape'; 'U-shape'};
mazeShape.Kevin = {'linear'};


for rr = 1%:length(rats) % loop over all the rats and sessions

    
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

mainDir = fullfile(currDir, 'analysesResults', fileinfo.name);
mkdir(mainDir)



%%% behavior 

if strcmp(rat, 'Kevin') % there are two wake periods, we need to concatenate them
   behavior.time(2,2) = behavior.time(3,2);
   behavior.time(3,1) = behavior.time(4,1);
   behavior.time(3,2) = behavior.time(4,2);
   behavior.time(4,:) = [];
end

behavior.time(2,1) = behavior.time(2,1) + 1e6; % to make sure this time period is only limited to the maze (1e6 = 1 sec)
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
% the information will be loaded to the xyt2 field of fileinfo structure


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

% for the ROY sessions we need to make sure confining the rat on the platform
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

%%% 2D spatial tuning

directory = fullfile(mainDir, 'PlaceFields');
mkdir(directory)


% spatialTunings2D_LR = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, directory);
% spatialTunings2D_RL = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, directory);



%%% 1D spatial tuning (using linearized positions)

% the place fields are calculated separately for the left and right
% direction

[spatialTunings_LR, PF_sorted_LR, PF_peak_LR, activeUnits_sorted_LR, sparsity_LR] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'LR', 2, runSpeedThresh, fileinfo.Fs, directory);
[spatialTunings_RL, PF_sorted_RL, PF_peak_RL, activeUnits_sorted_RL, sparsity_RL] = spatialTuning_1D(spikeStruct, [1 2 3], fileinfo, behavior, speed, 'RL', 2, runSpeedThresh, fileinfo.Fs, directory);


%% Determining PBEs during different behavioral periods


%%% Detemining population burst periods

% For detection of the population bursts, Kamran suggested using all the clusters
% even those inclusing MUA and low amplitude spikes. Therefore, it's better to load the spike timestamps from
% res files. I don't have access to the res files for now, so will
% continure for now with the sorted spikes


% Detection configurations

time_resolution = 0.001; % in second
threshZ = 3; % sdf with 3 std deviation above the mean


%%% PRE PBEs

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(1,2);

directory = fullfile(mainDir, 'PBEs', 'PRE');
mkdir(directory)


% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4


[primaryPBEs_pre, sdat_pre] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, spikeStruct, fileinfo);

savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, fileinfo, directory)

[poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE); % generate poission and within-PBE-time-swap surrogate datasets

plotSurroundMultiunitRate(sdat_pre, secondaryPBEs_pre, fileinfo, directory)




%%% RUN PBEs


fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);

directory = fullfile(mainDir, 'PBEs', 'RUN');
mkdir(directory)

% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4


[primaryPBEs_run, sdat_run] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, spikeStruct, fileinfo);

savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, fileinfo, directory)

[poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN); 

plotSurroundMultiunitRate(sdat_run, secondaryPBEs_run, fileinfo, directory)



%%% POST PBEs


fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);

directory = fullfile(mainDir, 'PBEs', 'POST');
mkdir(directory)

% exclude the active wake and REM periods
exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4


[primaryPBEs_post, sdat_post] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, spikeStruct, fileinfo);

savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, fileinfo, directory)

[poissonEventsBinnedfiringPOST, timeSwapEventsBinnedfiringPOST] = genSurrogates(eventsBinnedfiringPOST); 

plotSurroundMultiunitRate(sdat_post, secondaryPBEs_post, fileinfo, directory)



%% Bayesian decoding (validation of encoding model)


%%% linearized positions


% Validate BD
% Decode the postions during animal's running and comapre it with the
% actual position of the animal. The less the eror the more reliable
% template to be used for replay detection during the PBEs

binDur      = 0.5; % duration of active RUN time bins
posBinSize  = 2; % in cm


directory = fullfile(mainDir, 'PlaceFields');


% LR

laps_LR_idx = 2:2:length(laps);
laps_LR_idx = find(ismember(laps(:,3), laps_LR_idx));

laps_LR    = laps(laps_LR_idx, :);

individualAndPopulationSpatialInfo(spikeStruct, laps_LR, spatialTunings_LR, PF_peak_LR, sparsity_LR, behavior, speed, runSpeedThresh, posBinSize, binDur, 'LR', fileinfo, directory);


close all

% RL

laps_RL_idx = 1:2:length(laps);
laps_RL_idx = find(ismember(laps(:,3), laps_RL_idx));

laps_RL    = laps(laps_RL_idx, :);

individualAndPopulationSpatialInfo(spikeStruct, laps_RL, spatialTunings_RL, PF_peak_RL, sparsity_RL, behavior, speed, runSpeedThresh, posBinSize, binDur, 'RL', fileinfo, directory);


%% Distribution of decoded (Bayesian decoding) track locations pooled across the PBEs 

binDur = 0.02; % in sec
posBinSize = 2; % in cm


% For each time bin consdiering the peak position from the most likely
% direction. The other way is to calculate the decoded positions separately
% for each direction: take a look at function PBEsLocationCoding 


% PRE PBEs

peakDecodedPositionPRE  = PBEsLocationCodingV2(eventsBinnedfiringPRE, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);


% RUN PBEs

peakDecodedPositionRUN  = PBEsLocationCodingV2(eventsBinnedfiringRUN, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);


% POST PBEs

peakDecodedPositionPOST = PBEsLocationCodingV2(eventsBinnedfiringPOST, spatialTunings_LR, spatialTunings_RL, [], posBinSize, binDur, directory);



%%% plot

PosBinWidth = 10;

plotDecodedPositionDuringPBEs(peakDecodedPositionPRE, peakDecodedPositionRUN, peakDecodedPositionPOST, PosBinWidth, fileinfo, directory)


%% Replay evaluation using Bayesian decoding


%%% PBEs

binDur = 0.02;

% PRE

directory = fullfile(mainDir, 'Bayesian', 'PRE');
mkdir(directory)


[posteriorProbMatrixPRE, weightedCorrPRE, jumpDistPRE, BDseqscorePRE, sigMatrixPRE] = BDreplayDetect_v2(eventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data');
% BDreplayDetect(eventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data_twoShuffleMethods');

% 
% % time swap
% eventsBinnedfiringPRE_ts = cell(size(eventsBinnedfiringPRE));
% eventsBinnedfiringPRE_ts(:, 2) = genTimeSwap(eventsBinnedfiringPRE(:, 2));
% 

% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'shuffle');
% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'pooledTimeSwap');
[posteriorProbMatrixPRE_p, weightedCorrPRE_p, jumpDistPRE_p, BDseqscorePRE_p, sigMatrixPRE_p]      = BDreplayDetect_v2(timeSwapEventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'withinPBETimeSwap');
[posteriorProbMatrixPRE_ts, weightedCorrPRE_ts, jumpDistPRE_ts, BDseqscorePRE_ts, sigMatrixPRE_ts] = BDreplayDetect_v2(poissonEventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'poisson');



%[BDprctilePRE_ts, PosteriorProbMatrixPRE_ts, begPositionPRE_ts, endPositionPRE_ts, slopePRE_ts, replaydirectionPRE_ts, replayOrderPRE_ts] = 


% RUN

directory = fullfile(mainDir, 'Bayesian', 'RUN');
mkdir(directory)


[posteriorProbMatrixRUN, weightedCorrRUN, jumpDistRUN, BDseqscoreRUN, sigMatrixRUN] = BDreplayDetect_v2(eventsBinnedfiringRUN, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data');
% BDreplayDetect(eventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data_twoShuffleMethods');

% 
% % time swap
% eventsBinnedfiringPRE_ts = cell(size(eventsBinnedfiringPRE));
% eventsBinnedfiringPRE_ts(:, 2) = genTimeSwap(eventsBinnedfiringPRE(:, 2));
% 

% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'shuffle');
% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'pooledTimeSwap');
[posteriorProbMatrixRUN_p, weightedCorrRUN_p, jumpDistRUN_p, BDseqscoreRUN_p, sigMatrixRUN_p]      = BDreplayDetect_v2(timeSwapEventsBinnedfiringRUN, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'withinPBETimeSwap');
[posteriorProbMatrixRUN_ts, weightedCorrRUN_ts, jumpDistRUN_ts, BDseqscoreRUN_ts, sigMatrixRUN_ts] = BDreplayDetect_v2(poissonEventsBinnedfiringRUN, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'poisson');




% 
% BDreplayDetect(eventsBinnedfiringRUN, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'data');
% 
% 
% 
% % time swap
% eventsBinnedfiringRUN_ts = cell(size(eventsBinnedfiringRUN));
% eventsBinnedfiringRUN_ts(:, 2) = genTimeSwap(eventsBinnedfiringRUN(:, 2));
% 
% 
% BDreplayDetect(eventsBinnedfiringRUN_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'shuffle');
% 
% 
% POST

directory = fullfile(mainDir, 'Bayesian', 'POST');
mkdir(directory)


[posteriorProbMatrixPOST, weightedCorrPOST, jumpDistPOST, BDseqscorePOST, sigMatrixPOST] = BDreplayDetect_v2(eventsBinnedfiringPOST, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data');
% BDreplayDetect(eventsBinnedfiringPRE, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'data_twoShuffleMethods');

% 
% % time swap
% eventsBinnedfiringPRE_ts = cell(size(eventsBinnedfiringPRE));
% eventsBinnedfiringPRE_ts(:, 2) = genTimeSwap(eventsBinnedfiringPRE(:, 2));
% 

% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'shuffle');
% BDreplayDetect(eventsBinnedfiringPRE_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'pooledTimeSwap');
[posteriorProbMatrixPOST_p, weightedCorrPOST_p, jumpDistPOST_p, BDseqscorePOST_p, sigMatrixPOST_p]      = BDreplayDetect_v2(timeSwapEventsBinnedfiringPOST, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'withinPBETimeSwap');
[posteriorProbMatrixPOST_ts, weightedCorrPOST_ts, jumpDistPOST_ts, BDseqscorePOST_ts, sigMatrixPOST_ts] = BDreplayDetect_v2(poissonEventsBinnedfiringPOST, spatialTunings_RL, spatialTunings_LR, [], binDur, directory, 'poisson');



% BDreplayDetect(eventsBinnedfiringPOST, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'data');
% 
% 
% 
% % time swap
% eventsBinnedfiringPOST_ts = cell(size(eventsBinnedfiringPOST));
% eventsBinnedfiringPOST_ts(:, 2) = genTimeSwap(eventsBinnedfiringPOST(:, 2));
% 
% 
% BDreplayDetect(eventsBinnedfiringPOST_ts, spatialTunings_LR, spatialTunings_RL, [], binDur, directory, 'shuffle');


%% plot the BD results for example PBEs


close all
binDur = 0.02;


%%% PRE

% data

directory = fullfile(mainDir, 'Bayesian', 'PRE', 'data');

load(fullfile(directory, 'BDresults_wc_jd.mat'))

BDprctile = BDseqscore.prctilescore;

% BDprctile = BDprctile(:, 1, :);

% BDprctile = squeeze(BDprctile);
hiBDprctile = max(BDprctile(:, 1:2), [], 2);

[~, idx] = sort(hiBDprctile, 'descend');


PBEidx2plot = idx([1:50 60:10:length(idx)]);

% tempPlotBD(eventsBinnedfiringPRE, PosteriorProbMatrix, BDprctile, begPosition, endPosition, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)

tempPlotBD(eventsBinnedfiringPRE, posteriorProbMatrixPRE, BDprctile, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)
% tempPlotBD([], PosteriorProbMatrix, BDprctile, begPosition, endPosition, [], [], PBEidx2plot, binDur, directory)


%%% RUN

% data

directory = fullfile(mainDir, 'Bayesian', 'RUN', 'data');

load(fullfile(directory, 'BDresults_wc_jd.mat'))

BDprctile = BDseqscore.prctilescore;

% BDprctile = BDprctile(:, 1, :);
% 
% BDprctile = squeeze(BDprctile);
hiBDprctile = max(BDprctile(:, 1:2), [], 2);

[~, idx] = sort(hiBDprctile, 'descend');


PBEidx2plot = idx([1:50 60:10:length(idx)]);

% tempPlotBD(eventsBinnedfiringRUN, PosteriorProbMatrix, BDprctile, begPosition, endPosition, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)

tempPlotBD(eventsBinnedfiringRUN, posteriorProbMatrixRUN, BDprctile, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)
% tempPlotBD([], PosteriorProbMatrix, BDprctile, begPosition, endPosition, [], [], PBEidx2plot, binDur, directory)



%%% POST

% data

directory = fullfile(mainDir, 'Bayesian', 'POST', 'data');

load(fullfile(directory, 'BDresults_wc_jd.mat'))

BDprctile = BDseqscore.prctilescore;


% BDprctile = BDprctile(:, 1, :);
% 
% BDprctile = squeeze(BDprctile);
hiBDprctile = max(BDprctile(:, 1:2), [], 2);

[~, idx] = sort(hiBDprctile, 'descend');


PBEidx2plot = idx([1:50 60:10:length(idx)]);

% tempPlotBD(eventsBinnedfiringPOST, PosteriorProbMatrix, BDprctile, begPosition, endPosition, activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)

tempPlotBD(eventsBinnedfiringPOST, posteriorProbMatrixPOST, BDprctile, [], [], activeUnits_sorted_LR, activeUnits_sorted_RL, PBEidx2plot, binDur, directory)

% tempPlotBD([], PosteriorProbMatrix, BDprctile, begPosition, endPosition, [], [], PBEidx2plot, binDur, directory)


%%% plot the (cumulative) distribution of replay score

plotcdfsforBD
% plotcdfsforBD2

%% Example Hidden markov Models trained separately on PRE, RUN, and POST periods


% PRE
numofStates = 40; % dimension of the model should be optimized 

noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);

[transmat0, lambda0, prior0] = initialModel(numofStates, noActiveUnits);
[~, prior, transmatPRE, lambdaPRE] = phmm_em(eventsBinnedfiringPRE(:,2), prior0, transmat0, lambda0, 'max_iter', 200);
[transmatPRE, lambdaPRE] = sortStates(transmatPRE, lambdaPRE, prior);


% RUN
numofStates = 40;
% noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);

[transmat0, lambda0, prior0] = initialModel(numofStates, noActiveUnits);
[~, prior, transmatRUN, lambdaRUN] = phmm_em(eventsBinnedfiringRUN(:,2), prior0, transmat0, lambda0, 'max_iter', 200);
[transmatRUN, lambdaRUN] = sortStates(transmatRUN, lambdaRUN, prior);


% POST
numofStates = 40;
% noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);

[transmat0, lambda0, prior0] = initialModel(numofStates, noActiveUnits);
[~, prior, transmatPOST, lambdaPOST] = phmm_em(eventsBinnedfiringPOST(:,2), prior0, transmat0, lambda0, 'max_iter', 200);
[transmatPOST, lambdaPOST] = sortStates(transmatPOST, lambdaPOST, prior);




%%% plot the example transition and observation matrices from each period

% set the color limits

tmat_temp = [transmatPRE(:); transmatRUN(:); transmatPOST(:)];
tmat_clim = [min(tmat_temp) max(tmat_temp)];


omat_temp = [lambdaPRE(:); lambdaRUN(:); lambdaPOST(:)];
omat_clim = [min(omat_temp) max(omat_temp)];


figure;
set(gcf, 'units', 'centimeters', 'position', [0 0 9 12])


% PRE
subplot(3,2,1)
imagesc(transmatPRE, tmat_clim)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'PRE';'';'State i'}, 'fontsize', 8)
title('Transition matrix', 'fontsize', 8)

subplot(3,2,2)
imagesc(lambdaPRE, omat_clim)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)
title('Observation matrix', 'fontsize', 8)


% RUN
subplot(3,2,3)
imagesc(transmatRUN, tmat_clim)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'RUN';'';'State i'}, 'fontsize', 8)


subplot(3,2,4)
imagesc(lambdaRUN, omat_clim)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)


% POST
subplot(3,2,5)
imagesc(transmatPOST, tmat_clim)
% set(gca, 'YDir', 'normal')
xlabel('State j', 'fontsize', 8)
ylabel({'POST';'';'State i'}, 'fontsize', 8)


subplot(3,2,6)
imagesc(lambdaPOST, omat_clim)
set(gca, 'YDir', 'normal')
xlabel('State', 'fontsize', 8)
ylabel('Unit', 'fontsize', 8)

colormap('jet')

directory = fullfile(mainDir, 'HMM', 'exampleHMMs');
mkdir(directory)

savepdf(gcf, fullfile(directory, [fileinfo.name '-exampleHMMs']), '-dpng')
savepdf(gcf, fullfile(directory, [fileinfo.name '-exampleHMMs']), '-dpdf')


%% Test the sparisity of the models using Gini coefficients

% 250 HMMs trained on actual data with different initializtions compared with 250 models trained on shuffle datasets (temporal shuffle, time swap, and poisson
% simulation)


%% Cross-validation (training and testing models on same periods) to characterize consistency of sequential structures within each period


numofFolds = 5; %% 5-fold cross-validation
numofStates = 40;
noShuffle = 500;

directory = fullfile(mainDir, 'HMM', 'cross-validation');
mkdir(directory)



%%% PRE cross-valodation

pbePeriod = 'PRE';
subDirectory = fullfile(directory, pbePeriod);
mkdir(subDirectory)


rndgen                             = randperm(size(eventsBinnedfiringPRE, 1)); % randomizing the temporal order of the PBEs 
eventsBinnedfiringPRE_rnd          = eventsBinnedfiringPRE(rndgen, :); 

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPRE_rnd, numofStates, numofFolds); 

[HMMprctilePRE, statesProbDistsPRE, binResolvedScoresPRE] = HMMcongruence_crossvalid(eventsBinnedfiringPRE_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', subDirectory);

HMMprctilePRE        = HMMprctilePRE(rndgen);
statesProbDistsPRE   = statesProbDistsPRE(rndgen);
binResolvedScoresPRE = binResolvedScoresPRE(rndgen);

HMMprctilePRE_pts    = HMMcongruence_crossvalid_pts(eventsBinnedfiringPRE_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', subDirectory); % pts: pooled time swap

HMMprctilePRE_pts    = HMMprctilePRE_pts(rndgen);



%%% RUN cross-validation

pbePeriod = 'RUN'; 
subDirectory = fullfile(directory, pbePeriod);
mkdir(subDirectory)

rndgen                             = randperm(size(eventsBinnedfiringRUN, 1));
eventsBinnedfiringRUN_rnd          = eventsBinnedfiringRUN(rndgen, :);


[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN_rnd, numofStates, numofFolds); 

[HMMprctileRUN, statesProbDistsRUN, binResolvedScoresRUN] = HMMcongruence_crossvalid(eventsBinnedfiringRUN_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', subDirectory);

HMMprctileRUN        = HMMprctileRUN(rndgen);
statesProbDistsRUN   = statesProbDistsRUN(rndgen);
binResolvedScoresRUN = binResolvedScoresRUN(rndgen);

HMMprctileRUN_pts    = HMMcongruence_crossvalid_pts(eventsBinnedfiringRUN_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', subDirectory); 

HMMprctileRUN_pts    = HMMprctileRUN_pts(rndgen);



%%% POST cross-validation

pbePeriod = 'POST';
subDirectory = fullfile(directory, pbePeriod);
mkdir(subDirectory)

rndgen                             = randperm(size(eventsBinnedfiringPOST, 1));
eventsBinnedfiringPOST_rnd         = eventsBinnedfiringPOST(rndgen, :);


[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPOST_rnd, numofStates, numofFolds); 

[HMMprctilePOST, statesProbDistsPOST, binResolvedScoresPOST] = HMMcongruence_crossvalid(eventsBinnedfiringPOST_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', subDirectory);

HMMprctilePOST        = HMMprctilePOST(rndgen);
statesProbDistsPOST   = statesProbDistsPOST(rndgen);
binResolvedScoresPOST = binResolvedScoresPOST(rndgen);

HMMprctilePOST_pts                 = HMMcongruence_crossvalid_pts(eventsBinnedfiringPOST_rnd, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', subDirectory);

HMMprctilePOST_pts    = HMMprctilePOST_pts(rndgen);





%%% plot the distributions (CDFs) of HMM cross-valiation scores

figure;
set(gcf, 'Units', 'centimeters', 'position', [0 0 17.6 9])


% PRE
subplot(1,3,1)
congruencePlot2(HMMprctilePRE, HMMprctilePRE_pts, 1, sprintf('%s\n  (n=%d)', 'PRE', length(HMMprctilePRE)))

% RUN
subplot(1,3,2)
congruencePlot2(HMMprctileRUN, HMMprctileRUN_pts, 0, sprintf('%s\n  (n=%d)', 'RUN', length(HMMprctileRUN)))

% POST
subplot(1,3,3)
congruencePlot2(HMMprctilePOST, HMMprctilePOST_pts, 0, sprintf('%s\n  (n=%d)', 'POST', length(HMMprctilePOST)))


savepdf(gcf, fullfile(directory, 'cdfs'), '-dpng')
savepdf(gcf, fullfile(directory, 'cdfs'), '-dpdf')


%% plot the HMM cross-validation results for example PBEs


% For each significantly congruent PBE the virtual trajectory of the
% animal over the track location is calculated using the lsPFs (here it would be better to use the 2D lsPFs)
% (or otherwise, I can use the bayesian decoding)

% latent state place fields (lsPFs)

% 5-fold cross-validation: 

% (1) calculate lsPFs using the training set,
% (2) decode the state in each time bin for the PBEs in the test set,
% (3) calculate the weighted sum of the lsPFs with the weights as the
%     probability distribution over the states,
% (4) calculate the center of mass of the final distribution over the
%     position as the most likely decoded position in the corresponding time
%     bin
% (5) finally to calculate position decoding error by finding the distance
%     between the decoded and actual positions.


mainDir = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
mkdir(mainDir)


binDur = 0.1; % in sec for active wake periods
numofStates = 40;
posBinSize = 0.5; % in cm


% active periods during run

activePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < behavior.time(2,2), :); 
% ActivePeriods = ActivePeriods(50:end, :);


% Binnig the active run period

[runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(spikeStruct, activePeriods, behavior, binDur, posBinSize, fileinfo); % the positionIdx is 2D here (binRunData_1D does the linear)

activeUnits = 1:size(eventsBinnedfiringPRE{1,2}, 1);




%% PRE


pbePeriod = 'PRE';
% numofStates = 40;
% Train an HMM on PBEs
[transmatPRE, lambdaPRE] = trainHMM(eventsBinnedfiringPRE(:,2), activeUnits, numofStates);

% normalise_lsPFs = 1;

lsPFsPRE = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmatPRE, lambdaPRE, posbinIdx, numofStates, xposcenters, yposcenters);


lspfPeakPositionsPRE = zeros(numofStates, 2);
for ii = 1: numofStates
   
    currlspf = lsPFsPRE(:,:,ii);
    [r, c] = find(currlspf == max(currlspf(:)));
    
    lspfPeakPositionsPRE(ii, :) = [c r]*posBinSize;
   
end


%% RUN

pbePeriod = 'RUN';
% numofStates = 30;
% Train an HMM on PBEs
[transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);

% normalise_lsPFs = 1;
lsPFsRUN = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmatRUN, lambdaRUN, posbinIdx, numofStates, xposcenters, yposcenters);


lspfPeakPositionsRUN = zeros(numofStates, 2);
for ii = 1: numofStates
   
    currlspf = lsPFsRUN(:,:,ii);
    [r, c] = find(currlspf == max(currlspf(:)));
    
    lspfPeakPositionsRUN(ii, :) = [c r]*posBinSize;
   
end


%% POST

pbePeriod = 'POST';
% numofStates = 50;
% Train an HMM on PBEs
[transmatPOST, lambdaPOST] = trainHMM(eventsBinnedfiringPOST(:,2), activeUnits, numofStates);

% normalise_lsPFs = 1;
lsPFsPOST = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmatPOST, lambdaPOST, posbinIdx, numofStates, xposcenters, yposcenters);


lspfPeakPositionsPOST = zeros(numofStates, 2);
for ii = 1: numofStates
   
    currlspf = lsPFsPOST(:,:,ii);
    [r, c] = find(currlspf == max(currlspf(:)));
    
    lspfPeakPositionsPOST(ii, :) = [c r]*posBinSize;
   
end





%% plot example PBEs from cross-validation

binDur = 0.02;

%PRE
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection3/cross-validation/PRE'];
mkdir(FileBase)

[~, tempSortIdx] = sort(HMMprctilePRE, 'descend');
PBEidx2plot = tempSortIdx([1:50 60:20:length(HMMprctilePRE)]);

plotExamplePBEs2(eventsBinnedfiringPRE, HMMprctilePRE,  statesProbDistsPRE, lspfPeakPositionsPRE, binResolvedScoresPRE, PBEidx2plot,  binDur, fileinfo, behavior, FileBase)


%RUN
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection3/cross-validation/RUN'];
mkdir(FileBase)

[~, tempSortIdx] = sort(HMMprctileRUN, 'descend');
PBEidx2plot = tempSortIdx([1:50 60:20:length(HMMprctileRUN)]);

plotExamplePBEs2(eventsBinnedfiringRUN, HMMprctileRUN,  statesProbDistsRUN, lspfPeakPositionsRUN, binResolvedScoresRUN, PBEidx2plot,  binDur, fileinfo, behavior, FileBase)


%POST
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection3/cross-validation/POST'];
mkdir(FileBase)

[~, tempSortIdx] = sort(HMMprctilePOST, 'descend');
PBEidx2plot = tempSortIdx([1:50 60:20:length(HMMprctilePOST)]);

plotExamplePBEs2(eventsBinnedfiringPOST, HMMprctilePOST,  statesProbDistsPOST, lspfPeakPositionsPOST, binResolvedScoresPOST, PBEidx2plot,  binDur, fileinfo, behavior, FileBase)



% % % % %% Correspondence between the states from models trained on different periods(PRE, RUN, POST) 
% % % % 
% % % % % The idea is to train models separately on different periods,
% % % % % then decode the probability over states in the test period(s)
% % % % 
% % % % nPBEs_PRE = size(eventsBinnedfiringPRE, 1);
% % % % nPBEs_RUN = size(eventsBinnedfiringRUN, 1);
% % % % nPBEs_POST = size(eventsBinnedfiringPOST, 1);
% % % % 
% % % % 
% % % % PBES_PRE_RUN = [eventsBinnedfiringPRE; eventsBinnedfiringRUN];
% % % % PBES_PRE_POST = [eventsBinnedfiringPRE; eventsBinnedfiringPOST];
% % % % PBES_RUN_POST = [eventsBinnedfiringRUN; eventsBinnedfiringPOST];
% % % % 
% % % % 
% % % % % generate pooled time swap PBEs and the concatenated data
% % % % 
% % % % eventsBinnedfiringPRE_ts = genTimeSwap(eventsBinnedfiringPRE(:,2)); % we don't need the one milisecond binned data here
% % % % eventsBinnedfiringRUN_ts = genTimeSwap(eventsBinnedfiringRUN(:,2));
% % % % eventsBinnedfiringPOST_ts = genTimeSwap(eventsBinnedfiringPOST(:,2));
% % % % 
% % % % PBES_PRE_RUN_ts = [eventsBinnedfiringPRE_ts; eventsBinnedfiringRUN_ts];
% % % % PBES_PRE_POST_ts = [eventsBinnedfiringPRE_ts; eventsBinnedfiringPOST_ts];
% % % % PBES_RUN_POST_ts = [eventsBinnedfiringRUN_ts; eventsBinnedfiringPOST_ts];
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %% find the correspondences between pairs of states 
% % % % 
% % % % % (1) inner products of the population vectors (from observation matrix) 
% % % % 
% % % % 
% % % % % PRE and RUN
% % % % [RUN_PRE_PopVecSim, RUN_PRE_lambdaSortIdx] = PopVecSimilarity(lambdaRUN, lambdaPRE);
% % % %     
% % % % % RUN and POST
% % % % [POST_RUN_PopVecSim, POST_RUN_lambdaSortIdx] = PopVecSimilarity(lambdaPOST, lambdaRUN);
% % % % 
% % % % % PRE and POST
% % % % [POST_PRE_PopVecSim, POST_PRE_lambdaSortIdx] = PopVecSimilarity(lambdaPOST, lambdaPRE);
% % % % 
% % % % %%
% % % % % (2_A) By decoding the states within a test period (concatenating the
% % % % % periods the models were trained on) and calculate for each pair of
% % % % % states their joint decoding probability within each time bin and
% % % % % averaging over all of the bins
% % % % 
% % % % 
% % % % uniformTransMat = 0;
% % % % 
% % % % % if 0 -> The states correspondences is based on observation probabilties
% % % % 
% % % % % if 1 -> the correspondences are based on the gamma distributions, so the temporal information are being involved. 
% % % % % But the null distribution is calculated by shuffling after calculation of
% % % % % the gamma distributions. The focus isn't exclusively on the sequential contents.  
% % % % 
% % % % % if 2 -> controling for coactivities to what extent the sequence
% % % % % correspond to each other just based on the sequential contents. Pooled
% % % % % time swap is used for generating null distributions. This method does not
% % % % % control for the self-transitions.
% % % % 
% % % % % if 3 -> correspondences based on sequential contents. in this case,
% % % % % transition matrix shuffle (self-transitions spared) is used for the null
% % % % % distribution.
% % % % 
% % % % 
% % % % seqContentCorrespondence = 0; 
% % % % 
% % % % % PRE and RUN
% % % % [RUN_PRE, sortRUN_PRE, gammagivenRUN1, gammagivenPRE1, obsProbgivenRUN1, obsProbgivenPRE1] = stateCorrespondence(PBES_PRE_RUN(:, 2), nPBEs_PRE,  lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % % RUN and POST
% % % % [POST_RUN, sortPOST_RUN, gammagivenPOST2, gammagivenRUN2, obsProbgivenPOST2, obsProbgivenRUN2] = stateCorrespondence(PBES_RUN_POST(:, 2), nPBEs_RUN, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % % PRE and POST
% % % % [POST_PRE, sortPOST_PRE, gammagivenPOST3, gammagivenPRE3, obsProbgivenPOST3, obsProbgivenPRE3] = stateCorrespondence(PBES_PRE_POST(:, 2), nPBEs_PRE, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % 
% % % % %%
% % % % uniformTransMat = 0;
% % % % seqContentCorrespondence = 0; % In this section we are not focusing on the sequential content
% % % % 
% % % % % PRE and RUN
% % % % RUN_PRE_ts = stateCorrespondence(PBES_PRE_RUN_ts, nPBEs_PRE, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % % RUN and POST
% % % % POST_RUN_ts = stateCorrespondence(PBES_RUN_POST_ts, nPBEs_RUN, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % % PRE and POST
% % % % POST_PRE_ts = stateCorrespondence(PBES_PRE_POST_ts, nPBEs_PRE, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %%
% % % % % plot the correspondence matrices
% % % % 
% % % % figure;
% % % % set(gcf, 'position', [-1405 166 665 979])
% % % % %
% % % % 
% % % % currTitle = {'Population vectors inner products';''};
% % % % 
% % % % allValues = [RUN_PRE_PopVecSim(:); POST_RUN_PopVecSim(:); POST_PRE_PopVecSim(:)];
% % % % ccodeRange = [min(allValues) max(allValues)];
% % % % 
% % % % h1 = subplot(3,3,1);
% % % % plotPanel(90-RUN_PRE_PopVecSim, sortRUN_PRE, 'RUN', 'PRE', [0 90], currTitle)
% % % % colormap(h1, 'copper')
% % % % 
% % % % h2 = subplot(3,3,2);
% % % % plotPanel(90-POST_RUN_PopVecSim, sortPOST_RUN, 'POST', 'RUN', [0 90], '')
% % % % colormap(h2, 'copper')
% % % % 
% % % % h3 = subplot(3,3,3);
% % % % plotPanel(90-POST_PRE_PopVecSim, sortPOST_PRE, 'POST', 'PRE', [0 90], '')
% % % % colormap(h3, 'copper')
% % % % 
% % % % %
% % % % 
% % % % currTitle = {'Correspondence between states';''};
% % % % 
% % % % allValues = [RUN_PRE(:); POST_RUN(:); POST_PRE(:); RUN_PRE_ts(:); POST_RUN_ts(:); POST_PRE_ts(:)];
% % % % ccodeRange = [prctile(allValues, 0.1) prctile(allValues, 99.9)];
% % % % 
% % % % subplot(3,3,4)
% % % % plotPanel(RUN_PRE, sortRUN_PRE, 'RUN', 'PRE', [0 110], currTitle)
% % % % 
% % % % subplot(3,3,5)
% % % % plotPanel(POST_RUN, sortPOST_RUN, 'POST', 'RUN', [0 110], '')
% % % % 
% % % % subplot(3,3,6)
% % % % plotPanel(POST_PRE, sortPOST_PRE, 'POST', 'PRE', [0 110], '')
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % currTitle = {'Correspondence between states time swap';''};
% % % % 
% % % % subplot(3,3,7)
% % % % plotPanel(RUN_PRE_ts, sortRUN_PRE, 'RUN', 'PRE', [0 110], currTitle)
% % % % 
% % % % subplot(3,3,8)
% % % % plotPanel(POST_RUN_ts, sortPOST_RUN, 'POST', 'RUN', [0 110], '')
% % % % 
% % % % subplot(3,3,9)
% % % % plotPanel(POST_PRE_ts, sortPOST_PRE, 'POST', 'PRE', [0 110], '')
% % % % 
% % % % colormap('copper')
% % % % 
% % % % % FileBase = [currDir '/' fileinfo.name  '/HMM/StateCorrespondence/'];
% % % % % mkdir(FileBase)
% % % % 
% % % % % savepdf(gcf, fullfile(FileBase, [fileinfo.name '_stateCorrespondence_seqContent_uniformTransmat'])); % _seqContent_uniformTransmat
% % % % 
% % % % 
% % % % 
% % % % %% Visualize the correspondences and binned fRates together
% % % % 
% % % % % for each bin plot the binned frates in the first row and gamma in the
% % % % % second row
% % % % noPBE2plot = 100;
% % % % visCorresp(PBES_RUN_POST, POST_RUN, obsProbgivenPOST2, obsProbgivenRUN2, noPBE2plot) 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% 
% 
% %% HMM sequence detection
% 
% % Measuring the congruence of PBEs from a period (e.g., POST) with the model trained on a different period (e.g., RUN)
% 
% 
% noShuffle = 1000;
% % 
% % %%% train PRE___ test RUN
% % 
% % curr_transmat = transmatPRE;
% % curr_lambda = lambdaPRE;
% % 
% % 
% % 
% % trainPeriod = 'PRE';
% % testPeriod = 'RUN';
% % 
% % HMMprctile_PRE_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% % % 
% % % % _____________surrogates
% % % HMMprctile_PRE_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% % % 
% 
% 
% %%% train RUN___ test POST
% 
% curr_transmat = transmatRUN;
% curr_lambda = lambdaRUN;
% 
% 
% 
% trainPeriod = 'RUN';
% testPeriod = 'POST'; 
% 
% [HMMprctile_RUN_POST, gamma_RUN_POST] = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_RUN_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% % 
% % %%% train PRE___ test POST
% % 
% % curr_transmat = transmatPRE;
% % curr_lambda = lambdaPRE;
% % 
% % 
% % 
% % trainPeriod = 'PRE';
% % testPeriod = 'POST'; 
% % 
% % HMMprctile_PRE_POST = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% % % 
% % % % _____________surrogates
% % % HMMprctile_PRE_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% % % 
% 
% %%% train RUN___ test PRE
% 
% curr_transmat = transmatRUN;
% curr_lambda = lambdaRUN;
% 
% 
% 
% trainPeriod = 'RUN';
% testPeriod = 'PRE'; 
% 
% [HMMprctile_RUN_PRE, gamma_RUN_PRE] = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_RUN_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% % 
% % 
% % %%% train POST___ test PRE
% % 
% % curr_transmat = transmatPOST;
% % curr_lambda = lambdaPOST;
% % 
% % 
% % 
% % trainPeriod = 'POST';
% % testPeriod = 'PRE'; 
% % 
% % HMMprctile_POST_PRE = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% % 
% % % _____________surrogates
% % HMMprctile_POST_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% % 
% % %%% train POST___ test RUN
% % 
% % curr_transmat = transmatPOST;
% % curr_lambda = lambdaPOST;
% % 
% % 
% % 
% % trainPeriod = 'POST';
% % testPeriod = 'RUN'; 
% % 
% % HMMprctile_POST_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% % 
% % % _____________surrogates
% % HMMprctile_POST_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% % 
% 
% 
% % plot example PBEs from cross-validation
% 
% 
% %%
% % binDur = 0.02;
% % 
% % %POST
% % 
% % [~, tempSortIdx] = sort(HMMprctile_RUN_POST, 'descend');
% % PBEidx2plot = tempSortIdx([1:200]);
% % 
% % plotExamplePBEs2(eventsBinnedfiringPOST, HMMprctile_RUN_POST,  gamma_RUN_POST, [], [], PBEidx2plot,  binDur, fileinfo, behavior, directory)
% % 
% 
% 
% %%
% 
% 
% figure; 
% 
% set(gcf, 'position', [0 0 482 198])
% 
% 
% subplot(1,2,1)
% congruencePlot2(HMMprctile_RUN_PRE, HMMprctile_RUN_PRE_pts, 1, sprintf('test on %s', 'PRE'))
% 
% 
% subplot(1,2,2)
% congruencePlot2(HMMprctile_RUN_POST, HMMprctile_RUN_POST_pts, 0, sprintf('test on %s', 'POST'))
% 
% 
% 
% % 
% 
% folderBase = [FileBase  '/HMM/sequenceDetection/congruence/'];
% mkdir(folderBase)
% 
% 
% print(gcf, [folderBase fileinfo.name '_congruence'], '-dpng')
% savepdf(gcf, [folderBase fileinfo.name '_congruence'], '-dpdf')
% 
% 
% 
% 
% % figure;
% % set(gcf, 'position', [129 64 2408 1232])
% 
% % 
% % % Cross-validations 
% % 
% % subplot(2,6,1)
% % congruencePlot2(HMMprctilePRE, HMMprctilePRE_pts, 1, sprintf('%s Cross-validation (N = %d)', 'PRE', length(HMMprctilePRE)))
% % 
% % subplot(2,6,2)
% % congruencePlot2(HMMprctileRUN, HMMprctileRUN_pts, 1, sprintf('%s Cross-validation (N = %d)', 'RUN', length(HMMprctileRUN)))
% % 
% % subplot(2,6,3)
% % congruencePlot2(HMMprctilePOST, HMMprctilePOST_pts, 1, sprintf('%s Cross-validation (N = %d)', 'POST', length(HMMprctileRUN)))
% 
% % 
% % 
% % % Congruence of PBEs with other periods' models
% % 
% % subplot(2,6,7)
% % congruencePlot2(HMMprctile_PRE_RUN, HMMprctile_PRE_RUN_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'RUN'))
% % 
% % subplot(2,6,8)
% % congruencePlot2(HMMprctile_RUN_POST, HMMprctile_RUN_POST_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'POST'))
% % 
% % subplot(2,6,9)
% % congruencePlot2(HMMprctile_PRE_POST, HMMprctile_PRE_POST_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'POST'))
% % 
% % subplot(2,6,10)
% % congruencePlot2(HMMprctile_RUN_PRE, HMMprctile_RUN_PRE_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'PRE'))
% % 
% % subplot(2,6,11)
% % congruencePlot2(HMMprctile_POST_PRE, HMMprctile_POST_PRE_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'PRE'))
% % 
% % subplot(2,6,12)
% % congruencePlot2(HMMprctile_POST_RUN, HMMprctile_POST_RUN_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'RUN'))
% % 
% 
% % 
% % % % % currDir = pwd;
% % % % % FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/congruence/'];
% % % % % mkdir(FileBase)
% % % % % 
% % % % % 
% % % % % print(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
% % % % % savepdf(gcf, [FileBase fileinfo.name '_congruence'])
% % % % % 
% % % % % 
% % % % % save([FileBase fileinfo.name '_congruence.mat'], ...
% % % % %     'HMMprctilePRE', 'HMMprctilePRE_pts', ...
% % % % %     'HMMprctileRUN', 'HMMprctileRUN_pts', ...
% % % % %     'HMMprctilePOST', 'HMMprctilePOST_pts', ...
% % % % %     'HMMprctile_PRE_RUN', 'HMMprctile_PRE_RUN_pts', ...
% % % % %     'HMMprctile_RUN_POST', 'HMMprctile_RUN_POST_pts', ...
% % % % %     'HMMprctile_PRE_POST', 'HMMprctile_PRE_POST_pts', ...
% % % % %     'HMMprctile_RUN_PRE', 'HMMprctile_RUN_PRE_pts', ...
% % % % %     'HMMprctile_POST_PRE', 'HMMprctile_POST_PRE_pts', ...
% % % % %     'HMMprctile_POST_RUN', 'HMMprctile_POST_RUN_pts')
% % % % %  
% % % % 
% % % % 
% % % % %% 2D histograms of sequence scores 
% % % % 
% % % % % example: likelihood(POST|PRE) vs likelihood(POST|RUN)
% % % % 
% % % % % The extent to which RUN signifant PBEs are congruent with model traind of PRE
% % % % 
% % % % figure;
% % % % set(gcf, 'position', [796 125 1658 1184])
% % % % 
% % % % 
% % % % % plot 2D histogram
% % % % 
% % % % data = [HMMprctile_PRE_RUN HMMprctileRUN];
% % % % h = hist3(data, [10 10])/length(data);
% % % % 
% % % % 
% % % % subplot(2,2,1)
% % % % 
% % % % imagesc(h); set(gca, 'YDir', 'normal'); 
% % % % colormap(jet)
% % % % % colormap(flipud(colormap('hot')))
% % % % set(gca, 'fontsize', 16, 'linewidth', 3)
% % % % 
% % % % set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% % % % 
% % % % xlabel({'Conguence of RUN PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
% % % % ylabel({'Conguence of RUN PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)
% % % % 
% % % % title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)
% % % % 
% % % % 
% % % % axis square
% % % % colorbar
% % % % 
% % % % 
% % % % 
% % % % % plot stacked histogram
% % % % 
% % % % runIdx = HMMprctileRUN > 0.95;
% % % % prerunIdx = HMMprctile_PRE_RUN > 0.95;
% % % % 
% % % % noncong = [length(find(~runIdx & ~prerunIdx));
% % % %            length(find(~runIdx & prerunIdx))];
% % % %        
% % % % cong = [length(find(runIdx & ~prerunIdx));
% % % %         length(find(runIdx & prerunIdx))];      
% % % % 
% % % %        
% % % % noncong_cumcount = cumsum(noncong);
% % % % noncong_cumcount = noncong_cumcount (end:-1:1);
% % % % 
% % % % cong_cumcount = cumsum(cong);
% % % % cong_cumcount = cong_cumcount(end:-1:1);
% % % % 
% % % % thecolors =  [50 150 255 ;255 255 150]/255;
% % % % 
% % % % 
% % % % 
% % % % subplot(2,2,2)
% % % % 
% % % % 
% % % % hold on
% % % % for ii = 1:2
% % % %     bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
% % % % end
% % % % hold off
% % % % 
% % % % 
% % % % set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
% % % % set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% % % % 
% % % % xlabel('Congruence with RUN model', 'fontsize', 16)
% % % % ylabel('Number of PBEs', 'fontsize', 20)
% % % % 
% % % % % axis square
% % % % 
% % % % legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
% % % % legend boxoff 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % % Congruence of post PBEs with PRE and RUN models: how the two are related 
% % % % 
% % % % % plot 2D histogram
% % % % 
% % % % data = [HMMprctile_PRE_POST HMMprctile_RUN_POST];
% % % % h = hist3(data, [10 10])/length(data);
% % % % 
% % % % subplot(2,2,3)
% % % % 
% % % % imagesc(h); set(gca, 'YDir', 'normal'); 
% % % % colormap(jet)
% % % % % colormap(flipud(colormap('hot')))
% % % % set(gca, 'fontsize', 16, 'linewidth', 3)
% % % % 
% % % % set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% % % % 
% % % % xlabel({'Conguence of POST PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
% % % % ylabel({'Conguence of POST PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)
% % % % 
% % % % title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)
% % % % 
% % % % 
% % % % axis square
% % % % colorbar
% % % % 
% % % % 
% % % % 
% % % % % plot stacked histogram
% % % % 
% % % % runIdx = HMMprctile_RUN_POST > 0.95;
% % % % prerunIdx = HMMprctile_PRE_POST > 0.95;
% % % % 
% % % % noncong = [length(find(~runIdx & ~prerunIdx));
% % % %            length(find(~runIdx & prerunIdx))];
% % % %        
% % % % cong = [length(find(runIdx & ~prerunIdx));
% % % %         length(find(runIdx & prerunIdx))];      
% % % % 
% % % %        
% % % % noncong_cumcount = cumsum(noncong);
% % % % noncong_cumcount = noncong_cumcount (end:-1:1);
% % % % 
% % % % cong_cumcount = cumsum(cong);
% % % % cong_cumcount = cong_cumcount(end:-1:1);
% % % % 
% % % % thecolors =  [50 150 255 ;255 255 150]/255;
% % % % 
% % % % 
% % % % subplot(2,2,4)
% % % % 
% % % % 
% % % % hold on
% % % % for ii = 1:2
% % % %     bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
% % % % end
% % % % hold off
% % % % 
% % % % 
% % % % set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
% % % % set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% % % % 
% % % % xlabel('Congruence of POST with RUN model', 'fontsize', 16)
% % % % ylabel('Number of PBEs', 'fontsize', 20)
% % % % 
% % % % % axis square
% % % % 
% % % % legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
% % % % legend boxoff 
% % % % 
% % % % 
% % % % FileBase = [currDir  '/HMM/sequenceDetection/overlap/'];
% % % % mkdir(FileBase)
% % % % 
% % % % print(gcf, [FileBase  '_' '2Dhistograms'], '-dpng')
% % % % savepdf(gcf, [FileBase  '_' '2Dhistograms'])
% % % % 
% % % % 
% % 
% % 
% % 
% % % 
% % % %% latent state place fields (lsPFs); Decoding of position using lsPFs and decoding error
% % % 
% % % % 5-fold cross-validation: 
% % % 
% % % % (1) calculate lsPFs using the training set,
% % % % (2) decode the state in each time bin for the PBEs in the test set,
% % % % (3) calculate the weighted sum of the lsPFs with the weights as the
% % % %     probability distribution over the states,
% % % % (4) calculate the center of mass of the final distribution over the
% % % %     position as the most likely decoded position in the corresponding time
% % % %     bin
% % % % (5) finally to calculate position decoding error by finding the distance
% % % %     between the decoded and actual positions.
% % % 
% % % 
% % % FileBase = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
% % % mkdir(FileBase)
% % % 
% % % 
% % % binDur = 0.1; % in sec for active wake periods
% % % numofStates = 40;
% % % posBinSize = 0.5; % in cm
% % % 
% % % 
% % % % active periods during run
% % % 
% % % activePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < behavior.time(2,2), :); 
% % % % ActivePeriods = ActivePeriods(50:end, :);
% % % 
% % % 
% % % % Binnig the active run period
% % % 
% % % [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(spikeStruct, activePeriods, behavior, binDur, posBinSize, fileinfo);
% % % 
% % % activeUnits = 1:size(eventsBinnedfiringPRE{1,2}, 1);
% % % 
% % % 
% % % 
% % % %% PRE
% % % 
% % % 
% % % pbePeriod = 'PRE';
% % % 
% % % % Train an HMM on PBEs
% % % [transmatPRE, lambdaPRE] = trainHMM(eventsBinnedfiringPRE(:,2), activeUnits, numofStates);
% % % 
% % % % normalise_lsPFs = 1;
% % % plot_lsPF(runBinnedfiring, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% % % 
% % % savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% % % print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% % % 
% % % [PosDecErrorPRE, PosDecErrorPRE_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates);
% % % 
% % % 
% % % 
% % % %% RUN
% % % 
% % % pbePeriod = 'RUN';
% % % 
% % % % Train an HMM on PBEs
% % % [transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);
% % % 
% % % % normalise_lsPFs = 1;
% % % plot_lsPF(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% % % 
% % % savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% % % print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% % % 
% % % 
% % % [PosDecErrorRUN, PosDecErrorRUN_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates);
% % % 
% % % 
% % % %% POST
% % % 
% % % pbePeriod = 'POST';
% % % 
% % % % Train an HMM on PBEs
% % % [transmatPOST, lambdaPOST] = trainHMM(eventsBinnedfiringPOST(:,2), activeUnits, numofStates);
% % % 
% % % % normalise_lsPFs = 1;
% % % plot_lsPF(runBinnedfiring, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% % % 
% % % 
% % % savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% % % print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% % % 
% % % 
% % % [PosDecErrorPOST, PosDecErrorPOST_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates);
% % % 
% % % 
% % % 
% % % 
% % % %% plot decoding error
% % % 
% % % 
% % % figure;
% % % set(gcf, 'Units', 'pixels', 'position', [814 399 1648 606])
% % % 
% % % subplot(1,3,1)
% % % plotErrorCdf(PosDecErrorPRE, PosDecErrorPRE_shuffledlsPFs, 1)
% % % title('PRE', 'fontsize', 16)
% % % 
% % % 
% % % subplot(1,3,2)
% % % plotErrorCdf(PosDecErrorRUN, PosDecErrorRUN_shuffledlsPFs, 0)
% % % title('RUN', 'fontsize', 16)
% % % 
% % % 
% % % subplot(1,3,3)
% % % plotErrorCdf(PosDecErrorPOST, PosDecErrorPOST_shuffledlsPFs, 0)
% % % title('POST', 'fontsize', 16)
% % % 
% % % 
% % % % suptitle([fileinfo.name])
% % % 
% % % savepdf(gcf, [FileBase fileinfo.name '_decodingError'])
% % % print(gcf, [FileBase fileinfo.name  '_decodingError'], '-dpng')
% % % 
% % % save([FileBase fileinfo.name '_decodingError.mat'], 'PosDecErrorPRE', 'PosDecErrorPRE_shuffledlsPFs', ...
% % %           'PosDecErrorRUN', 'PosDecErrorRUN_shuffledlsPFs', 'PosDecErrorPOST', 'PosDecErrorPOST_shuffledlsPFs')
% % % 
% % % 
% % 
% % 
% % 
% % %%
% % % % %% Calculating the sequence score over time 
% % % % 
% % % % 
% % % % bvrTimeList = behavior.list(:,[1 2])/Fs;
% % % % bvrState = behavior.list(:,3);
% % % % 
% % % % 
% % % % figure;
% % % % set(gca, 'position', [995 783 796 538])
% % % % hold on
% % % % 
% % % % 
% % % % 
% % % % % (1) time courses during RUN (for PRE and RUN(CV))
% % % % 
% % % % 
% % % % tbeginRUN = behavior.time(2,1)/Fs;
% % % % tendRUN = behavior.time(2,2)/Fs;
% % % % 
% % % % PBEtimesRUN = secondaryPBEs_run(:,1)/Fs;
% % % % 
% % % % 
% % % % 
% % % % %QW periods during RUN
% % % % 
% % % % qw = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% % % % qwRUN = qw(qw(:,1) > tbeginRUN & qw(:,2) < tendRUN, :);
% % % % qwRUN = qwRUN - tbeginRUN;
% % % % 
% % % % 
% % % % subplot(2,1,1)
% % % % PlotSeqFreqvsTime(PBEtimesRUN, HMMprctileRUN, HMMprctile_PRE_RUN, 'RUN', 'PRE', 'RUN', qwRUN, 12, tbeginRUN, tendRUN)
% % % % 
% % % % 
% % % % 
% % % % % (2) time courses during POST (for PRE and RUN(CV))
% % % % 
% % % % 
% % % % tbeginPOST = behavior.time(3,1)/Fs;
% % % % tendPOST = behavior.time(3,2)/Fs;
% % % % 
% % % % PBEtimesPOST = secondaryPBEs_post(:,1)/Fs;
% % % % 
% % % % %nrem periods during POST
% % % % 
% % % % nrem = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% % % % nremPOST = nrem(nrem(:,1) > tbeginPOST & nrem(:,2) < tendPOST, :);
% % % % nremPOST = nremPOST - tbeginPOST;
% % % % 
% % % % %
% % % % subplot(2,1,2)
% % % % PlotSeqFreqvsTime(PBEtimesPOST, HMMprctile_RUN_POST, HMMprctile_PRE_POST, 'RUN', 'PRE', 'POST', nremPOST, 12, tbeginPOST, tendPOST)
% % % % 
% % % % 
% % % % FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceScoreOverTime/'];
% % % % mkdir(FileBase)
% % % % 
% % % % 
% % % % 
% % % % savepdf(gcf, [FileBase '/' fileinfo.name '_seqvstime'])
% % % % print(gcf, [FileBase '/' fileinfo.name  '_seqvstime'], '-dpng')
% % % % 
% % % % 
% % % % 
% % % 
% % % 
% % % %%
% % % 
% % % % Train on the first 15 min epoch in the post and test on the following
% % % % epochs
% % % % We are extending this to train on all individual epochs on at a time and test on the
% % % % remaining ones
% % % 
% % % 
% % % % Train and test on non-overlapping POST epochs
% % % 
% % % tbegin = behavior.time(3,1); 
% % % tend = behavior.time(3,2);
% % % 
% % % PBEtimesPOST = secondaryPBEs_post(:,1);
% % % 
% % % noEpochs = 12; % epochs with duration of ~ 15 min
% % % 
% % % epochs = linspace(tbegin, tend, noEpochs+1);
% % % epochDur = epochs(2) - epochs(1);
% % % 
% % % epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;
% % % 
% % % [~, pbeEpochInd] = histc(PBEtimesPOST, epochs);
% % % 
% % % 
% % % numofStates = 40;
% % % PBELikelihood = sequenceScorevstime_CV(eventsBinnedfiringPOST, pbeEpochInd, numofStates, fileinfo);
% % % 
% % % 
% % % figure; 
% % % hold on
% % % 
% % % cmap = colormap('jet');
% % % plotcolors = cmap(floor(linspace(20,64,12)), :);
% % % 
% % % meanofmean = zeros(noEpochs, 1);
% % % 
% % % for ii = 1: noEpochs
% % %     
% % %     currLikelihood = PBELikelihood{ii, 1};
% % %     
% % %     meanLL = nan(noEpochs, 1);
% % %     stdLL = nan(noEpochs, 1);
% % %     
% % %     for jj = 1: noEpochs
% % %         
% % %         meanLL(jj) = mean(currLikelihood(pbeEpochInd == jj));
% % %         stdLL(jj) = std(currLikelihood(pbeEpochInd == jj))/ numel(currLikelihood(pbeEpochInd == jj));
% % %         
% % %     end
% % %     
% % %     meanofmean(ii) = nanmean(currLikelihood);
% % %     
% % %     errorbar(epochCenters/3600, meanLL, stdLL, 'color', plotcolors(ii, :), 'linewidth', 2)
% % %     
% % % end
% % %     
% % %     
% % % 
% % % % Train on PRE and test on POST epochs
% % % 
% % % 
% % % 
% % % % Train on RUN and test on POST epochs
% % % 
% % % 
% % 
% % 

end

end










