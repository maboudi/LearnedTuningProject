clear; clc; close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/data/GrossmarkData_BapunsClusters';
cd(currDir)


%% loading session data

VarList = {'spikes','behavior','position'};

% note that the spikes within the .mat file are from just the pyramidal
% units (only those that are stable), for MUA I need to load the clu and res files again


% Although no difference regarding the behavior and position data

for var = 1 : length(VarList)
    load([currDir '/CRCNSReclustered-' VarList{var} '.mat'])
end
spikes1 = spikes;
behavior1 = behavior;
position1 = Position;

sessionNames = fieldnames(spikes1);
noSessions = numel(sessionNames);


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100
mazeLimits = [-inf inf -inf 13; -inf inf -inf inf; -inf inf -140 inf; -inf inf -inf inf; -20 inf -inf 110];


for currSession = 1


sessionName = sessionNames{currSession};


spikes = eval(['spikes1.' sessionName]);
behavior = eval(['behavior1.' sessionName]);
position = eval(['position1.' sessionName]);

behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity



%% sessioninfo

fileBase = [currDir '/' sessionName '/' sessionName];

fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


% Load the recording configurations from the xml file

Par = LoadXml(fileBase);
fileinfo.lfpSampleRate = Par.lfpSampleRate; 
fileinfo.nCh = Par.nChannels;
fileinfo.Fs = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording




% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;


fileinfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track

fileinfo.xyt(find(ypos > mazeLimits(currSession, 4) | ypos < mazeLimits(currSession, 3) | xpos > mazeLimits(currSession, 2) | xpos < mazeLimits(currSession, 1)), 1:2) = NaN;



%% Behavioral states

% nrem = 1, Drowsy = 3.5; rem = 2, Intermediate = 1.5, wake = 4


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Can we consider them as NREM???



bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ... % including both the active and quiet wake periods
            1 *   ones(length(behavior.Drowsy), 1); ... % considering Drowsy as a part of NREM, for now 
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering Intermediate as a part of NREM, for now 
            2   *   ones(length(behavior.REM), 1)];



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

% mySpikes = spikeBehaviorAnalysis(spikes, laps, ripple.time, speed, qual2consider, fileinfo, Fs);
mySpikes = spikeBehaviorAnalysis2(spikes, speed, qual2consider, fileinfo);



% We need the multiunit (any detected spike) as well for defining the population burst events

numberofShanks = length(dir([fileBase '.res.*']));

MUA.t = [];
for shank = 1:numberofShanks
    
    currSpikes = load([fileBase '.res.' num2str(shank)]);
    MUA.t = [MUA.t; currSpikes]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');
% The originla smapling frequency is 20000

MUA.t = MUA.t./20000;


% %% Ripple detection

 
% % Here we intend to detect only ripples that have on average relatively large
% % amplitudes over the all four shanks resided in CA1. Alternatively, we can
% % include also ripples (with usually smaller amplitude) having a high
% amplitude on just one or few shanks.
% 
% 
% 
% best_channels = BigRippleChannels(fileinfo,1);
% fileinfo.bestch = best_channels;
% 
% rippPowThresh = 3; 
% rippleEvents = RippleDetect(fileinfo, rippPowThresh); % in lfpSampleRate (beg peak end normalized amplitude)
% 
% Filename = [FileBase '.spl.evt'];
% MakeEvtFile(rippleEvents(:, 1:3), Filename, {'beg', 'peak', 'end'}, fileinfo.lfpSampleRate, 1)
% 
% 
% eventVelocity = interp1(speed.t, velocity, (rippleEvents(:,2)/fileinfo.lfpSampleRate *1e6) + periods(1,1));
% rippleEvents2 = rippleEvents(eventVelocity < 5, :); % refining the ripples events; excluding the ones that happen during mobile state
% 
% 

%% PBEs


%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the meanposbinIdx(:,1)

thetaPeriods = load([fileBase '.theta.1']);
thetaPeriods = floor(thetaPeriods./fileinfo.lfpSampleRate);



%% PRE

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(1,2);


exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, wake=4


% add to the exclude set all the theta periods (they may overlap with the existing period but doesn't matter)

% exclude = [exclude; thetaPeriods];
[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude = exclude(sortIdx, :);


fileinfo.lfpSampleRate = 1;
primaryPBEs_pre = PBPeriods(MUA, fileinfo, [], time_resolution, threshZ, exclude);


[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, mySpikes, fileinfo);

savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, 'PRE', fileinfo)

% [poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE, temporalEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE);




%% RUN


fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [1 2]),:); % nrem=1, rem=2, quiet=3, wake=4


% add theta periods to exclude
% exclude = [exclude; thetaPeriods]; 


[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude = exclude(sortIdx, :);


primaryPBEs_run = PBPeriods(MUA, fileinfo, [], time_resolution, 3, exclude);


[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);

savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, 'RUN', fileinfo)

% [poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, temporalEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);




%% POST

fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);


exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4


% add theta periods to exclude
% exclude = [exclude; thetaPeriods];
[~, sortIdx] = sort(exclude(:,1), 'ascend');
exclude = exclude(sortIdx, :);


primaryPBEs_post = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);


[eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, mySpikes, fileinfo);

savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, 'POST', fileinfo)

% [poissonEventsBinnedfiringPOST, timeSwapEventsBinnedfiringPOST, temporalEventsBinnedfiringPOST] = genSurrogates(eventsBinnedfiringPOST);



%% latent state place fields (lsPFs); Decoding of position using lsPFs and decoding error

% 5-fold cross-validation: 

% (1) calculate lsPFs using the training set,

% (2) decode the state in each time bin for the PBEs in the test set,

% (3) calculate the weighted sum of the lsPFs with the weights as the
% probability distribution over the states,

% (4) calculate the center of mass of the final distribution over the
% position as the most likely decoded position in the corresponding time
% bin

% (5) finally to calculate position decoding error by finding the distance
% between the decoded position using lsPFs and actual position of the animal on the track 


currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
mkdir(FileBase)


binDur = 0.1; % in sec 
numofStates = 40;
posBinSize = 0.5; % in cm 


% active periods during run
activePeriods = thetaPeriods(find(thetaPeriods(:,1) > behavior.time(2,1) & thetaPeriods(:,2) < behavior.time(2,2)), :); 


% Binnig the active run period

[runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo);

activeUnits = 1:size(eventsBinnedfiringRUN{1,2}, 1);



%% PRE

pbePeriod = 'PRE';


% Train an HMM on PBEs
[transmatPRE, lambdaPRE] = trainHMM(eventsBinnedfiringPRE(:,2), activeUnits, numofStates);


% normalise_lsPFs = 1;
plot_lsPF(runBinnedfiring, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);

savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')


[PosDecErrorPRE, PosDecErrorPRE_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates);



%% RUN

pbePeriod = 'RUN';

% Train an HMM on PBEs
[transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);


% normalise_lsPFs = 1;
plot_lsPF(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);

savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')


[PosDecErrorRUN, PosDecErrorRUN_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates);



% %% POST
% 
% pbePeriod = 'POST';
% 
% % Train an HMM on PBEs
% [transmatPOST, lambdaPOST] = trainHMM(eventsBinnedfiringPOST(:,2), activeUnits, numofStates);
% 
% 
% % normalise_lsPFs = 1;
% plot_lsPF(runBinnedfiring, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% 
% 
% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% 
% [PosDecErrorPOST, PosDecErrorPOST_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates);
% 
% 
% 
% 
% % 
% % posbinIdx2 = posbinIdx;
% % 
% % 
% % excludeBins = find(posbinIdx(:,1).* posbinIdx(:,2) == 0);
% % posbinIdx2(excludeBins, :) = [];
% % 
% % posbinIdx2 = posbinIdx2 - 1;
% % posbinIdx2 (posbinIdx2 <= 0) = 1;
% % 
% % 
% % probabilityOverPosition(:, :, excludeBins) = [];
% % 
% % for ii = 1:length(posbinIdx)
% %     
% %     probabilityOverPosition(posbinIdx2(ii,2), posbinIdx2(ii,1), ii) = 1;
% %     
% % end
% % 
% 
% 
% % maxGraylevelPerFrame = max(reshape(probabilityOverPosition, numel(probabilityOverPosition(:,:,1)), size(probabilityOverPosition, 3)));
% % normArray = repmat(permute(maxGraylevelPerFrame, [1 3 2]), [size(probabilityOverPosition, 1) size(probabilityOverPosition, 2) 1]);
% % 
% % norm_probabilityOverPosition = probabilityOverPosition ./ normArray;
% % indexedImage = gray2ind(norm_probabilityOverPosition);
% % 
% % 
% % cmap = colormap('jet');
% % movie = immovie(permute(indexedImage, [1 2 4 3]), cmap);
% % implay(movie, 10)
% 
% 
% %% plot decoding error
% 
% 
% figure;
% set(gcf, 'Units', 'pixels', 'position', [814 399 1648 606])
% 
% 
% subplot(1,3,1)
% plotErrorCdf(PosDecErrorPRE, PosDecErrorPRE_shuffledlsPFs, 1)
% title('PRE', 'fontsize', 16)
% 
% 
% subplot(1,3,2)
% plotErrorCdf(PosDecErrorRUN, PosDecErrorRUN_shuffledlsPFs, 0)
% title('RUN', 'fontsize', 16)
% 
% 
% subplot(1,3,3)
% plotErrorCdf(PosDecErrorPOST, PosDecErrorPOST_shuffledlsPFs, 0)
% title('POST', 'fontsize', 16)
% 
% 
% % suptitle([fileinfo.name])
% 
% savepdf(gcf, [FileBase fileinfo.name '_decodingError'])
% print(gcf, [FileBase fileinfo.name  '_decodingError'], '-dpng')
% 
% save([FileBase fileinfo.name '_decodingError.mat'], 'PosDecErrorPRE', 'PosDecErrorPRE_shuffledlsPFs', ...
%           'PosDecErrorRUN', 'PosDecErrorRUN_shuffledlsPFs', 'PosDecErrorPOST', 'PosDecErrorPOST_shuffledlsPFs')
% 
% 
% 
% 
% 
% %% HMM sequence detection
% 
% % In this step we test the extent to which an HMM trained on the
% % candidate events is able to detect sequence events among non-sequence
% % events
% 
% % We do this by k-fold cross-validation: using k-1 folds for training and
% % and the remaining fold for sequence detection. Randomizing the events so the events picked up from random periods to
% % train and test the models
% 
% 
% numofFolds = 5; %% 5-fold cross-validation
% numofStates = 40;
% noShuffle = 500;
% 
% 
% 
% %% PRE
% 
% pbePeriod = 'PRE';
% 
% eventsBinnedfiringPRE = eventsBinnedfiringPRE(randperm(size(eventsBinnedfiringPRE, 1)), :);
% 
% 
% [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPRE, numofStates, numofFolds, fileinfo, pbePeriod); 
% 
% HMMprctilePRE = HMMcongruence_crossvalid(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
% HMMprctilePRE_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);
% 
% 
% %% RUN
% 
% pbePeriod = 'RUN';
% 
% eventsBinnedfiringRUN = eventsBinnedfiringRUN(randperm(size(eventsBinnedfiringRUN, 1)), :);
% 
% 
% [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN, numofStates, numofFolds, fileinfo, pbePeriod); 
% 
% HMMprctileRUN = HMMcongruence_crossvalid(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
% HMMprctileRUN_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);
% 
% 
% %% POST
% 
% pbePeriod = 'POST';
% 
% eventsBinnedfiringPOST = eventsBinnedfiringPOST(randperm(size(eventsBinnedfiringPOST, 1)), :);
% 
% 
% [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPOST, numofStates, numofFolds, fileinfo, pbePeriod); 
% 
% HMMprctilePOST = HMMcongruence_crossvalid(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
% HMMprctilePOST_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);
% 
% 
% 


%% Measuring the congruence of PBEs within a state (e.g., POST) with the model trained on a different period (e.g., RUN)


% PRE
numofStates = 50;
noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);
prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);

[~, prior, transmatPRE, lambdaPRE] = phmm_em(eventsBinnedfiringPRE(:,2), prior0, transmat0, lambda0, 'max_iter', 200);

[~, startIdx] = max(prior);
sortIdx = sortStates(transmatPRE, startIdx);
transmatPRE = transmatPRE(sortIdx, sortIdx);
lambdaPRE = lambdaPRE(:, sortIdx);




% RUN
numofStates = 30;
noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);
prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);

[~, prior, transmatRUN, lambdaRUN] = phmm_em(eventsBinnedfiringRUN(:,2), prior0, transmat0, lambda0, 'max_iter', 200);

[~, startIdx] = max(prior);
sortIdx = sortStates(transmatRUN, startIdx);
transmatRUN = transmatRUN(sortIdx, sortIdx);
lambdaRUN = lambdaRUN(:, sortIdx);



% POST
numofStates = 50;
noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);
prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);

[~, prior, transmatPOST, lambdaPOST] = phmm_em(eventsBinnedfiringPOST(:,2), prior0, transmat0, lambda0, 'max_iter', 200);

[~, startIdx] = max(prior);
sortIdx = sortStates(transmatPOST, startIdx);
transmatPOST = transmatPOST(sortIdx, sortIdx);
lambdaPOST = lambdaPOST(:, sortIdx);




%%
% noShuffle = 500;
% 
% 
% %%% train PRE___ test RUN
% 
% curr_transmat = transmatPRE;
% curr_lambda = lambdaPRE;
% 
% trainPeriod = 'PRE';
% testPeriod = 'RUN';
% 
% HMMprctile_PRE_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_PRE_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% %%% train RUN___ test PRE
% 
% curr_transmat = transmatRUN;
% curr_lambda = lambdaRUN;
% 
% trainPeriod = 'RUN';
% testPeriod = 'PRE'; 
% 
% HMMprctile_RUN_PRE = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_RUN_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% %%% train RUN___ test POST
% 
% curr_transmat = transmatRUN;
% curr_lambda = lambdaRUN;
% 
% trainPeriod = 'RUN';
% testPeriod = 'POST'; 
% 
% HMMprctile_RUN_POST = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_RUN_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% %%% train PRE___ test POST
% 
% curr_transmat = transmatPRE;
% curr_lambda = lambdaPRE;
% 
% 
% trainPeriod = 'PRE';
% testPeriod = 'POST'; 
% 
% HMMprctile_PRE_POST = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_PRE_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% %%% train POST___ test PRE
% 
% curr_transmat = transmatPOST;
% curr_lambda = lambdaPOST;
% 
% 
% trainPeriod = 'POST';
% testPeriod = 'PRE'; 
% 
% HMMprctile_POST_PRE = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_POST_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% 
% %%% train POST___ test RUN
% 
% curr_transmat = transmatPOST;
% curr_lambda = lambdaPOST;
% 
% 
% trainPeriod = 'POST';
% testPeriod = 'RUN'; 
% 
% HMMprctile_POST_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);
% 
% % _____________surrogates
% HMMprctile_POST_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);
% 
% 
% 
% %%
% figure;
% set(gcf, 'position', [129 64 2408 1232])
% 
% % Cross-validations (within-periods)
% 
% subplot(2,6,1)
% congruencePlot2(HMMprctilePRE, HMMprctilePRE_pts, 1, sprintf('%s Cross-validation', 'PRE'))
% 
% subplot(2,6,2)
% congruencePlot2(HMMprctileRUN, HMMprctileRUN_pts, 1, sprintf('%s Cross-validation', 'RUN'))
% 
% subplot(2,6,3)
% congruencePlot2(HMMprctilePOST, HMMprctilePOST_pts, 1, sprintf('%s Cross-validation', 'POST'))
% 
% 
% 
% % Congruence of PBEs with other periods' models
% 
% subplot(2,6,7)
% congruencePlot2(HMMprctile_PRE_RUN, HMMprctile_PRE_RUN_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'RUN'))
% 
% subplot(2,6,8)
% congruencePlot2(HMMprctile_RUN_POST, HMMprctile_RUN_POST_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'POST'))
% 
% subplot(2,6,9)
% congruencePlot2(HMMprctile_PRE_POST, HMMprctile_PRE_POST_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'POST'))
% 
% subplot(2,6,10)
% congruencePlot2(HMMprctile_RUN_PRE, HMMprctile_RUN_PRE_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'PRE'))
% 
% subplot(2,6,11)
% congruencePlot2(HMMprctile_POST_PRE, HMMprctile_POST_PRE_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'PRE'))
% 
% subplot(2,6,12)
% congruencePlot2(HMMprctile_POST_RUN, HMMprctile_POST_RUN_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'RUN'))
% 
% 
% 
% 
% %%
% 
% currDir = pwd;
% FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/congruence(wi_and_bw)/'];
% mkdir(FileBase)
% 
% 
% print(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
% savepdf(gcf, [FileBase fileinfo.name '_congruence'])
% 
% 
% save([FileBase fileinfo.name '_congruence.mat'], ...
%     'HMMprctilePRE', 'HMMprctilePRE_pts', ...
%     'HMMprctileRUN', 'HMMprctileRUN_pts', ...
%     'HMMprctilePOST', 'HMMprctilePOST_pts', ...
%     'HMMprctile_PRE_RUN', 'HMMprctile_PRE_RUN_pts', ...
%     'HMMprctile_RUN_POST', 'HMMprctile_RUN_POST_pts', ...
%     'HMMprctile_PRE_POST', 'HMMprctile_PRE_POST_pts', ...
%     'HMMprctile_RUN_PRE', 'HMMprctile_RUN_PRE_pts', ...
%     'HMMprctile_POST_PRE', 'HMMprctile_POST_PRE_pts', ...
%     'HMMprctile_POST_RUN', 'HMMprctile_POST_RUN_pts')
%  
%  
% 
% %% Joint distribution (histogram) of congruence score (for example the congruence of POST PBEs with PRE and POST models)
% 
% 
% % Joint distribution of congruence 
% % Whether RUN signifant sequences are congruent with model traind of PRE
% % PBEs
% 
% 
% figure;
% set(gcf, 'position', [796 125 1658 1184])
% 
% 
% % plot 2D histogram
% 
% data = [HMMprctile_PRE_RUN HMMprctileRUN];
% h = hist3(data, [10 10])/length(data);
% 
% subplot(2,2,1)
% 
% imagesc(h); set(gca, 'YDir', 'normal'); 
% colormap(jet)
% % colormap(flipud(colormap('hot')))
% set(gca, 'fontsize', 16, 'linewidth', 3)
% 
% set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% 
% xlabel({'Conguence of RUN PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
% ylabel({'Conguence of RUN PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)
% 
% title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)
% 
% 
% axis square
% colorbar
% 
% 
% 
% % plot stacked histogram
% 
% runIdx = HMMprctileRUN > 0.95;
% prerunIdx = HMMprctile_PRE_RUN > 0.95;
% 
% noncong = [length(find(~runIdx & ~prerunIdx));
%            length(find(~runIdx & prerunIdx))];
%        
% cong = [length(find(runIdx & ~prerunIdx));
%         length(find(runIdx & prerunIdx))];      
% 
%        
% noncong_cumcount = cumsum(noncong);
% noncong_cumcount = noncong_cumcount (end:-1:1);
% 
% cong_cumcount = cumsum(cong);
% cong_cumcount = cong_cumcount(end:-1:1);
% 
% thecolors =  [50 150 255 ;255 255 150]/255;
% 
% 
% subplot(2,2,2)
% 
% 
% hold on
% for ii = 1:2
%     bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
% end
% hold off
% 
% 
% set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
% set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% 
% xlabel('Congruence with RUN model', 'fontsize', 16)
% ylabel('Number of PBEs', 'fontsize', 20)
% 
% % axis square
% 
% legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
% legend boxoff 
% 
% 
% 
% 
% 
% 
% % Joint distribution of Congruence of post PBEs with PRE and RUN models
% 
% 
% % plot 2D histogram
% 
% data = [HMMprctile_PRE_POST HMMprctile_RUN_POST];
% h = hist3(data, [10 10])/length(data);
% 
% subplot(2,2,3)
% 
% imagesc(h); set(gca, 'YDir', 'normal'); 
% colormap(jet)
% % colormap(flipud(colormap('hot')))
% set(gca, 'fontsize', 16, 'linewidth', 3)
% 
% set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% 
% xlabel({'Conguence of POST PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
% ylabel({'Conguence of POST PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)
% 
% title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)
% 
% 
% axis square
% colorbar
% 
% 
% % plot stacked histogram
% 
% runIdx = HMMprctile_RUN_POST > 0.95;
% prerunIdx = HMMprctile_PRE_POST > 0.95;
% 
% noncong = [length(find(~runIdx & ~prerunIdx));
%            length(find(~runIdx & prerunIdx))];
%        
% cong = [length(find(runIdx & ~prerunIdx));
%         length(find(runIdx & prerunIdx))];      
% 
%        
% noncong_cumcount = cumsum(noncong);
% noncong_cumcount = noncong_cumcount (end:-1:1);
% 
% cong_cumcount = cumsum(cong);
% cong_cumcount = cong_cumcount(end:-1:1);
% 
% thecolors =  [50 150 255 ;255 255 150]/255;
% 
% 
% subplot(2,2,4)
% 
% 
% hold on
% for ii = 1:2
%     bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
% end
% hold off
% 
% 
% set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
% set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% 
% xlabel('Congruence of POST with RUN model', 'fontsize', 16)
% ylabel('Number of PBEs', 'fontsize', 20)
% 
% % axis square
% 
% legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
% legend boxoff 
% 
% 
% print(gcf, [FileBase fileinfo.name '_Overlap_of_wibw'], '-dpng')
% savepdf(gcf, [FileBase fileinfo.name '_Overlap_of_wibw'])
% 
% 

%% Correspondence between the states from models trained on different periods(PRE, RUN, POST) 

% The idea is to train models separately on different periods,
% then decode the probability over states 


% concatenate all the time bins within each period

cntPBES_PRE = cell2mat(eventsBinnedfiringPRE(:,2)');
cntPBES_RUN = cell2mat(eventsBinnedfiringRUN(:,2)');
cntPBES_POST = cell2mat(eventsBinnedfiringPOST(:,2)');




% if we need shuffle data ...
% 
% cntPBES_PRE = cntPBES_PRE(:, randperm(size(cntPBES_PRE, 2)));
% cntPBES_RUN = cntPBES_RUN(:, randperm(size(cntPBES_RUN, 2)));
% cntPBES_POST = cntPBES_POST(:, randperm(size(cntPBES_POST, 2)));


uniformTransMat = 0;
normalizeROW = 0;



% Decode in PRE

[PRE_RUN_testonPRE, sortPRE_RUN_PRE] = stateCorrespondence(cntPBES_PRE, lambdaPRE, transmatPRE, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);
[RUN_PRE_testonPRE, sortRUN_PRE_PRE] = stateCorrespondence(cntPBES_PRE, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);

[RUN_POST_testonPRE, sortRUN_POST_PRE] = stateCorrespondence(cntPBES_PRE, lambdaRUN, transmatRUN, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_RUN_testonPRE, sortPOST_RUN_PRE] = stateCorrespondence(cntPBES_PRE, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);

[PRE_POST_testonPRE, sortPRE_POST_PRE] = stateCorrespondence(cntPBES_PRE, lambdaPRE, transmatPRE, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_PRE_testonPRE, sortPOST_PRE_PRE] = stateCorrespondence(cntPBES_PRE, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);



% Decode in RUN

[PRE_RUN_testonRUN, sortPRE_RUN_RUN] = stateCorrespondence(cntPBES_RUN, lambdaPRE, transmatPRE, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);
[RUN_PRE_testonRUN, sortRUN_PRE_RUN] = stateCorrespondence(cntPBES_RUN, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);

[RUN_POST_testonRUN, sortRUN_POST_RUN] = stateCorrespondence(cntPBES_RUN, lambdaRUN, transmatRUN, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_RUN_testonRUN, sortPOST_RUN_RUN] = stateCorrespondence(cntPBES_RUN, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);

[PRE_POST_testonRUN, sortPRE_POST_RUN] = stateCorrespondence(cntPBES_RUN, lambdaPRE, transmatPRE, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_PRE_testonRUN, sortPOST_PRE_RUN] = stateCorrespondence(cntPBES_RUN, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);



% Decode in POST

[PRE_RUN_testonPOST, sortPRE_RUN_POST] = stateCorrespondence(cntPBES_POST, lambdaPRE, transmatPRE, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);
[RUN_PRE_testonPOST, sortRUN_PRE_POST] = stateCorrespondence(cntPBES_POST, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);

[RUN_POST_testonPOST, sortRUN_POST_POST] = stateCorrespondence(cntPBES_POST, lambdaRUN, transmatRUN, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_RUN_testonPOST, sortPOST_RUN_POST] = stateCorrespondence(cntPBES_POST, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, normalizeROW);

[PRE_POST_testonPOST, sortPRE_POST_POST] = stateCorrespondence(cntPBES_POST, lambdaPRE, transmatPRE, lambdaPOST, transmatPOST, uniformTransMat, normalizeROW);
[POST_PRE_testonPOST, sortPOST_PRE_POST] = stateCorrespondence(cntPBES_POST, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, normalizeROW);



% temp = [PRE_RUN_testonPRE(:); RUN_POST_testonPRE(:);  PRE_POST_testonPRE(:); PRE_RUN_testonRUN(:); RUN_POST_testonRUN(:); PRE_POST_testonRUN(:); PRE_RUN_testonPOST(:); RUN_POST_testonPOST(:); PRE_POST_testonPOST(:)];
% 
% bottom = min(temp);
% top = max(temp);



%%

figure;
set(gcf, 'position', [90 416 2408 877])


% PRE
subplot(3,6,1); plotPanel(PRE_RUN_testonPRE, sortPRE_RUN_PRE, 'PRE', 'RUN', 'PRE') 

subplot(3,6,2); plotPanel(RUN_PRE_testonPRE, sortRUN_PRE_PRE, 'RUN', 'PRE', 'PRE')

subplot(3,6,3); plotPanel(RUN_POST_testonPRE, sortRUN_POST_PRE, 'RUN', 'POST', 'PRE')

subplot(3,6,4); plotPanel(POST_RUN_testonPRE, sortPOST_RUN_PRE, 'POST', 'RUN', 'PRE')

subplot(3,6,5); plotPanel(PRE_POST_testonPRE, sortPRE_POST_PRE, 'PRE', 'POST', 'PRE')

subplot(3,6,6); plotPanel(POST_PRE_testonPRE, sortPOST_PRE_PRE, 'POST', 'PRE', 'PRE')


% RUN
subplot(3,6,7); plotPanel(PRE_RUN_testonRUN, sortPRE_RUN_RUN , 'PRE', 'RUN', 'RUN')

subplot(3,6,8); plotPanel(RUN_PRE_testonRUN, sortRUN_PRE_RUN, 'RUN', 'PRE', 'RUN')

subplot(3,6,9); plotPanel(RUN_POST_testonRUN, sortRUN_POST_RUN, 'RUN', 'POST', 'RUN')

subplot(3,6,10); plotPanel(POST_RUN_testonRUN, sortPOST_RUN_RUN, 'POST', 'RUN', 'RUN')

subplot(3,6,11); plotPanel(PRE_POST_testonRUN, sortPRE_POST_RUN, 'PRE', 'POST', 'RUN')

subplot(3,6,12); plotPanel(POST_PRE_testonRUN, sortPOST_PRE_RUN, 'POST', 'PRE', 'RUN')


% POST
subplot(3,6,13); plotPanel(PRE_RUN_testonPOST, sortPRE_RUN_POST, 'PRE', 'RUN', 'POST')

subplot(3,6,14); plotPanel(RUN_PRE_testonPOST, sortRUN_PRE_POST, 'RUN', 'PRE', 'POST')

subplot(3,6,15); plotPanel(RUN_POST_testonPOST, sortRUN_POST_POST, 'RUN', 'POST', 'POST')

subplot(3,6,16); plotPanel(POST_RUN_testonPOST, sortPOST_RUN_POST, 'POST', 'RUN', 'POST')

subplot(3,6,17); plotPanel(PRE_POST_testonPOST, sortPRE_POST_POST, 'PRE', 'POST', 'POST')

subplot(3,6,18); plotPanel(POST_PRE_testonPOST, sortPOST_PRE_POST, 'POST', 'PRE', 'POST')





%%
currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/StateCorrespondence/'];
mkdir(FileBase)

set(gcf,'Units','points');
pos = get(gcf,'Position');

set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(gcf, [FileBase 'correspondenceBetweenStates'],'-dpng','-r0')

%% Calculating the sequence score over time 
% 
% 
% % bvrTimeList = behavior.list(:,[1 2])/Fs;
% % bvrState = behavior.list(:,3);
% 
% 
% 
% 
% figure;
% set(gca, 'position', [995 783 796 538])
% hold on
% 
% 
% 
% % (1) time courses during RUN (for PRE and RUN(CV))
% 
% 
% % tbeginRUN = behavior.time(2,1)/Fs;
% % tendRUN = behavior.time(2,2)/Fs;
% % 
% % PBEtimesRUN = secondaryPBEs_run(:,1)/Fs;
% % 
% % 
% % 
% % %QW periods during RUN
% % 
% % qw = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% % qwRUN = qw(qw(:,1) > tbeginRUN & qw(:,2) < tendRUN, :);
% % qwRUN = qwRUN - tbeginRUN;
% 
% % 
% % subplot(2,1,1)
% % PlotSeqFreqvsTime(PBEtimesRUN, HMMprctile_RUN, HMMprctile_PRE_RUN, 'RUN', 'PRE', 'RUN', qwRUN, 12, tbeginRUN, tendRUN)
% 
% 
% 
% % (2) time courses during RUN (for PRE and RUN(CV))
% Fs = fileinfo.Fs;
% 
% tbeginPOST = behavior.time(3,1)/Fs;
% tendPOST = behavior.time(3,2)/Fs;
% 
% PBEtimesPOST = secondaryPBEs_post(:,1)/Fs;
% 
% %nrem periods during POST
% 
% % nrem = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% 
% 
% nrem = [behavior.Drowsy;behavior.NREM;behavior.Intermediate]/Fs;
% nremPOST = nrem(nrem(:,1) > tbeginPOST & nrem(:,2) < tendPOST, :);
% nremPOST = nremPOST - tbeginPOST;
% 
% %
% subplot(2,1,1)
% PlotSeqFreqvsTime(PBEtimesPOST, HMMprctile_RUN_POST, HMMprctile_PRE_POST, 'RUN', 'PRE', 'POST', nremPOST, 12, tbeginPOST, tendPOST)
% 
% 
% FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceScoreOverTime/'];
% mkdir(FileBase)
% 
% 
% 
% savepdf(gcf, [FileBase '/' fileinfo.name '_seqvstime'])
% print(gcf, [FileBase '/' fileinfo.name  '_seqvstime'], '-dpng')
% 


end






%% Functions

function PlotSeqFreqvsTime(PBEtimes, PBEscores1, PBEscores2, trainPeriod1, trainPeriod2, testPeriod, offPeriods, noEpochs, tbegin, tend)

% noEpochs = 12; % epochs with duration of ~ 15 min

epochs = linspace(tbegin, tend, noEpochs+1);
epochDur = epochs(2) - epochs(1);

epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;

[~, epochIdx] = histc(PBEtimes, epochs);


% model #1
epochMean1 = zeros(noEpochs, 1);
epochErr1 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores1(epochIdx == epoch);
    epochMean1(epoch) = mean(epochPrctls);
    epochErr1(epoch) = std(epochPrctls)/length(epochPrctls);
end





% model #2
epochMean2 = zeros(noEpochs, 1);
epochErr2 = zeros(noEpochs, 1);

for epoch = 1 : noEpochs
    epochPrctls = PBEscores2(epochIdx == epoch);
    epochMean2(epoch) = mean(epochPrctls);
    epochErr2(epoch) = std(epochPrctls)/length(epochPrctls);
end



hold on

minY = floor(min([epochMean1; epochMean2])*10)/10;
maxY = ceil(max([epochMean1; epochMean2])*10)/10;

for ii = 1:length(offPeriods)
    patch([offPeriods(ii,1) offPeriods(ii,2) offPeriods(ii,2) offPeriods(ii,1)]/3600, [minY minY maxY maxY], ...
           [204 229 255]/255, 'EdgeColor', 'none');
%         set(p,'FaceAlpha', 0.5)
end

h1 = errorbar(epochCenters/3600, epochMean1, epochErr1, 'color', 'k', 'linewidth', 2);

h2 = errorbar(epochCenters/3600, epochMean2, epochErr2, 'color', 'r', 'linewidth', 2);




legend([h1, h2], ['train on ' trainPeriod1], ['train on ' trainPeriod2], 'Location','southeast')
legend boxoff 

set(gca, 'fontsize', 12, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'box', 'off', 'Layer', 'Top')
xlim([0 3])
ylim()

xlabel(['Time in ' testPeriod '(hr)'], 'fontsize', 14)
ylabel('PBE sequence score', 'fontsize', 12)



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


% smoothing the speed

sigma = 20; %%% smoothing the speed, the positon sampling rate is around 40 Hz (0.0256 sec sampling period), so duration of the sigma is about 25.6 ms times sigma (25.6*20 ~ 512 ms)
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
velocity = conv(velocity, smoothwin, 'same'); 


speed.v = velocity;
speed.t = timepnts;


end


function [eventsBinnedfiring, secondaryPBEs, idx_acceptedEvents] = finalBinningResult(pbEvents2, mySpikes, fileinfo)

%% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events


%%% Calculating the number of firing pyramidal units within each event (this is redundant; This calculation was merged into the binnig function)

qclus = [1 2 3]; % pyramidal and interneurons

% noFiringUnits = zeros(size(pbEvents2, 1), 1);
% for evt = 1 : size(pbEvents2, 1)
%     
%     evtSpkInd = find(mySpikes.t > pbEvents2(evt, 1) & mySpikes.t < pbEvents2(evt, 2) & ismember(mySpikes.qclu, qclus));
%     noFiringUnits(evt) = length(unique(mySpikes.unit(evtSpkInd)));
% end


%%% binning the spikes within each event

binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) 
[eventsBinnedfiring2, noFiringUnits] = timeBinning(pbEvents2, mySpikes, qclus, binDur, fileinfo);


% Remove the flanking zero-firing periods(silent bins) from the beginnig
% and end of each event. Since the interneurons were involved in calculating
% the PBEs' boundaries and here we are binning just pyramidal units' spike
% trains, some silent bins are expected. 

[eventsBinnedfiring, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring2, pbEvents2(:, 1), binDur, fileinfo);

% Qualify the events based on the minumum active units and length criteria
% Considering only the event with number of active units firing above 10 percent of the total or 5 (whichever is greater) and
% duration of at least 4 20ms-time bins.


% activeUnits = unique(mySpikes.unit);
activeUnits = 1:size(eventsBinnedfiring2{1,1}, 1);

% idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4 & eventLen <= 25)); % maximum length 500 ms (???)
idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4 & eventLen <= 25)); % maximum length 500 ms (???)
secondaryPBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:4)]; 


eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


end


function [poissonEventsBinnedfiring, timeSwapEventsBinnedfiring, temporalEventsBinnedfiring] = genSurrogates(eventsBinnedfiring)


%% Generating shuffle surrogate PBE datasets

% generate poisson simulated dataset

binDur = 0.02;
noEvents = size(eventsBinnedfiring, 1);

poissonEventsBinnedfiring = poissonSpikeTrain(eventsBinnedfiring, binDur); 
% cmprPoisson2Raw(poissonEventsBinnedfiring, eventsBinnedfiring, longORshort, fileinfo, binDur); %%% compare the poisson simulated data with the actual data distribution

% generate time swap dataset (coherent shuffle within each event, keeping the cofirings the same)

timeSwapEventsBinnedfiring = timeswap(eventsBinnedfiring, binDur);


% time bins temporal shuffle (independently for each unit, incoherent shuffle)

temporalEventsBinnedfiring = cell(noEvents, 1);
for evt = 1 : noEvents
    
    currEvt = eventsBinnedfiring(evt, :);
    
    randShift = randi(size(currEvt{2}, 2)-1, 1, size(currEvt{2}, 1)); % generate shift amount for each row
    
    
    % Because we are doing row cycle shuffle, we need to transpose the data
    % matrix and after doing column cycle shuffles transpose it again
    
    temporalEventsBinnedfiring{evt, 2} = column_cycle_shuffle(currEvt{2}', randShift)';
    temporalEventsBinnedfiring{evt, 1} = column_cycle_shuffle(currEvt{1}', randShift*20)'; % note the shift here
    
end

end



function HMMprctile = modelCongruence(testData, transmat, lambda, noShuffle, fileinfo, trainPeriod, testPeriod, ifShuffle)




cnctPBEs = cell2mat(testData(:,2)');

% resample each surrogate PBE (or better to say the bins) from the concatenated PBEs

rndgen = randperm(size(cnctPBEs, 2));
shuffled_cnctPBEs = cnctPBEs(:, rndgen); 



noEvents = size(testData, 1);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, noShuffle); %% log likelihood of shuffle data



HMMprctile = zeros(noEvents, 1);


setaside = [];
shuffleTransmats = genshuffles(transmat, noShuffle, setaside);



for evt = 1 : noEvents
    
    currEvent = testData{evt, 2};
    noTimebins = size(currEvent, 2);
    
    if ifShuffle
        currEvent = shuffled_cnctPBEs(:, 1:noTimebins);
        shuffled_cnctPBEs(:, 1:noTimebins) = [];
    end
    
      
    B = poisson_prob(currEvent, lambda,1);
    prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the actual data
    [~, ~, ~,  dataLL(evt), ~] = fwdback(prior, transmat, B);
    
    %%%% null distribution of observation likelihhod
    %%%% shuffling the transmat matrix

%     shuffled_transmat = zeros(size(transmat));

    for sn = 1 : noShuffle

%         for s = 1 : numofStates
% 
%             %%% redistribute the transition probabilities within each 
%             %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
% 
%             randgen = randperm(numofStates)';
%             randgen(randgen == s) = [];
% 
%             shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 
% 
% %             shuffled_transmat(s,:) = transmat(s, shuffleInd); 
%             
%             shuffled_transmat(:,s) = transmat(shuffleInd, s);
% 
%         end             

%         [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffled_transmat, B); %%% B and Prior are the same as in the case of raw model
        [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffleTransmats(:,:,sn), B); %%% B and Prior are the same as in the case of raw model
        
    end




    for method = 1%:2

        HMMprctile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/noShuffle;

    end
    
    
    if ifShuffle
        dataType = 'pooledTimeSwap';
    else
        dataType = 'actual';
    end
    
    
    if mod(evt, 20) == 0
        fprintf(1, [trainPeriod '-' testPeriod '-' dataType '__event %d, HMMPercentile = %f\n'], evt, HMMprctile(evt, 1));
    end
    
    
end



currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/bw_periods_congruence/'];
mkdir(FileBase)




save([FileBase '/Train_' trainPeriod '___Test_' testPeriod '__' dataType '.mat'], 'HMMprctile', 'dataLL', 'nullLL')

end


function congruencePlot(HMMprctile, HMMprctile_p, HMMprctile_ts, HMMprctile_tc, okLegend, curr_title)

hold on
[h, bins] = hist(HMMprctile, 100);
cumhist = cumsum(h);
threshBin = find(bins > 0.95, 1, 'first');
crosspnt = cumhist(threshBin);

h_p = hist(HMMprctile_p, bins);
cumhist_p = cumsum(h_p);
crosspnt_p = cumhist_p(threshBin);


h_ts = hist(HMMprctile_ts, bins);
cumhist_ts = cumsum(h_ts);
crosspnt_ts = cumhist_ts(threshBin);


h_tc = hist(HMMprctile_tc, bins);
cumhist_tc = cumsum(h_tc);
crosspnt_tc = cumhist_tc(threshBin);

p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9], 'linewidth', 4);

p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

p2 = plot(bins, cumhist_ts, 'color', [100 180 255]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_ts crosspnt_ts], 'color', [100 180 255]/255,'LineStyle', '--', 'linewidth', 2)

p3 = plot(bins, cumhist_p, 'color', [255 180 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_p crosspnt_p], 'color', [255 180 100]/255,'LineStyle', '--', 'linewidth', 2)

p4 = plot(bins, cumhist_tc, 'color', [255 100 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_tc crosspnt_tc], 'color', [255 100 100]/255,'LineStyle', '--', 'linewidth', 2)



noEvents = length(HMMprctile);
line([bins(threshBin) bins(threshBin)], [0 noEvents], 'color', 'k','LineStyle', '--', 'linewidth', 2)


if okLegend
    legend([p1 p2 p3 p4 p0],{'actual','time swap','poisson','temporal shuffle', 'chance line'}, 'Location', 'northwest')
    legend boxoff 
end

set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])
xlabel('percentile score', 'fontsize', 16)

if okLegend
    ylabel('cumulative number of PBEs', 'fontsize', 16)
end

title(curr_title, 'fontsize', 14)

ylim([0 length(HMMprctile)])
axis square

end


function [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo)


%% binning the run data


qclus = [1 2 3]; % pyramidal and interneurons

runBinnedfiring2 = timeBinning(activePeriods , mySpikes, qclus, binDur, fileinfo.Fs); % the theta periods


% Calculate the track positions corresponding to each time bin

% periods = activePeriods * 1/fileinfo.lfpSampleRate;

periods = activePeriods;

runBinPos = cell(size(periods, 1), 1);

for ii = 1 : size(periods, 1)
    
%     runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/fileinfo.lfpSampleRate, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3), fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin

    runBinPos{ii}(end, :) = [];
end
    


% runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields

% posBinSize = 0.5; % in cm3
runidx = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2));

xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);


% noXPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
% noYPosBins = floor((max(ypos) - min(ypos))/posBinSize);

xposBins = min(xpos): posBinSize: max(xpos); % center of the position bins
xposBins(end) = max(xpos);

xposcenters = xposBins + posBinSize/2;
xposcenters(end) = [];


yposBins = min(ypos): posBinSize: max(ypos); 
yposBins(end) = max(ypos);

yposcenters = yposBins + posBinSize/2;
yposcenters(end) = [];


posbinIdx2 = cell(size(periods, 1), 1);
for ii = 1 : size(periods, 1)
    
    [~, posbinIdx2{ii}(:,1)] = histc(runBinPos{ii}(:,1), xposBins);
    [~, posbinIdx2{ii}(:,2)] = histc(runBinPos{ii}(:,2), yposBins);
end


runBinnedfiring{1,1} = cell2mat(runBinnedfiring2(:,1)'); % concatenating all the run binned spikes
runBinnedfiring{1,2} = cell2mat(runBinnedfiring2(:,2)'); % the same as the first column but with milisecond bins

posbinIdx = cell2mat(posbinIdx2);


end


function lsPFs = plot_lsPF(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod)


% numofStates = 40;
lsPFs = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmat, lambda, posbinIdx, [], numofStates, xposcenters, yposcenters);


% gamma_avg = permute(gamma_avg, [2 1 3]); 



% plot the lsPFs

figure;

numStatesperRow = 3;

runidx = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2));
xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);

for ii = 1: numofStates
    
    subplot(ceil(numofStates/numStatesperRow), numStatesperRow, ii)
    
    plot(fileinfo.xyt(runidx, 1), fileinfo.xyt(runidx, 2), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
    
    
    hold on
    h1 = imagesc(xposcenters, yposcenters, lsPFs(:,:, ii));
    set(h1, 'AlphaData', 0.7)

    
    colormap('hot')
    set(gca,'YDir','normal')
    
    cornerstate = (ceil(numofStates/numStatesperRow)-1) * numStatesperRow + 1;
    
    if ii == cornerstate
        
        set(gca, 'XTick', floor([min(xpos) max(xpos)]), 'YTick', floor([min(ypos) max(ypos)]), 'box', 'off')
        xlabel('X (cm)', 'fontsize', 14)
        ylabel('Y (cm)', 'fontsize', 14)
        
    else
        
        set(gca, 'XTick', [], 'YTick', [], 'box', 'off')
    end

    xlim(floor([min(xpos) max(xpos)]))
    ylim(floor([min(ypos) max(ypos)]))

end

% suptitle([fileinfo.name '-' pbePeriod])
set(gcf, 'position', [1 1 800 1500])

end



function [PosDecError, PosDecError_shuffledlsPFs, actualPositions, decodedPositions_actual, probabilityOverPosition] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates)
% [PosDecError, PosDecError_shuffledOrder, PosDecError_shuffledlsPFs] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates)


excludeBins = find(posbinIdx(:,1).* posbinIdx(:,2) == 0);

posbinIdx(excludeBins, :) = [];
runBinnedfiring(:, excludeBins) = [];


noTimeBins = length(posbinIdx);


% 5-fold cross-validation

noFolds = 5;
foldSize = floor(noTimeBins/noFolds);


decodedPositions_actual = zeros(length(posbinIdx), 2);
% decodedPositions_shuffledOrder = zeros(length(posbinIdx), 2);
decodedPositions_shuffledPosition = zeros(length(posbinIdx), 2);

probabilityOverPosition = zeros(length(yposcenters), length(xposcenters), length(posbinIdx));


for k = 1:noFolds
    
    
    if k == noFolds
        testSet = (k-1)*foldSize+1 : noTimeBins;
    else
        testSet = (k-1)*foldSize+1 : k*foldSize;
    end
    
    trainSet = setdiff(1:noTimeBins, testSet);
    
    
    %% calculate lsPFs for the train set
    
    
    train_timeBins = runBinnedfiring(:,trainSet);
    
    train_posbinIdx = posbinIdx(trainSet, :); % x and y coordiantes in columns
    
    train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_posbinIdx, [], numofStates, xposcenters, yposcenters);
    
    
    train_shuffle_posbinIdx = train_posbinIdx(randperm(length(trainSet)), :);
    
    shuffle_train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_shuffle_posbinIdx, [], numofStates, xposcenters, yposcenters);
    
    
    % Decode the states in RUN data 
    
    test_timeBins = runBinnedfiring(:,testSet);
    
    obsProb = poisson_prob(test_timeBins, lambda,1);
      
    prior = 1/numofStates * ones(numofStates,1); % uniform prior
    
    
    %(1) Decoding track positions based on actual lsPFs
    
    decodedPositions_actual(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters, yposcenters); 
    
%     probabilityOverPosition(:,:,testSet) = rr;
    %--- shuffled time bin orders
    
%     decodedPositions_shuffledOrder(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledOrder'); 
    
    
    %--- shuffled gamma: it supposed to take away all the order and
    %coactivity information regarding the states
    
%     decodedPositions_shuffledPosition(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledPosition'); 
    

    %(2) Decoding track positions based on shuffled lsPFs (lsPFs calculated after shuffling the position indices)
    
    decodedPositions_shuffledPosition(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, shuffle_train_lsPFs, xposcenters, yposcenters); 
    

end


posbinIdx(find(posbinIdx(:,2) > length(yposcenters)), 2) = length(yposcenters);
posbinIdx(find(posbinIdx(:,1) > length(xposcenters)), 1) = length(xposcenters);


% excludeBins = find(posbinIdx(:,1).* posbinIdx(:,2) == 0);
% 
% posbinIdx(excludeBins, :) = [];
% 

% decodedPositions_actual(excludeBins, :) = [];
% decodedPositions_shuffledPosition(excludeBins, :) = [];


actualPositions = [yposcenters(posbinIdx(:,2))' xposcenters(posbinIdx(:,1))'];


PosDecError = sqrt((decodedPositions_actual(:, 1) - actualPositions(:, 1)).^2  + (decodedPositions_actual(:, 2) - actualPositions(:, 2)).^2);
% PosDecError(isnan(PosDecError)) = [];


% PosDecError_shuffledOrder = sqrt((decodedPositions_shuffledOrder(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledOrder(:, 2) - xposcenters(posbinIdx(:,1))').^2);

PosDecError_shuffledlsPFs = sqrt((decodedPositions_shuffledPosition(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledPosition(:, 2) - xposcenters(posbinIdx(:,1))').^2);    
% PosDecError_shuffledlsPFs(isnan(PosDecError_shuffledlsPFs)) = [];



end


function decodedPositions = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters)
% decodedPositions = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, shuffleMode)


noTestBins = size(obsProb, 2);
% numofStates = size(transmat, 1);
% shuffledOrder = randperm(noTestBins);

    
% if strcmp(shuffleMode, 'shuffledOrder')
%     [~, ~, gamma(:, shuffledOrder)] = fwdback(prior, transmat, obsProb(:, shuffledOrder));
% else
%     [~, ~, gamma] = fwdback(prior, transmat, obsProb);
% end

[~, ~, gamma] = fwdback(prior, transmat, obsProb);

decodedPositions = zeros(noTestBins, 2);
for jj = 1: noTestBins
    
        
%     if strcmp(shuffleMode, 'shuffledPosition')
%         stateProbability = gamma(randperm(numofStates), jj);
%     else
%         stateProbability = gamma(:, jj);
%     end
    
    stateProbability = gamma(:, jj);

    probabilityOverPosition = sum(repmat(permute(stateProbability, [3 2 1]), length(yposcenters), length(xposcenters)) .* train_lsPFs, 3);

    decodedPositions(jj, 1) = yposcenters(ceil(sum((1:length(yposcenters))' .* sum(probabilityOverPosition, 2))));
    decodedPositions(jj, 2) = xposcenters(ceil(sum((1:length(xposcenters)) .* sum(probabilityOverPosition, 1))));
end


end


function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-depsc','-r0')

end


function plotErrorCdf(data, shuffledPosition, ylabelneeded)
% plotErrorCdf(data, shuffledOrder, shuffledPosition, ylabelneeded)

hold on

% actual data
[h1, bins] = hist(data, 100);
h1c = cumsum(h1)/sum(h1);

dataMed = median(data);
 
curve1 = plot(bins, h1c, 'linewidth', 4, 'color', [150, 150, 255]/255);
line([dataMed dataMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
line([0 dataMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
text(dataMed, 0.03, sprintf('%.1f cm', dataMed), 'color', [150, 150, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')

% % shuffled order
% [h2, bins] = hist(shuffledOrder, 100);
% h2c = cumsum(h2)/sum(h2);
% 
% shuffleMed = median(shuffledOrder);
% 
% curve2 = plot(bins, h2c, 'linewidth', 4, 'color', [150, 200, 255]/255);
% line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
% line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
% text(shuffleMed, 0.07, sprintf('%.1f cm', shuffleMed), 'color', [150, 200, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')
% 

% shuffled Gamma
[h3, bins] = hist(shuffledPosition, 100);
h3c = cumsum(h3)/sum(h3);

shuffleMed = median(shuffledPosition);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(shuffleMed, 0.1, sprintf('%.1f cm', shuffleMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')




xlim([0 max([data; shuffledPosition])])
set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])


xlabel('Decoding Error(cm)', 'fontsize', 16)

if ylabelneeded
    ylabel('Cumulative ratio', 'fontsize', 24)
    
%     legend('actual', 'shuffle state probability', 'shuffle position index', 'Location', 'northwest'); 
%     legend([curve1, curve2, curve3], 'actual', 'shuffled run bins', 'shuffled gamma dist.', 'Location', 'northwest'); 
        legend([curve1, curve3], 'actual', 'shuffled lsPFs', 'Location', 'northwest'); 


    legend boxoff 
end

axis square
end


function shuffleTransmats = genshuffles(transmat, noShuffle, setaside)


numofStates = size(transmat, 1);
shuffleTransmats = zeros(numofStates, numofStates, noShuffle);

for sn = 1 : noShuffle
    
    
    shuffle_transmat = zeros(numofStates, numofStates);
    for s = 1 : numofStates

        %%% redistribute the transition probabilities within each 
        %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
        %
        
        exclude = s;
        
%         exclude = [setaside s];
        initialIdx = setdiff(1:numofStates, exclude);

        shuffleSet = initialIdx(randperm(length(initialIdx)));

        shuffleInd = zeros(numofStates, 1);

        shuffleInd(exclude) = exclude; % keep the same indices for states within the exclude set
        shuffleInd(initialIdx) = shuffleSet; % use the randomized indices for the remained


%                 randgen = randperm(numofStates)';
%                 randgen(randgen == s) = [];
% 
%                 shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 

%         shuffleTransmats(s,:, sn) = transmat(s, shuffleInd); % shuffling within the rows

        shuffle_transmat(:, s) = transmat(shuffleInd, s); % shuffling within the columns
    end
    
    % if doing column-wise shuffle we need to normalize every row again to
    % sum up to one
    
    shuffleTransmats(:, :, sn) = shuffle_transmat./repmat(sum(shuffle_transmat, 2), [1, numofStates]);
    
end

end


function savePBEs(primaryPBEs, secondaryPBEs, acceptedEvtsIdx, analysisPeriod, fileinfo) %#ok<INUSL>


currDir = pwd;
Folderbase = [currDir '/' fileinfo.name '/PopulationBurstEvents/' analysisPeriod];
mkdir(Folderbase)

save([Folderbase '/' fileinfo.name  '-' 'PBEs.mat'], 'primaryPBEs', 'secondaryPBEs', 'acceptedEvtsIdx') 

Filename = [Folderbase '/' fileinfo.name '.prm.evt'];
MakeEvtFile(primaryPBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)

Filename = [Folderbase '/' fileinfo.name '.scn.evt'];
MakeEvtFile(secondaryPBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)


end


function HMMPercentile = HMMcongruence_crossvalid_pts(data, testEvts, transmat, lambda, noShuffle, dataType, fileinfo, period)


currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/cross-validation/' period '/'];

mkdir(FileBase)




noEvents = size(data, 1);
numofFolds = length(testEvts);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, noShuffle); %% log likelihood of shuffle data

% HMMScoreZ = zeros(noEvents, 1); %% 2 if we add again time swap
HMMPercentile = zeros(noEvents, 1);

% partialLikeli = cell(1, noEvents); 
% bins_nullDistPercentiles = cell(1, noEvents);
% binsCongSignificance = cell(1, noEvents);
% maxCongruenceLength = zeros(1, noEvents);
% totalCongruenceLength = zeros(1, noEvents);

for fold = 1 : numofFolds
    
    curr_transmat = transmat(:, :, fold);
    curr_lambda = lambda(:, :, fold);
    
    %
%     arrivalSparsity = ginicoeff(curr_transmat, 1);
%     departureSparsity = ginicoeff(curr_transmat, 2);
%     
%     diffSparsity = abs(departureSparsity' - arrivalSparsity);
%     
%     setaside = find(diffSparsity > prctile(diffSparsity, 75)); % keep aside the states with diffsparisties in the 4th quartile
%     
%     setaside = randi(numofStates, 1, length(setaside)); % if we are keeping aside some random states
    
    setaside = [];
    shuffleTransmats = genshuffles(curr_transmat, noShuffle, setaside);

    %
    
    % testing the model on shuffled test PBEs (pooled timeswap surrogate dataset)
    testPBEs = data(testEvts{fold}, :);
    
%     pts_testPBEs = cell(size(testPBEs));
%     pts_testPBEs(:,1) = []; % 1 ms bins, we don't need to compute it here
    
    cnctPBEs = cell2mat(testPBEs(:,2)');
    
    % resample each surrogate PBE (or better to say the bins) from the concatenated PBEs
    
    rndgen = randperm(size(cnctPBEs, 2));
    shuffled_cnctPBEs = cnctPBEs(:, rndgen); 
    
    for evt = testEvts{fold}

        
        currEvent = data{evt, 2};
        noTimebins = size(currEvent, 2);
        
        surrogateEvent = shuffled_cnctPBEs(:, 1:noTimebins);
        shuffled_cnctPBEs(:, 1:noTimebins) = []; % the array should exhaust when we are done with all of the events

        B = poisson_prob(surrogateEvent, curr_lambda,1);
        prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, ~,  dataLL(evt), ~] = fwdback(prior, curr_transmat, B);
        
        
        
        %%%% null distribution of observation likelihhod
        %%%% shuffling the transmat matrix
        
%         shuffled_transmat = zeros(size(curr_transmat));
        
        for sn = 1 : noShuffle

%             for s = 1 : numofStates
%                 
%                 %%% redistribute the transition probabilities within each 
%                 %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed
%                 %
%                 exclude = [setaside s];
%                 initialIdx = setdiff(1:numofStates, exclude);
%                 
%                 shuffleSet = initialIdx(randperm(length(initialIdx)));
%                 
%                 shuffleInd = zeros(numofStates, 1);
%                 
%                 shuffleInd(exclude) = exclude;
%                 shuffleInd(initialIdx) = shuffleSet;
%                 
%                 
% %                 randgen = randperm(numofStates)';
% %                 randgen(randgen == s) = [];
% % 
% %                 shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 
% 
%                 shuffled_transmat(s,:) = curr_transmat(s, shuffleInd); % shuffling within the rows
%                 
% %                 shuffled_transmat(:,s) = curr_transmat(shuffleInd, s); % shuffling within the columns
%             end             

%             [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffled_transmat, B); %%% B and Prior are the same as in the case of raw model
            [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffleTransmats(:, :, sn), B);
        end

        
%         
%         %%%%%%% shuffling the events (time swap)
% HMMPercentile
%         for sn = 1 : noShuffle
%             
%             %%% generating the time swapped data 
%             
%             rndgen = randperm(size(currEvent,2));
%             currEvent_sh = currEvent(:, rndgen); 
%     
% 
%             
%             B = poisson_prob(currEvent_sh, curr_lambda,1);
% 
%             [~, ~, ~, nullLL(evt,sn,2), ~] = fwdback(prior, curr_transmat, B);
% 
%         end
        
        
        for method = 1%:2

%             HMMScoreZ(evt, method) = rawScoreCompare2null(dataLL(evt), nullLL(evt,:,method)); 
            HMMPercentile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/noShuffle;

        end
        
        
        if mod(evt, 50) == 0
            fprintf(1, 'cross-validation_%s_event %d, HMMPercentile = %f\n', dataType, evt, HMMPercentile(evt, 1));
        end
        
        %%%% a time-resolved analysis
%         
%         percentile = 5;
%         
%         [partialLikeli{evt}, bins_nullDistPercentiles{evt}] = likelihoodProfile(currEvent, curr_transmat, curr_lambda, percentile, noShuffle);
%         
%         binsCongSignificance{evt} = (partialLikeli{evt}' - bins_nullDistPercentiles{evt}(:,2))./(bins_nullDistPercentiles{evt}(:,3) - bins_nullDistPercentiles{evt}(:,2)); 
%         
%         [maxCongruenceLength(evt), totalCongruenceLength(evt)] = measureCongLen(partialLikeli{evt}, bins_nullDistPercentiles{evt});
%         
    end

end

save([FileBase '/cvScores_' period '_' dataType '.mat'], 'HMMPercentile', 'dataLL', 'nullLL')

% save([FileBase '/' dataType '.mat'], 'HMMScoreZ', 'dataLL', 'nullLL', 'partialLikeli', 'bins_nullDistPercentiles', 'binsCongSignificance', 'maxCongruenceLength', 'totalCongruenceLength')


end



function congruencePlot2(HMMprctile, HMMprctile_pts,  okLegend, curr_title)

hold on

bins = 0:1:100;

h = hist(HMMprctile*100, bins);
h = h/sum(h)*100;

cumhist = 100-cumsum(h); %%  for the paper
threshBin = find(bins > 99, 1, 'first');
crosspnt = cumhist(threshBin);


h_pts = hist(HMMprctile_pts*100, bins);
h_pts = h_pts/sum(h_pts)*100;

cumhist_pts = 100-cumsum(h_pts); %% for the paper
crosspnt_pts = cumhist_pts(threshBin);



% p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9],
% 'linewidth', 4); % the chance line-what would be the chance distribution
% is not established yet

p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

p2 = plot(bins, cumhist_pts, 'color', 'b', 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt_pts crosspnt_pts], 'color', 'b','LineStyle', '--', 'linewidth', 2)


% noEvents = length(HMMprctile);
% line([bins(threshBin) bins(threshBin)], [0 1], 'color', 'k','LineStyle', '--', 'linewidth', 2)


if okLegend
    legend([p1 p2],{'actual', 'pooled time swap'}, 'Location', 'southwest')
    legend boxoff 
end

set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])
xlabel('percentile threshold', 'fontsize', 16)

set(gca, 'XTick', 0:20:100, 'YTick', 0:20:100, ...
    'XTickLabels', {'0', '', '', '', '', '100'}, 'YTickLabels', {'0', '', '', '', '', '100'})

if okLegend
    ylabel('% sig. PBEs', 'fontsize', 16)
end

title(curr_title, 'fontsize', 14)

ylim([-5 100])
xlim([-5 100])

axis square

end




function [A_B_corresp_norm,  sortAstates] = stateCorrespondence(testData, lambdaA, transmatA, lambdaB, transmatB, uniformTranprobs, normalizeROW)
        


% Decoding the states in the test period given model A and B (trained on the periods A and B, resp.)

noTimebins = size(testData, 2);

% Decoded states given model A
B = poisson_prob(testData, lambdaA, 1);

numofStatesA = size(lambdaA, 2);
prior = 1/numofStatesA * ones(numofStatesA,1);

if ~uniformTranprobs
    currTransmat = transmatA;
else
    currTransmat = 1/numofStatesA * ones(numofStatesA);% to make the transition probabilities irrelevant 
end
[~,~, gammagivenA] = fwdback(prior, currTransmat, B);


% Decoded states given model B
B = poisson_prob(testData, lambdaB, 1);

numofStatesB = size(lambdaB, 2);
prior = 1/numofStatesB * ones(numofStatesB,1);

if ~uniformTranprobs
    currTransmat = transmatB;
else
    currTransmat = 1/numofStatesB * ones(numofStatesB);
end
[~,~, gammagivenB] = fwdback(prior, currTransmat, B);

gammagivenA_chance = repmat(mean(gammagivenA, 2), [1, noTimebins]);
gammagivenB_chance = repmat(mean(gammagivenB, 2), [1, noTimebins]);


A_B_corresp = (gammagivenA * gammagivenB')/noTimebins; % A states in rows and B states in columns
A_B_corresp_chance = (gammagivenA_chance * gammagivenB_chance')/noTimebins;


A_B_corresp_norm = A_B_corresp - A_B_corresp_chance;


% Normalize (for each states from A what is the most correspondant state from B)

if normalizeROW
    A_B_corresp_norm = A_B_corresp_norm ./ repmat(sum(A_B_corresp_norm, 2), 1, size(A_B_corresp_norm, 2));
end


% sorting the states in one period based on the maximum corresponding
% states from anoter period

[~, maxBstates] = max(A_B_corresp_norm, [], 2);

[~, sortAstates] = sort(maxBstates);

end


function plotPanel(stateCorrespMat, rowSortInd, model, model2comp, decodingPeriod)


% 
% for t = 1:5
%     [i,j] = find(stateCorrespMat == max(stateCorrespMat(:)));
%     stateCorrespMat(i,j) = mean(stateCorrespMat(:));
% end
% cl = max([abs(prctile(stateCorrespMat, 99)) abs(prctile(stateCorrespMat, 1))]);
cl = 0.01;

imagesc(stateCorrespMat(rowSortInd, :), [-cl cl]); colormap('jet');

set(gca, 'YDir', 'normal', 'ytick', 1:length(rowSortInd), 'yticklabel', rowSortInd, 'fontsize', 16)
set(gca, 'XTick', [], 'YTick', [])

ylabel([model ' states'] , 'fontsize', 16)
xlabel([model2comp ' states'], 'fontsize', 16)

title(['test on ' decodingPeriod ' PBEs'], 'fontsize', 16, 'fontweight', 'normal')
colorbar

end




function [A_BC_corresp, sortAstates] = stateCorrespondence3(testData, lambdaA, transmatA, lambdaB, transmatB, lambdaC, transmatC, uniformTranprobs)


% Decoding the states in the test period given model A and B (trained on the periods A and B, resp.)


% Decoded states given model A
B = poisson_prob(testData, lambdaA, 1);

numofStatesA = size(lambdaA, 2);
prior = 1/numofStatesA * ones(numofStatesA,1);

if ~uniformTranprobs
    currTransmat = transmatA;
else
    currTransmat = 1/numofStatesA * ones(numofStatesA);% to make the transition probabilities irrelevant 
end
[~,~, gammagivenA] = fwdback(prior, currTransmat, B);


% Decoded states given model B
B = poisson_prob(testData, lambdaB, 1);

numofStatesB = size(lambdaB, 2);
prior = 1/numofStatesB * ones(numofStatesB,1);

if ~uniformTranprobs
    currTransmat = transmatB;
else
    currTransmat = 1/numofStatesB * ones(numofStatesB);
end
[~,~, gammagivenB] = fwdback(prior, currTransmat, B);



% Decoded states given model C
B = poisson_prob(testData, lambdaC, 1);

numofStatesC = size(lambdaC, 2);
prior = 1/numofStatesC * ones(numofStatesC,1);

if ~uniformTranprobs
    currTransmat = transmatC;
else
    currTransmat = 1/numofStatesC * ones(numofStatesC);% to make the transition probabilities irrelevant 
end
[~,~, gammagivenC] = fwdback(prior, currTransmat, B);


A_BC_corresp = gammagivenA * [gammagivenB; gammagivenC]'; 


% Normalize (for each states from A what is the most correspondant state from B)

A_BC_corresp = A_BC_corresp ./ repmat(sum(A_BC_corresp, 2), 1, size(A_BC_corresp, 2));


% sorting the states in one period based on the maximum corresponding
% states from anoter period

[~, maxBCstates] = max(A_BC_corresp, [], 2);

[~, sortAstates] = sort(maxBCstates);


end

