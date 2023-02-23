clear; clc; close all

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {1, 1};


for rr = 1%:length(rats) % loop for all the rats and sessions
    
rat = rats{rr};
sessionNumbers = allsessionNumbers{rr};
    
for sessionNumber = sessionNumbers

% 
% rat = input('\n Enter the animal name (in quotes) \n');
% rat(1) = upper(rat(1));
% 
% 
% sessionNumber = input('\n Enter the session number \n');


sessionName = [rat 'Maze' num2str(sessionNumber)];
sessionName2 = [rat '-maze' num2str(sessionNumber)];

VarList = {'spikes','behavior','position','speed','basics','ripple'};


for var = 1 : length(VarList)
    load([currDir '/pooled_includingKevin/wake-' VarList{var} '.mat'])
end


spikes = eval(['spikes.' sessionName]);

behavior = eval(['behavior.' sessionName]);


if strcmp(rat, 'Kevin') % there are two wake periods, we need to concatenate them
   behavior.time(2,2) = behavior.time(3,2);
   behavior.time(3,1) = behavior.time(4,1);
   behavior.time(3,2) = behavior.time(4,2);
   behavior.time(4,:) = [];
end

position = eval(['position.' sessionName]);
speed = eval(['speed.' sessionName]);
basics = eval(['basics.' sessionName]);
ripple = eval(['ripple.' sessionName]);


%% session info


% Fs = basics.SampleRate;
Fs = 1e6; % Since the spike timestamps are already in microsecond we don't need to use the original sampling rate. 
            %The behavioral timepoints are also in microsecond therefore easier to handle them.


% lfpSampleRate = basics.lfpSampleRate;
lfpSampleRate = 1e6; % agian we use the time unit after conversion


nCh = basics.nChannels;


fileinfo = struct('name', sessionName2, 'animal', rat, 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', Fs, 'nCh', nCh, 'lfpSampleRate', lfpSampleRate, 'pix2cm', 0.3861); 
fileinfo.pix2cm = 0.3861; % in cm


FileBase = [currDir '/' fileinfo.name];
mkdir(FileBase)



%% position info and animal wake behavior analysis


fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 


%%% calculate timing of the laps (animal travels from one end to the other end of the track)

laps = calculateLapTimings(fileinfo, speed, FileBase);


%%% see if we want to partition the data to halves (here based on the time gaps between the laps)

gapbtwLaps = laps(2:end, 1) - laps(1:end-1, 2);

figure; 
plot(2:length(laps), gapbtwLaps)

% breakLap = input('\n Enter the lap number for truncating the data\n(If no partitioning is needed just press enter) \n');  
breakLap = [];

if ~isempty(breakLap)
    firstORsecond = input('\nselect a part of data to process (one or two)\n');
else
    firstORsecond = [];
end


%% making a new spike data just to conform with the format of Kamran's data

qual2consider = 'all';

mySpikes = spikeBehaviorAnalysis(spikes, laps, ripple.time, speed, qual2consider, fileinfo, Fs);


bvrTimeList = behavior.list(:,[1 2]);
bvrState = behavior.list(:,3);
% activeTimes = bvrTimeList(bvrState == 4,:);
% noActEpochs = length(activeTimes);


%% PBEs

%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the meanposbinIdx(:,1)



% PRE

fileinfo.tbegin = behavior.time(1,1) ; 
fileinfo.tend = behavior.time(1,2);


% exclude all the periods which we are not interested in for now

exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_pre = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, mySpikes, fileinfo);

savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, 'PRE', fileinfo)

[poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE, temporalEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE);





% RUN

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [1 2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);

savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, 'RUN', fileinfo)

[poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, temporalEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);





% % POST
% 
% fileinfo.tbegin = behavior.time(3,1); 
% fileinfo.tend = behavior.time(3,2);
% 
% 
% exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4
% 
% primaryPBEs_post = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);
% 
% [eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, mySpikes, fileinfo);
% 
% savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, 'POST', fileinfo)
% 
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
% between the 


currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
mkdir(FileBase)


binDur = 0.25; % in sec 
numofStates = 40;
posBinSize = 0.5; % in cm


% active periods during run
activePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < behavior.time(2,2), :); 
% ActivePeriods = ActivePeriods(50:end, :);


% Binnig the active run period

[runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo);


activeUnits = 1:size(eventsBinnedfiringPRE{1,2}, 1);



% 
% % PRE
% pbePeriod = 'PRE';
% 
% 
% % Train an HMM on PBEs
% [transmatPRE, lambdaPRE] = trainHMM(eventsBinnedfiringPRE(:,2), activeUnits, numofStates);
% 
% 
% plot_lsPF(runBinnedfiring, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% 
% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% 
% [PosDecErrorPRE, PosDecErrorPRE_shuffledOrder, PosDecErrorPRE_shuffledGamma] = positionDecodingError(runBinnedfiring, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates);



% RUN
pbePeriod = 'RUN';

% Train an HMM on PBEs
[transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);


plot_lsPF(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);

savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')


[PosDecErrorRUN, PosDecErrorRUN_shuffledOrder, PosDecErrorRUN_shuffledGamma] = positionDecodingError(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates);


% 
% 
% % POST
% pbePeriod = 'POST';
% 
% % Train an HMM on PBEs
% [transmatPOST, lambdaPOST] = trainHMM(eventsBinnedfiringPOST(:,2), activeUnits, numofStates);
% 
% 
% plot_lsPF(runBinnedfiring, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% 
% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% 
% 
% [PosDecErrorPOST, PosDecErrorPOST_shuffledOrder, PosDecErrorPOST_shuffledGamma] = positionDecodingError(runBinnedfiring, transmatPOST, lambdaPOST, posbinIdx, xposcenters, yposcenters, numofStates);
% 
% 

% plot decoding error

figure;
set(gcf, 'Units', 'pixels', 'position', [814 399 1648 606])

subplot(1,3,1)
plotErrorCdf(PosDecErrorPRE, PosDecErrorPRE_shuffledOrder, PosDecErrorPRE_shuffledGamma, 1)
title('PRE', 'fontsize', 16)


subplot(1,3,2)
plotErrorCdf(PosDecErrorRUN, PosDecErrorRUN_shuffledOrder, PosDecErrorRUN_shuffledGamma, 0)
title('RUN', 'fontsize', 16)


subplot(1,3,3)
plotErrorCdf(PosDecErrorPOST, PosDecErrorPOST_shuffledOrder, PosDecErrorPOST_shuffledGamma, 0)
title('POST', 'fontsize', 16)

% suptitle([fileinfo.name])

savepdf(gcf, [FileBase fileinfo.name '_decodingError'])
print(gcf, [FileBase fileinfo.name  '_decodingError'], '-dpng')

save([FileBase fileinfo.name '_decodingError.mat'], 'PosDecErrorPRE', 'PosDecErrorPRE_shuffledOrder', 'PosDecErrorPRE_shuffledGamma', ...
          'PosDecErrorRUN', 'PosDecErrorRUN_shuffledOrder', 'PosDecErrorRUN_shuffledGamma', 'PosDecErrorPOST', 'PosDecErrorPOST_shuffledOrder', 'PosDecErrorPOST_shuffledGamma')



%% HMM sequence detection

% In this step we test the extent to which an HMM trained on the
% candidate events is able to detect sequence events among non-sequence
% events

% we do this by k-fold cross-validation: using k-1 folds for training and
% and the remaining fold for sequence detection %% this is temporary assigning. Later will change it when analyzing run/post-sleep data


numofFolds = 5; %% 5-fold cross-validation
numofStates = 40;
noShuffle = 2000;


% PRE

pbePeriod = 'PRE';

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPRE, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctile_pre = HMMcongruence_crossvalid(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctile_pre_p = HMMcongruence_crossvalid(poissonEventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
HMMprctile_pre_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
HMMprctile_pre_t = HMMcongruence_crossvalid(temporalEventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);



% RUN cross-validation

pbePeriod = 'RUN';

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctile_run = HMMcongruence_crossvalid(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctile_run_p = HMMcongruence_crossvalid(poissonEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
HMMprctile_run_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
HMMprctile_run_t = HMMcongruence_crossvalid(temporalEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);



% POST cross-validation

pbePeriod = 'POST';

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPOST, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctile_post = HMMcongruence_crossvalid(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctile_post_p = HMMcongruence_crossvalid(poissonEventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
HMMprctile_post_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
HMMprctile_post_t = HMMcongruence_crossvalid(temporalEventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);






%% Measuring the congruence of PBEs within a state (e.g., POST) with the model trained on a different period (e.g., RUN)

numofStates = 40;


% trained models for each part

% initialize the model parameters

noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);
prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);


[~, ~, transmatPRE, lambdaPRE] = phmm_em(eventsBinnedfiringPRE(:,2), prior0, transmat0, lambda0, 'max_iter', 200);

[~, ~, transmatRUN, lambdaRUN] = phmm_em(eventsBinnedfiringRUN(:,2), prior0, transmat0, lambda0, 'max_iter', 200);
     
[~, ~, transmatPOST, lambdaPOST] = phmm_em(eventsBinnedfiringPOST(:,2), prior0, transmat0, lambda0, 'max_iter', 200);



noShuffle = 2000;


%%% train PRE___ test RUN

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;

trainPeriod = 'PRE';
testPeriod = 'RUN';

HMMprctile_pre_run = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates

HMMprctile_pre_run_p = modelCongruence(poissonEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMprctile_pre_run_ts = modelCongruence(timeSwapEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMprctile_pre_run_t = modelCongruence(temporalEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');




%%% train RUN___ test POST

curr_transmat = transmatRUN;
curr_lambda = lambdaRUN;

trainPeriod = 'RUN';
testPeriod = 'POST'; 

HMMprctile_run_post = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates


HMMprctile_run_post_p = modelCongruence(poissonEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMprctile_run_post_ts = modelCongruence(timeSwapEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMprctile_run_post_t = modelCongruence(temporalEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');



%%% train PRE___ test POST

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;


trainPeriod = 'PRE';
testPeriod = 'POST'; 

HMMprctile_pre_post = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates

HMMprctile_pre_post_p = modelCongruence(poissonEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMprctile_pre_post_ts = modelCongruence(timeSwapEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMprctile_pre_post_t = modelCongruence(temporalEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');





figure;
set(gcf, 'position', [801 82 1648 1239])

% cross-validations (within-periods)

subplot(2,3,1)
congruencePlot(HMMprctile_pre, HMMprctile_pre_p, HMMprctile_pre_ts, HMMprctile_pre_t, 1, sprintf('Train on %s and test on %s (CV)', 'PRE', 'PRE'))

subplot(2,3,2)
congruencePlot(HMMprctile_run, HMMprctile_run_p, HMMprctile_run_ts, HMMprctile_run_t, 0, sprintf('Train on %s and test on %s (CV)', 'RUN', 'RUN'))

subplot(2,3,3)
congruencePlot(HMMprctile_post, HMMprctile_post_p, HMMprctile_post_ts, HMMprctile_post_t, 0, sprintf('Train on %s and test on %s (CV)', 'POST', 'POST'))



% congruence of PBEs with other periods' models

subplot(2,3,4)
congruencePlot(HMMprctile_pre_run, HMMprctile_pre_run_p, HMMprctile_pre_run_ts, HMMprctile_pre_run_t, 1, sprintf('Train on %s and test on %s', 'PRE', 'RUN'))

subplot(2,3,5)
congruencePlot(HMMprctile_run_post, HMMprctile_run_post_p, HMMprctile_run_post_ts, HMMprctile_run_post_t, 0, sprintf('Train on %s and test on %s', 'RUN', 'POST'))

subplot(2,3,6)
congruencePlot(HMMprctile_pre_post, HMMprctile_pre_post_p, HMMprctile_pre_post_ts, HMMprctile_pre_post_t, 0, sprintf('Train on %s and test on %s', 'PRE', 'POST'))





currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/congruence(wi_and_bw)/'];
mkdir(FileBase)


print(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
savepdf(gcf, [FileBase fileinfo.name '_congruence'])


save([FileBase fileinfo.name '_congruence.mat'], ...
    'HMMprctile_pre', 'HMMprctile_pre_p', 'HMMprctile_pre_ts', 'HMMprctile_pre_t', ...
    'HMMprctile_run', 'HMMprctile_run_p', 'HMMprctile_run_ts', 'HMMprctile_run_t', ...
    'HMMprctile_post', 'HMMprctile_post_p', 'HMMprctile_post_ts', 'HMMprctile_post_t', ...
    'HMMprctile_pre_run', 'HMMprctile_pre_run_p', 'HMMprctile_pre_run_ts', 'HMMprctile_pre_run_t', ...
    'HMMprctile_run_post', 'HMMprctile_run_post_p', 'HMMprctile_run_post_ts', 'HMMprctile_run_post_t', ...
    'HMMprctile_pre_post', 'HMMprctile_pre_post_p', 'HMMprctile_pre_post_ts', 'HMMprctile_pre_post_t')
 


%% Measure the extent to which PBEs congruent with the model trained on the 
%% same period as the PBEs (measured through cross validation) are congruent with the model trained on 
%% preceding periods 

% For example, one interesting question could be that the sequential activity during which behavioral state, 
% PRE or RUN, relates (explains) the most the significant sequences during
% POST


% Whether RUN signifant sequences are congruent with model traind of PRE
% PBEs


figure;
set(gcf, 'position', [796 125 1658 1184])


% plot 2D histogram

data = [HMMprctile_pre_run HMMprctile_run];
h = hist3(data, [10 10])/length(data);

subplot(2,2,1)

imagesc(h); set(gca, 'YDir', 'normal'); 
colormap(jet)
% colormap(flipud(colormap('hot')))
set(gca, 'fontsize', 16, 'linewidth', 3)

set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')

xlabel({'Conguence of RUN PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
ylabel({'Conguence of RUN PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)

title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)


axis square
colorbar



% plot stacked histogram

runIdx = HMMprctile_run > 0.95;
prerunIdx = HMMprctile_pre_run > 0.95;

noncong = [length(find(~runIdx & ~prerunIdx));
           length(find(~runIdx & prerunIdx))];
       
cong = [length(find(runIdx & ~prerunIdx));
        length(find(runIdx & prerunIdx))];      

       
noncong_cumcount = cumsum(noncong);
noncong_cumcount = noncong_cumcount (end:-1:1);

cong_cumcount = cumsum(cong);
cong_cumcount = cong_cumcount(end:-1:1);

thecolors =  [50 150 255 ;255 255 150]/255;


subplot(2,2,2)


hold on
for ii = 1:2
    bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
end
hold off


set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')

xlabel('Congruence with RUN model', 'fontsize', 16)
ylabel('Number of PBEs', 'fontsize', 20)

% axis square

legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
legend boxoff 






% Congruence of post PBEs with PRE and RUN models: how the two are related 

% plot 2D histogram

data = [HMMprctile_pre_post HMMprctile_run_post];
h = hist3(data, [10 10])/length(data);

subplot(2,2,3)

imagesc(h); set(gca, 'YDir', 'normal'); 
colormap(jet)
% colormap(flipud(colormap('hot')))
set(gca, 'fontsize', 16, 'linewidth', 3)

set(gca, 'XTick', [0.5 10.5], 'XTickLabel', {'0', '100'}, 'YTick', [0.5 10.5], 'YTickLabel', {'0', '100'}, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')

xlabel({'Conguence of POST PBEs'; 'with RUN model (percentile)'}, 'fontsize', 16)
ylabel({'Conguence of POST PBEs'; 'with PRE model (percentile)'}, 'fontsize', 16)

title(sprintf('n(PBEs) = %s', num2str(length(data))), 'fontsize', 14)


axis square
colorbar



% plot stacked histogram

runIdx = HMMprctile_run_post > 0.95;
prerunIdx = HMMprctile_pre_post > 0.95;

noncong = [length(find(~runIdx & ~prerunIdx));
           length(find(~runIdx & prerunIdx))];
       
cong = [length(find(runIdx & ~prerunIdx));
        length(find(runIdx & prerunIdx))];      

       
noncong_cumcount = cumsum(noncong);
noncong_cumcount = noncong_cumcount (end:-1:1);

cong_cumcount = cumsum(cong);
cong_cumcount = cong_cumcount(end:-1:1);

thecolors =  [50 150 255 ;255 255 150]/255;


subplot(2,2,4)


hold on
for ii = 1:2
    bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
end
hold off


set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'}, 'fontsize', 16)
set(gca, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')

xlabel('Congruence of POST with RUN model', 'fontsize', 16)
ylabel('Number of PBEs', 'fontsize', 20)

% axis square

legend({'cong. with PRE', 'not cong. with PRE','Location'}, 'Location', 'northoutside', 'fontsize', 14)
legend boxoff 



% % How POST significant sequences comply with the models during RUN and/or
% % PRE
% 
% 
% % Plot stacked histogram
% subplot(1,3,3)
% 
% postIdx = HMMprctile_post > 0.95;
% prepostIdx = HMMprctile_pre_post > 0.95;
% runpostIdx = HMMprctile_run_post > 0.95;
% 
% noncong = [length(find(~postIdx & ~prepostIdx & ~runpostIdx)); ...
%            length(find(~postIdx & prepostIdx & ~runpostIdx)); ...
%            length(find(~postIdx & ~prepostIdx & runpostIdx));...
%            length(find(~postIdx & prepostIdx & runpostIdx))];
%        
% cong = [length(find(postIdx & ~prepostIdx & ~runpostIdx)); ...
%            length(find(postIdx & prepostIdx & ~runpostIdx)); ...
%            length(find(postIdx & ~prepostIdx & runpostIdx));...
%            length(find(postIdx & prepostIdx & runpostIdx))];
%        
%        
% noncong_cumcount = cumsum(noncong);
% noncong_cumcount = noncong_cumcount (end:-1:1);
% 
% cong_cumcount = cumsum(cong);
% cong_cumcount = cong_cumcount(end:-1:1);
% 
% thecolors = [100 255 180; 255 100 100; 50 150 255 ;255 255 150]/255;
% 
% 
% hold on
% for ii = 1:4
%     bar(1:2, [noncong_cumcount(ii) cong_cumcount(ii)], 'FaceColor', thecolors(ii, :), 'EdgeColor', 'none')
% end
% 
% 
% set(gca, 'xlim', [0.5 2.5],'XTick', [1 2], 'XTickLabel',  {'non-cong.', 'congruent'})
% set(gca, 'fontsize', 16, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'Top')
% 
% xlabel('Congruence with Post BPE model', 'fontsize', 16)
% ylabel('Number of PBEs', 'fontsize', 20)
% 
% % axis square
% 
% legend({'both pre and run', 'run excl.', 'pre excl.', 'neither pre nor run'}, 'Location', 'northoutside','FontSize',14)
% legend boxoff 



FileBase = [currDir   '/HMM/sequenceDetection/overlap/'];
mkdir(FileBase)

print(gcf, [FileBase  '_' 'Overlap_of_wibw'], '-dpng')
savepdf(gcf, [FileBase  '_' 'Overlap_of_wibw'])


end

end





%% Functions

function [eventsBinnedfiring, secondaryPBEs, idx_acceptedEvents] = finalBinningResult(pbEvents2, mySpikes, fileinfo)

%% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events


%%% Calculating the number of firing pyramidal units within each event

qclus = [1 2 3]; % pyramidal and interneurons

noFiringUnits = zeros(size(pbEvents2, 1), 1);
for evt = 1 : size(pbEvents2, 1)
    
    evtSpkInd = find(mySpikes.t > pbEvents2(evt, 1) & mySpikes.t < pbEvents2(evt, 2) & ismember(mySpikes.qclu, qclus));
    noFiringUnits(evt) = length(unique(mySpikes.unit(evtSpkInd)));
end


%%% binning the spikes within each event

binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) 
eventsBinnedfiring2 = timeBinning(pbEvents2, mySpikes, qclus, binDur, fileinfo.Fs);


% Remove the flanking zero-firing periods(silent bins) from the beginnig
% and end of each event. Since the interneurons were involved in calculating
% the PBEs' boundaries and here we are binning just pyramidal units' spike
% trains, some silent bins are expected. 

[eventsBinnedfiring, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring2, pbEvents2(:, 1), binDur, fileinfo.Fs);

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



function HMMprctile = modelCongruence(testData, transmat, lambda, noShuffle, fileinfo, trainPeriod, testPeriod, surrogateData)



noEvents = size(testData, 1);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, noShuffle); %% log likelihood of shuffle data

HMMprctile = zeros(noEvents, 1);


    setaside = [];
    shuffleTransmats = genshuffles(transmat, noShuffle, setaside);




for evt = 1 : noEvents
    
    currEvent = testData{evt, 2};
%     noTimebins = size(currEvent, 2);

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
    
    if mod(evt, 20) == 0
        fprintf(1, [trainPeriod '-' testPeriod '-' surrogateData '__event %d, HMMPercentile = %f\n'], evt, HMMprctile(evt, 1));
    end
    
    
end



currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/bw_periods_congruence/'];
mkdir(FileBase)


save([FileBase '/Train_' trainPeriod '___Test_' testPeriod '__' surrogateData '.mat'], 'HMMprctile', 'dataLL', 'nullLL')

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

periods = activePeriods * 1/fileinfo.lfpSampleRate;
runBinPos = cell(size(periods, 1), 1);

for ii = 1 : size(periods, 1)
    
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/fileinfo.lfpSampleRate, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin
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


function [PosDecError, PosDecError_shuffledOrder, PosDecError_shuffledGamma] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates)


noTimeBins = length(posbinIdx);

% 5-fold cross-validation

noFolds = 5;
foldSize = floor(noTimeBins/ noFolds);

decodedPositions_actual = zeros(length(posbinIdx), 2);
decodedPositions_shuffledOrder = zeros(length(posbinIdx), 2);
decodedPositions_shuffledGamma = zeros(length(posbinIdx), 2);


for k = 1:noFolds
    
    
    if k == noFolds
        testSet = (k-1)*foldSize+1 : noTimeBins;
    else
        testSet = (k-1)*foldSize+1 : k*foldSize;
    end
    
    trainSet = setdiff(1:noTimeBins, testSet);
    
    
    %% calculate lfPFs for the train set
    
    
    train_timeBins = runBinnedfiring{2}(:,trainSet);
    train_posbinIdx = posbinIdx(trainSet, :); % x and y coordiantes as columns
    
    train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_posbinIdx, [], numofStates, xposcenters, yposcenters);

    
    % Decode the states in RUN data 
    
    test_timeBins = runBinnedfiring{2}(:,testSet);

    obsProb = poisson_prob(test_timeBins, lambda,1);
      
    prior = 1/numofStates * ones(numofStates,1); % a uniform prior
    %(1) actual time bin orders
    
    decodedPositions_actual(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'actual'); 
    
    
    %(2) shuffled time bin orders
    
    decodedPositions_shuffledOrder(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledOrder'); 
    
    %(3) shuffled gamma: it supposed to take away all the order and
    %coactivity information regarding the states
    
    decodedPositions_shuffledGamma(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledGamma'); 
    
    
    

        
end

posbinIdx(find(posbinIdx(:,2) > length(yposcenters)), 2) = length(yposcenters);
posbinIdx(find(posbinIdx(:,1) > length(xposcenters)), 1) = length(xposcenters);


PosDecError = sqrt((decodedPositions_actual(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_actual(:, 2) - xposcenters(posbinIdx(:,1))').^2);

PosDecError_shuffledOrder = sqrt((decodedPositions_shuffledOrder(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledOrder(:, 2) - xposcenters(posbinIdx(:,1))').^2);

PosDecError_shuffledGamma = sqrt((decodedPositions_shuffledGamma(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledGamma(:, 2) - xposcenters(posbinIdx(:,1))').^2);    
    
end


function decodedPositions = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, shuffleMode)



noTestBins = size(obsProb, 2);
numofStates = size(transmat, 1);
shuffledOrder = randperm(noTestBins, noTestBins);

    
if strcmp(shuffleMode, 'shuffledOrder')
    [~, ~, gamma(:, shuffledOrder)] = fwdback(prior, transmat, obsProb(:, shuffledOrder));
else
    [~, ~, gamma] = fwdback(prior, transmat, obsProb);
end

decodedPositions = zeros(noTestBins, 2);
for jj = 1: noTestBins
    
        
    if strcmp(shuffleMode, 'shuffledGamma')
        stateProbability = gamma(randperm(numofStates, numofStates), jj);
    else
        stateProbability = gamma(:, jj);
    end
    
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


function plotErrorCdf(data, shuffledOrder, shuffledGamma, ylabelneeded)

hold on

% actual data
[h1, bins] = hist(data, 100);
h1c = cumsum(h1)/sum(h1);

dataMed = median(data);
 
curve1 = plot(bins, h1c, 'linewidth', 4, 'color', [150, 150, 255]/255);
line([dataMed dataMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
line([0 dataMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
text(dataMed, 0.03, sprintf('%.1f cm', dataMed), 'color', [150, 150, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')

% shuffled order
[h2, bins] = hist(shuffledOrder, 100);
h2c = cumsum(h2)/sum(h2);

shuffleMed = median(shuffledOrder);

curve2 = plot(bins, h2c, 'linewidth', 4, 'color', [150, 200, 255]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
text(shuffleMed, 0.07, sprintf('%.1f cm', shuffleMed), 'color', [150, 200, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')


% shuffled Gamma
[h3, bins] = hist(shuffledGamma, 100);
h3c = cumsum(h3)/sum(h3);

shuffleMed = median(shuffledGamma);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(shuffleMed, 0.1, sprintf('%.1f cm', shuffleMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')




xlim([0 max([data; shuffledOrder; shuffledGamma])])
set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])


xlabel('Decoding Error(cm)', 'fontsize', 16)

if ylabelneeded
    ylabel('Cumulative ratio', 'fontsize', 24)
    
%     legend('actual', 'shuffle state probability', 'shuffle position index', 'Location', 'northwest'); 
    legend([curve1, curve2, curve3], 'actual', 'shuffled run bins', 'shuffled gamma dist.', 'Location', 'northwest'); 

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
    
    % if doing column-wise shuffle
    
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


% 
% function PosDecError = positionDecodingError2(gamma, posbinIdx, xposcenters, yposcenters, mode)
% 
% noTimeBins = length(posbinIdx);
% numofStates = size(gamma, 1);
% 
% noxPosBins = length(xposcenters);
% noyPosBins = length(yposcenters);
% 
% % 5-fold cross-validation
% 
% noFolds = 5;
% foldSize = floor(noTimeBins/ noFolds);
% decodedPositions = zeros(length(posbinIdx), 2);
% 
% for k = 1:noFolds
%     
%     
%     if k == noFolds
%         testSet = (k-1)*foldSize+1 : noTimeBins;
%     else
%         testSet = (k-1)*foldSize+1 : k*foldSize;
%     end
%     
%     trainSet = setdiff(1:noTimeBins, testSet);
%     
%     
%     % calculate lfPFs for the train set
%     
%     train_gamma = gamma(:, trainSet);
%     train_posbinIdx = posbinIdx(trainSet, :); % x and y coordiantes as columns
% 
% 
%     gamma_avg = zeros(noyPosBins, noxPosBins, numofStates); % doing a summation over the state probabilty distributions corresponding to each position bin
% 
% 
%     for ii  = 1 : noyPosBins
%         for jj = 1: noxPosBins
% 
%             idx = find(train_posbinIdx(:,2) == ii & train_posbinIdx(:,1) == jj);
% 
%             if ~isempty(idx)
%                 gamma_avg(ii, jj, :) = permute(sum(train_gamma(:,idx), 2)/length(idx), [3,2,1]);
%             end
% 
%         end
%     end
% 
%     
%     sigma = 3; % the same smoothing parameters used for place fields calculation of units 
%     halfwidth = 3 * sigma;
%     smoothwin = gausswindow(sigma, halfwidth);
%     smoothwin2D = smoothwin' * smoothwin;
% 
% 
%     for ii = 1 : numofStates
%         gamma_avg(:, :, ii) = conv2(gamma_avg(:, :, ii), smoothwin2D, 'same');
%     end
% 
%     train_lsPFs = gamma_avg ./ repmat(sum(sum(gamma_avg, 1),2), [size(gamma_avg, 1) size(gamma_avg, 2)]); % normalize for each state
%     
%     % Decode the position in the test time bins
%     
% %     test_gamma = gamma(:, testSet);
% %     test_posbinIdx = posbinIdx(testSet, :); % x and y coordiantes as columns
%     
%     for jj = 1:length(testSet)
%         
%             
%         if strcmp(mode, 'shuffleState')
%             stateProbability = gamma(randperm(numofStates, numofStates), testSet(jj));
%         else
%             stateProbability = gamma(:, testSet(jj));
%         end
%     
%         probabilityOverPosition = sum(repmat(permute(stateProbability, [3 2 1]), length(yposcenters), length(xposcenters)) .* train_lsPFs, 3);
% 
%         decodedPositions(testSet(jj), 1) = yposcenters(ceil(sum((1:length(yposcenters))' .* sum(probabilityOverPosition, 2))));
%         decodedPositions(testSet(jj), 2) = xposcenters(ceil(sum((1:length(xposcenters)) .* sum(probabilityOverPosition, 1))));
%     end
%     
% end
% 
% 
% % if strcmp(mode, 'shufflePosition')
% %     posbinIdx = posbinIdx(randperm(length(posbinIdx)), :);
% % end
% 
% 
% try
% 
%     PosDecError = sqrt((decodedPositions(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions(:, 2) - xposcenters(posbinIdx(:,1))').^2);
% 
% catch
%     
%     posbinIdx(find(posbinIdx(:,2) > length(yposcenters)), 2) = length(yposcenters);
% 
%     posbinIdx(find(posbinIdx(:,1) > length(xposcenters)), 1) = length(xposcenters);
% 
%     PosDecError = sqrt((decodedPositions(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions(:, 2) - xposcenters(posbinIdx(:,1))').^2);
% 
%     
% end
% 
% 
% 
% 
% end
% 

