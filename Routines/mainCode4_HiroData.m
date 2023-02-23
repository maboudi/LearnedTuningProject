clear; 
clc; 
% close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3], [1 2 3], 1};


for rr = 1 %:length(rats) % loop for all the rats and sessions

    
rat = rats{rr};

sessionNumbers = allsessionNumbers{rr};


for sessionNumber = 1 % sessionNumbers


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
Fs = 1e6; % Since in Hiro's data the spike timestamps are already in microsecond we don't need to use the original sampling rate. 
            %The behavioral timepoints are also in microsecond therefore easier to handle them.

% lfpSampleRate = basics.lfpSampleRate;
lfpSampleRate = 1e6; % agian we use the time unit after conversion


nCh = basics.nChannels;


fileinfo = struct('name', sessionName2, 'animal', rat, 'xyt', [], ...
    'tbegin', [], 'tend', [], 'Fs', Fs, 'nCh', nCh, 'lfpSampleRate', lfpSampleRate, 'pix2cm', 0.3861); 
fileinfo.pix2cm = 0.3861; % in cm


FileBase = fullfile(currDir, fileinfo.name);
mkdir(FileBase)




%% position info and animal wake behavior analysis


fileinfo.xyt = [position.x*fileinfo.pix2cm; position.y*fileinfo.pix2cm; position.t]'; 

%%% calculate timing of the laps (animal travels from one end to the other end of the track)

laps = calculateLapTimings(fileinfo, speed, FileBase);


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


%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

mySpikes = spikeBehaviorAnalysis(spikes, laps, ripple.time, speed, qual2consider, fileinfo, Fs);


bvrTimeList = behavior.list(:,[1 2]);
bvrState = behavior.list(:,3);

% activeTimes = bvrTimeList(bvrState == 4,:);
% noActEpochs = length(activeTimes);


%% PBEs

%% Detemining population burst periods

% For detection of the population bursts, Kamran suggested using all the cluster
% includeing the noise. Therefore, it's better to load the spiketimes from
% res files. I don't have access to the res files for now, so just keep up
% with the sorted spikes


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the meanposbinIdx(:,1)


%% PRE


fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(1,2);

% exclude all the periods which we are not interested in for now

exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

[PBEs, sdfZ] = DetectPopBursts(mySpikes, time_resolution, threshZ, exclude, fileinfo);


primaryPBEs_pre = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, mySpikes, fileinfo);

savePBEs(primaryPBEs_pre, secondaryPBEs_pre, acceptedprmEvts_pre, 'PRE', fileinfo)

% [poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE, temporalEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE);



%% RUN


fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);

savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, 'RUN', fileinfo)

% [poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, temporalEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);



%% POST


fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);


exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_post = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, mySpikes, fileinfo);

savePBEs(primaryPBEs_post, secondaryPBEs_post, acceptedprmEvts_post, 'POST', fileinfo)

% [poissonEventsBinnedfiringPOST, timeSwapEventsBinnedfiringPOST, temporalEventsBinnedfiringPOST] = genSurrogates(eventsBinnedfiringPOST);



%% Example Hidden markov Models trained separately on PRE, RUN, and POST periods


% PRE
numofStates = 40;
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




%% Correspondence between the states from models trained on different periods(PRE, RUN, POST) 

% The idea is to train models separately on different periods,
% then decode the probability over states in the test period(s)

nPBEs_PRE = length(eventsBinnedfiringPRE(:,2));
nPBEs_RUN = length(eventsBinnedfiringRUN(:,2));
nPBEs_POST = length(eventsBinnedfiringPOST(:,2));


PBES_PRE_RUN = [eventsBinnedfiringPRE(:,2); eventsBinnedfiringRUN(:,2)];
PBES_PRE_POST = [eventsBinnedfiringPRE(:,2); eventsBinnedfiringPOST(:,2)];
PBES_RUN_POST = [eventsBinnedfiringRUN(:,2); eventsBinnedfiringPOST(:,2)];


% generate pooled time swap PBEs and the concatenated data

eventsBinnedfiringPRE_ts = genTimeSwap(eventsBinnedfiringPRE(:,2));
eventsBinnedfiringRUN_ts = genTimeSwap(eventsBinnedfiringRUN(:,2));
eventsBinnedfiringPOST_ts = genTimeSwap(eventsBinnedfiringPOST(:,2));

PBES_PRE_RUN_ts = [eventsBinnedfiringPRE_ts; eventsBinnedfiringRUN_ts];
PBES_PRE_POST_ts = [eventsBinnedfiringPRE_ts; eventsBinnedfiringPOST_ts];
PBES_RUN_POST_ts = [eventsBinnedfiringRUN_ts; eventsBinnedfiringPOST_ts];




%% find the correspondences between pairs of states 

% (1) inner products of the population vectors (from observation matrix) 


% PRE and RUN
[RUN_PRE_PopVecSim, RUN_PRE_lambdaSortIdx] = PopVecSimilarity(lambdaRUN, lambdaPRE);
    
% RUN and POST
[POST_RUN_PopVecSim, POST_RUN_lambdaSortIdx] = PopVecSimilarity(lambdaPOST, lambdaRUN);

% PRE and POST
[POST_PRE_PopVecSim, POST_PRE_lambdaSortIdx] = PopVecSimilarity(lambdaPOST, lambdaPRE);


% (2_A) By decoding the states within a test period (concatenating the
% periods the models were trained on) and calculate for each pair of
% states their joint decoding probability within each time bin and
% averaging over all of the bins


uniformTransMat = 0;
seqContentCorrespondence = 2; 

% if 0 -> the correspondences are mostly based on coactivities 

% if 1 -> controling for coactivities to what extent the sequence
% correspond to each other just based on the sequential contents. Pooled
% time swap is used for generating null distributions. This method does not
% control for the self-transitions.

% if 2 -> correspondences based on sequential contents. in this case,
% transition matrix shuffle (self-transitions spared) is used for the null
% distribution.


% PRE and RUN
[RUN_PRE, sortRUN_PRE] = stateCorrespondence(PBES_PRE_RUN, nPBEs_PRE,  lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);

% RUN and POST
[POST_RUN, sortPOST_RUN] = stateCorrespondence(PBES_RUN_POST, nPBEs_RUN, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, seqContentCorrespondence);

% PRE and POST
[POST_PRE, sortPOST_PRE] = stateCorrespondence(PBES_PRE_POST, nPBEs_PRE, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);


uniformTransMat = 0;
seqContentCorrespondence = 2; % In this section we are not focusing on the sequential content

% PRE and RUN
[RUN_PRE_ts] = stateCorrespondence(PBES_PRE_RUN_ts, nPBEs_PRE, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);

% RUN and POST
[POST_RUN_ts] = stateCorrespondence(PBES_RUN_POST_ts,  nPBEs_RUN, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, seqContentCorrespondence);

% PRE and POST
[POST_PRE_ts] = stateCorrespondence(PBES_PRE_POST_ts, nPBEs_PRE, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);



% 
% % (2_B) The extent to which the correspondences remain in the third period:
% % e.g., how the correspondence between the RUN and POST period (replay) is
% % reflected in PRE 
% 
% uniformTransMat = 0;
% seqContentCorrespondence = 1;
% 
% % RUN_PRE correspondence within POST
% [RUN_PRE_withinPOST, sortRUN_PRE_withinPOST] = stateCorrespondence(cntPBES_POST, lambdaRUN, transmatRUN, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% RUN_PRE_withinPOST_subset = RUN_PRE_withinPOST(RUN_PRE > -log10(0.05)); % /numofStates/numofStates
% 
% % POST_RUN correspondence within POST
% [POST_RUN_withinPRE, sortPOST_RUN_withinPRE] = stateCorrespondence(cntPBES_PRE, lambdaPOST, transmatPOST, lambdaRUN, transmatRUN, uniformTransMat, seqContentCorrespondence);
% POST_RUN_withinPRE_subset = POST_RUN_withinPRE(POST_RUN > -log10(0.05));
% 
% 
% % POST_PRE correspondence within POST
% [POST_PRE_withinRUN, sortPOST_PRE_withinRUN] = stateCorrespondence(cntPBES_RUN, lambdaPOST, transmatPOST, lambdaPRE, transmatPRE, uniformTransMat, seqContentCorrespondence);
% POST_PRE_withinRUN_subset = POST_PRE_withinRUN(POST_PRE > -log10(0.05));
% 

%%
% plot the correspondence matrices
figure;
set(gcf, 'position', [-1405 166 665 979])
%

currTitle = {'Population vectors inner products';''};

allValues = [RUN_PRE_PopVecSim(:); POST_RUN_PopVecSim(:); POST_PRE_PopVecSim(:)];
ccodeRange = [prctile(allValues, 0.1) prctile(allValues, 99.9)];

subplot(3,3,1)
plotPanel(RUN_PRE_PopVecSim, sortRUN_PRE, 'RUN', 'PRE', ccodeRange, currTitle)

subplot(3,3,2)
plotPanel(POST_RUN_PopVecSim, sortPOST_RUN, 'POST', 'RUN', ccodeRange, '')

subplot(3,3,3)
plotPanel(POST_PRE_PopVecSim, sortPOST_PRE, 'POST', 'PRE', ccodeRange, '')


%

currTitle = {'Correspondence between states';''};

allValues = [RUN_PRE(:); POST_RUN(:); POST_PRE(:); RUN_PRE_ts(:); POST_RUN_ts(:); POST_PRE_ts(:)];
ccodeRange = [prctile(allValues, 0.1) prctile(allValues, 99.9)];

subplot(3,3,4)
plotPanel(RUN_PRE, sortRUN_PRE, 'RUN', 'PRE', ccodeRange, currTitle)

subplot(3,3,5)
plotPanel(POST_RUN, sortPOST_RUN, 'POST', 'RUN', ccodeRange, '')

subplot(3,3,6)
plotPanel(POST_PRE, sortPOST_PRE, 'POST', 'PRE', ccodeRange, '')





currTitle = {'Correspondence between states time swap';''};

subplot(3,3,7)
plotPanel(RUN_PRE_ts, sortRUN_PRE, 'RUN', 'PRE', ccodeRange, currTitle)

subplot(3,3,8)
plotPanel(POST_RUN_ts, sortPOST_RUN, 'POST', 'RUN', ccodeRange, '')

subplot(3,3,9)
plotPanel(POST_PRE_ts, sortPOST_PRE, 'POST', 'PRE', ccodeRange, '')


% 
% 
% currTitle = {'Correspondence w/in the third period';''};
% 
% subplot(4,3,7)
% plotPanel(RUN_PRE_withinPOST, sortRUN_PRE, 'RUN', 'PRE', [], currTitle)
% 
% subplot(4,3,8)
% plotPanel(POST_RUN_withinPRE, sortPOST_RUN, 'POST', 'RUN', [], '')
% 
% subplot(4,3,9)
% plotPanel(POST_PRE_withinRUN, sortPOST_PRE, 'POST', 'PRE', [], '')
% 
% 
% %
% 
% currTitle = 'Histogram of correspondences';
% 
% subplot(4,3,10)
% plotHistCorrespondence(RUN_PRE_withinPOST_subset(:), RUN_PRE_withinPOST(:), 'RUN', 'PRE', currTitle)
% 
% subplot(4,3,11)
% plotHistCorrespondence(POST_RUN_withinPRE_subset(:), POST_RUN_withinPRE(:), 'POST', 'RUN', '')
% 
% subplot(4,3,12)
% plotHistCorrespondence(POST_PRE_withinRUN_subset(:), POST_PRE_withinRUN(:), 'POST', 'PRE', '')

% 

% FileBase = [currDir '/' fileinfo.name  '/HMM/StateCorrespondence/'];
% mkdir(FileBase)
% 
% savepdf(gcf, fullfile(FileBase, [fileinfo.name '_stateCorrespondence_seqContent_uniformTransmat'])); % _seqContent_uniformTransmat

%% 


[POSTstates_run, RUNstates] = find(POST_RUN > -log10(0.05/40/40));
[POSTstates_pre, PREstates_AlsoinPOST] = find(POST_PRE > -log10(0.05/40/40));
[RUNstates_pre, PREstates_AlsoinRUN] = find(RUN_PRE > -log10(0.05/40/40));


POST_RUN_signifCorresps = POST_RUN(POST_RUN > -log10(0.05/40/40));
POST_PRE_signifCorresps = POST_PRE(POST_PRE > -log10(0.05/40/40));
RUN_PRE_signifCorresps = RUN_PRE(RUN_PRE > -log10(0.05/40/40));


noRandomization = 5000;


POST_RUN_sim = Lambdasimilarity(lambdaPOST(:, POSTstates_run), lambdaRUN(:, RUNstates));
POST_RUN_simChance = Lambdasimilarity(lambdaPOST(:, randi(40, 1, noRandomization)), lambdaRUN(:, randi(40, 1, noRandomization)));
% POST_RUN_all = Lambdasimilarity(lambdaPOST, lambdaRUN);


POST_PRE_sim = Lambdasimilarity(lambdaPOST(:, POSTstates_pre), lambdaPRE(:, PREstates_AlsoinPOST));
POST_PRE_simChance = Lambdasimilarity(lambdaPOST(:, randi(40, 1, noRandomization)), lambdaPRE(:, randi(40, 1, noRandomization)));
% POST_PRE_all = Lambdasimilarity(lambdaPOST, lambdaPRE);


RUN_PRE_sim = Lambdasimilarity(lambdaRUN(:, RUNstates_pre), lambdaPRE(:, PREstates_AlsoinRUN));
RUN_PRE_simChance = Lambdasimilarity(lambdaRUN(:, randi(40, 1, noRandomization)), lambdaPRE(:, randi(40, 1, noRandomization)));
% RUN_PRE_all = Lambdasimilarity(lambdaRUN, lambdaPRE);



bins1 = linspace(min([RUN_PRE_sim; RUN_PRE_simChance]), max([RUN_PRE_sim; RUN_PRE_simChance]), 50);
h1 = hist(RUN_PRE_sim, bins1)/length(RUN_PRE_sim);
h1_chance = hist(RUN_PRE_simChance, bins1)/length(RUN_PRE_simChance);
pval1 = ranksum(RUN_PRE_sim, RUN_PRE_simChance, 'tail', 'right');


bins2 = linspace(min([POST_RUN_sim; POST_RUN_simChance]), max([POST_RUN_sim; POST_RUN_simChance]), 50);
h2 = hist(POST_RUN_sim, bins2)/length(POST_RUN_sim);
h2_chance = hist(POST_RUN_simChance, bins1)/length(POST_RUN_simChance);
pval2 = ranksum(POST_RUN_sim, POST_RUN_simChance, 'tail', 'right');


bins3 = linspace(min([POST_PRE_sim; POST_PRE_simChance]), max([POST_PRE_sim; POST_PRE_simChance]), 50);
h3 = hist(POST_PRE_sim, bins3)/length(POST_PRE_sim);
h3_chance = hist(POST_PRE_simChance, bins3)/length(POST_PRE_simChance);
pval3 = ranksum(POST_PRE_sim, POST_PRE_simChance, 'tail', 'right');





figure;
set(gcf, 'position', [-2129 886 885 320])

subplot(1,3,1)
set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.005, 0.01])
hold on

bar(bins1, h1_chance, 'FaceColor', 'none', 'EdgeColor', 'k')
bar(bins1, h1, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)


h1_chance_sm = conv(h1_chance, gausswindow(2,5), 'same');
h1_sm = conv(h1, gausswindow(2,5), 'same');
plot(bins1, h1_chance_sm, 'color', [0.7 0.7 0.7], 'linewidth', 3)
plot(bins1, h1_sm, 'k', 'linewidth', 3)


ylim=get(gca,'ylim');
xlim=get(gca,'xlim');

if pval1 < 0.001
    text(mean(xlim),mean(ylim), 'p < 0.001', 'fontsize', 14)
else
    text(mean(xlim),mean(ylim), sprintf('p = %.2f', pval1), 'fontsize', 14)
end

xlabel('Pop. vector product', 'fontsize', 14)
ylabel('probability', 'fontsize', 14)
title('PRE and RUN', 'fontsize', 14)





subplot(1,3,2)
set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.005, 0.01])

hold on

bar(bins2, h2_chance, 'FaceColor', 'none', 'EdgeColor', 'k') 
bar(bins2, h2, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)


h2_chance_sm = conv(h2_chance, gausswindow(2,5), 'same');
h2_sm = conv(h2, gausswindow(2,5), 'same');
plot(bins2, h2_chance_sm, 'color', [0.7 0.7 0.7], 'linewidth', 3)
plot(bins2, h2_sm, 'k', 'linewidth', 3)


ylim=get(gca,'ylim');
xlim=get(gca,'xlim');

if pval2 < 0.001
    text(mean(xlim),mean(ylim), 'p < 0.001', 'fontsize', 14)
else
    text(mean(xlim),mean(ylim), sprintf('p = %.2f', pval2), 'fontsize', 14)
end

xlabel('Pop. vector product', 'fontsize', 14)
title({'population vector similarity between correspondent states';'';'POST and RUN'}, 'fontsize', 14)




subplot(1,3,3)
set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.005, 0.01])

hold on

bar(bins3, h3_chance, 'FaceColor', 'none', 'EdgeColor', 'k')
bar(bins3, h3, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)

h3_chance_sm = conv(h3_chance, gausswindow(2,5), 'same');
h3_sm = conv(h3, gausswindow(2,5), 'same');
pl2 = plot(bins3, h3_chance_sm, 'color', [0.7 0.7 0.7], 'linewidth', 3);
pl1 = plot(bins3, h3_sm, 'k', 'linewidth', 3);


ylim=get(gca,'ylim');
xlim=get(gca,'xlim');

if pval3 < 0.001
    text(mean(xlim),mean(ylim), 'p < 0.001', 'fontsize', 14)
else
    text(mean(xlim),mean(ylim), sprintf('p = %.2f', pval3), 'fontsize', 14)
end

xlabel('Pop. vector product', 'fontsize', 14)
title('POST and PRE', 'fontsize', 14)

legend([pl1, pl2], 'corresp. states', 'random states')

% savepdf(gcf, [fileinfo.name '_PopVecSimilarity'])



%% HMM sequence detection

% In this step we test the extent to which an HMM trained on the
% PBEs is able to detect the PBEs with significant sequential activitiy
% patterns


% (1) Cross-validation: the model trained on a period is tested on the same period 

numofFolds = 5; %% 5-fold cross-validation
numofStates = 40;
noShuffle = 500;


%% PRE cross-valodation

pbePeriod = 'PRE';

eventsBinnedfiringPRE = eventsBinnedfiringPRE(randperm(size(eventsBinnedfiringPRE, 1)), :); % randomizing the temporal order of events across the period


[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPRE, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctilePRE = HMMcongruence_crossvalid(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctilePRE_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod); % pts: pooled time swap



%% RUN cross-validation

pbePeriod = 'RUN';

eventsBinnedfiringRUN = eventsBinnedfiringRUN(randperm(size(eventsBinnedfiringRUN, 1)), :);


[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctileRUN = HMMcongruence_crossvalid(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctileRUN_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod); 



%% POST cross-validation

pbePeriod = 'POST';

eventsBinnedfiringPOST = eventsBinnedfiringPOST(randperm(size(eventsBinnedfiringPOST, 1)), :);


[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPOST, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctilePOST = HMMcongruence_crossvalid(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctilePOST_pts = HMMcongruence_crossvalid_pts(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);



%% Measuring the congruence of PBEs from a period (e.g., POST) with the model trained on a different period (e.g., RUN)


noShuffle = 500;

%%% train PRE___ test RUN

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;



trainPeriod = 'PRE';
testPeriod = 'RUN';

HMMprctile_PRE_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_PRE_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);



%%% train RUN___ test POST

curr_transmat = transmatRUN;
curr_lambda = lambdaRUN;



trainPeriod = 'RUN';
testPeriod = 'POST'; 

HMMprctile_RUN_POST = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_RUN_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);


%%% train PRE___ test POST

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;



trainPeriod = 'PRE';
testPeriod = 'POST'; 

HMMprctile_PRE_POST = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_PRE_POST_pts = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);


%%% train RUN___ test PRE

curr_transmat = transmatRUN;
curr_lambda = lambdaRUN;



trainPeriod = 'RUN';
testPeriod = 'PRE'; 

HMMprctile_RUN_PRE = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_RUN_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);



%%% train POST___ test PRE

curr_transmat = transmatPOST;
curr_lambda = lambdaPOST;



trainPeriod = 'POST';
testPeriod = 'PRE'; 

HMMprctile_POST_PRE = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_POST_PRE_pts = modelCongruence(eventsBinnedfiringPRE, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);


%%% train POST___ test RUN

curr_transmat = transmatPOST;
curr_lambda = lambdaPOST;



trainPeriod = 'POST';
testPeriod = 'RUN'; 

HMMprctile_POST_RUN = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 0);

% _____________surrogates
HMMprctile_POST_RUN_pts = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 1);




figure;
set(gcf, 'position', [129 64 2408 1232])


% Cross-validations 

subplot(2,6,1)
congruencePlot2(HMMprctilePRE, HMMprctilePRE_pts, 1, sprintf('%s Cross-validation (N = %d)', 'PRE', length(HMMprctilePRE)))

subplot(2,6,2)
congruencePlot2(HMMprctileRUN, HMMprctileRUN_pts, 1, sprintf('%s Cross-validation (N = %d)', 'RUN', length(HMMprctileRUN)))

subplot(2,6,3)
congruencePlot2(HMMprctilePOST, HMMprctilePOST_pts, 1, sprintf('%s Cross-validation (N = %d)', 'POST', length(HMMprctileRUN)))



% Congruence of PBEs with other periods' models

subplot(2,6,7)
congruencePlot2(HMMprctile_PRE_RUN, HMMprctile_PRE_RUN_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'RUN'))

subplot(2,6,8)
congruencePlot2(HMMprctile_RUN_POST, HMMprctile_RUN_POST_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'POST'))

subplot(2,6,9)
congruencePlot2(HMMprctile_PRE_POST, HMMprctile_PRE_POST_pts, 1, sprintf('Train on %s and test on %s', 'PRE', 'POST'))

subplot(2,6,10)
congruencePlot2(HMMprctile_RUN_PRE, HMMprctile_RUN_PRE_pts, 1, sprintf('Train on %s and test on %s', 'RUN', 'PRE'))

subplot(2,6,11)
congruencePlot2(HMMprctile_POST_PRE, HMMprctile_POST_PRE_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'PRE'))

subplot(2,6,12)
congruencePlot2(HMMprctile_POST_RUN, HMMprctile_POST_RUN_pts, 1, sprintf('Train on %s and test on %s', 'POST', 'RUN'))



currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/congruence/'];
mkdir(FileBase)


print(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
savepdf(gcf, [FileBase fileinfo.name '_congruence'])


save([FileBase fileinfo.name '_congruence.mat'], ...
    'HMMprctilePRE', 'HMMprctilePRE_pts', ...
    'HMMprctileRUN', 'HMMprctileRUN_pts', ...
    'HMMprctilePOST', 'HMMprctilePOST_pts', ...
    'HMMprctile_PRE_RUN', 'HMMprctile_PRE_RUN_pts', ...
    'HMMprctile_RUN_POST', 'HMMprctile_RUN_POST_pts', ...
    'HMMprctile_PRE_POST', 'HMMprctile_PRE_POST_pts', ...
    'HMMprctile_RUN_PRE', 'HMMprctile_RUN_PRE_pts', ...
    'HMMprctile_POST_PRE', 'HMMprctile_POST_PRE_pts', ...
    'HMMprctile_POST_RUN', 'HMMprctile_POST_RUN_pts')
 


%% 2D histograms of sequence scores 

% example: likelihood(POST|PRE) vs likelihood(POST|RUN)

% The extent to which RUN signifant PBEs are congruent with model traind of PRE

figure;
set(gcf, 'position', [796 125 1658 1184])


% plot 2D histogram

data = [HMMprctile_PRE_RUN HMMprctileRUN];
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

runIdx = HMMprctileRUN > 0.95;
prerunIdx = HMMprctile_PRE_RUN > 0.95;

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

data = [HMMprctile_PRE_POST HMMprctile_RUN_POST];
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

runIdx = HMMprctile_RUN_POST > 0.95;
prerunIdx = HMMprctile_PRE_POST > 0.95;

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


FileBase = [currDir  '/HMM/sequenceDetection/overlap/'];
mkdir(FileBase)

print(gcf, [FileBase  '_' '2Dhistograms'], '-dpng')
savepdf(gcf, [FileBase  '_' '2Dhistograms'])




% %% latent state place fields (lsPFs); Decoding of position using lsPFs and decoding error
% 
% % 5-fold cross-validation: 
% 
% % (1) calculate lsPFs using the training set,
% % (2) decode the state in each time bin for the PBEs in the test set,
% % (3) calculate the weighted sum of the lsPFs with the weights as the
% %     probability distribution over the states,
% % (4) calculate the center of mass of the final distribution over the
% %     position as the most likely decoded position in the corresponding time
% %     bin
% % (5) finally to calculate position decoding error by finding the distance
% %     between the decoded and actual positions.
% 
% 
% FileBase = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
% mkdir(FileBase)
% 
% 
% binDur = 0.1; % in sec for active wake periods
% numofStates = 40;
% posBinSize = 0.5; % in cm
% 
% 
% % active periods during run
% 
% activePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < behavior.time(2,2), :); 
% % ActivePeriods = ActivePeriods(50:end, :);
% 
% 
% % Binnig the active run period
% 
% [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo);
% 
% activeUnits = 1:size(eventsBinnedfiringPRE{1,2}, 1);
% 
% 
% 
% %% PRE
% 
% 
% pbePeriod = 'PRE';
% 
% % Train an HMM on PBEs
% [transmatPRE, lambdaPRE] = trainHMM(eventsBinnedfiringPRE(:,2), activeUnits, numofStates);
% 
% % normalise_lsPFs = 1;
% plot_lsPF(runBinnedfiring, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% 
% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% 
% [PosDecErrorPRE, PosDecErrorPRE_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatPRE, lambdaPRE, posbinIdx, xposcenters, yposcenters, numofStates);
% 
% 
% 
% %% RUN
% 
% pbePeriod = 'RUN';
% 
% % Train an HMM on PBEs
% [transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);
% 
% % normalise_lsPFs = 1;
% plot_lsPF(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);
% 
% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% 
% [PosDecErrorRUN, PosDecErrorRUN_shuffledlsPFs] = positionDecodingError(runBinnedfiring{2}, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates);
% 
% 
% %% POST
% 
% pbePeriod = 'POST';
% 
% % Train an HMM on PBEs
% [transmatPOST, lambdaPOST] = trainHMM(eventsBinnedfiringPOST(:,2), activeUnits, numofStates);
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
% %% plot decoding error
% 
% figure;
% set(gcf, 'Units', 'pixels', 'position', [814 399 1648 606])
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


%%
% %% Calculating the sequence score over time 
% 
% 
% bvrTimeList = behavior.list(:,[1 2])/Fs;
% bvrState = behavior.list(:,3);
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
% tbeginRUN = behavior.time(2,1)/Fs;
% tendRUN = behavior.time(2,2)/Fs;
% 
% PBEtimesRUN = secondaryPBEs_run(:,1)/Fs;
% 
% 
% 
% %QW periods during RUN
% 
% qw = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% qwRUN = qw(qw(:,1) > tbeginRUN & qw(:,2) < tendRUN, :);
% qwRUN = qwRUN - tbeginRUN;
% 
% 
% subplot(2,1,1)
% PlotSeqFreqvsTime(PBEtimesRUN, HMMprctileRUN, HMMprctile_PRE_RUN, 'RUN', 'PRE', 'RUN', qwRUN, 12, tbeginRUN, tendRUN)
% 
% 
% 
% % (2) time courses during POST (for PRE and RUN(CV))
% 
% 
% tbeginPOST = behavior.time(3,1)/Fs;
% tendPOST = behavior.time(3,2)/Fs;
% 
% PBEtimesPOST = secondaryPBEs_post(:,1)/Fs;
% 
% %nrem periods during POST
% 
% nrem = bvrTimeList(ismember(bvrState, [1,3]),:); % nrem=1, rem=2, quiet=3, wake=4
% nremPOST = nrem(nrem(:,1) > tbeginPOST & nrem(:,2) < tendPOST, :);
% nremPOST = nremPOST - tbeginPOST;
% 
% %
% subplot(2,1,2)
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
% 
% 


%%

% Train on the first 15 min epoch in the post and test on the following
% epochs
% We are extending this to train on all individual epochs on at a time and test on the
% remaining ones


% Train and test on non-overlapping POST epochs

tbegin = behavior.time(3,1); 
tend = behavior.time(3,2);

PBEtimesPOST = secondaryPBEs_post(:,1);

noEpochs = 12; % epochs with duration of ~ 15 min

epochs = linspace(tbegin, tend, noEpochs+1);
epochDur = epochs(2) - epochs(1);

epochCenters = epochs(1:end-1)+epochDur/2 - tbegin;

[~, pbeEpochInd] = histc(PBEtimesPOST, epochs);


numofStates = 40;
PBELikelihood = sequenceScorevstime_CV(eventsBinnedfiringPOST, pbeEpochInd, numofStates, fileinfo);


figure; 
hold on

cmap = colormap('jet');
plotcolors = cmap(floor(linspace(20,64,12)), :);

meanofmean = zeros(noEpochs, 1);

for ii = 1: noEpochs
    
    currLikelihood = PBELikelihood{ii, 1};
    
    meanLL = nan(noEpochs, 1);
    stdLL = nan(noEpochs, 1);
    
    for jj = 1: noEpochs
        
        meanLL(jj) = mean(currLikelihood(pbeEpochInd == jj));
        stdLL(jj) = std(currLikelihood(pbeEpochInd == jj))/ numel(currLikelihood(pbeEpochInd == jj));
        
    end
    
    meanofmean(ii) = nanmean(currLikelihood);
    
    errorbar(epochCenters/3600, meanLL, stdLL, 'color', plotcolors(ii, :), 'linewidth', 2)
    
end
    
    

% Train on PRE and test on POST epochs



% Train on RUN and test on POST epochs





end

end



%% Functions

function sequenceScores = sequenceScorevstime_CV(data, pbeEpochInd, numofStates, fileinfo)



noEpochs = max(pbeEpochInd);

noEvents = size(data, 1);
wholeEvts = 1 : noEvents;

noActiveUnits = size(data{1,1}, 1);

sequenceScores = cell(noEpochs, 1);
% PBELikelihood = cell(noEpochs, 1);

% initialize model parameters

noShuffle = 500;

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);


for epoch = 1: noEpochs
    
    
    % train the model on a single epoch
    
    trainPBEset = find(pbeEpochInd == epoch);
    testPBEset = setdiff(wholeEvts, trainPBEset);
    
    trainData = data(trainPBEset, 2);
    testData = data(testPBEset, :);
    
    [~, ~, curr_transmat, curr_lambda] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 200, 'verbose', 1);
    
    
    
    
    % calcualte the likelihood of the test events (the train events are excluded)
    
%     PBELikelihood{epoch} = nan(noEvents, 1);
%     
%     for evt = 1: length(testData)
%         
%         currEvent = testData{evt};
%         B = poisson_prob(currEvent, curr_lambda,1);
%         prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
% 
%         [~, ~, ~,  PBELikelihood{epoch}(testPBEset(evt)), ~] = fwdback(prior, curr_transmat, B);
%     end


    sequenceScores{epoch} = nan(noEvents, 1);
    sequenceScores{epoch}(testPBEset) = modelCongruence(testData, curr_transmat, curr_lambda, noShuffle, fileinfo, [], [], 0);
    
end


end



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


% For generating pooled time swap data

cnctPBEs = cell2mat(testData(:,2)');

rndgen = randperm(size(cnctPBEs, 2)); 
shuffled_cnctPBEs = cnctPBEs(:, rndgen); % randomizing the order of time bins within the concatenated data  


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

runBinnedfiring2 = timeBinning(activePeriods , mySpikes, qclus, binDur, fileinfo); % the theta periods



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
lsPFs = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmat, lambda, posbinIdx, numofStates, xposcenters, yposcenters);



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
foldSize = floor(noTimeBins/ noFolds);

decodedPositions_actual = zeros(length(posbinIdx), 2);
% decodedPositions_shuffledOrder = zeros(length(posbinIdx), 2);
decodedPositions_shuffledlsPFs = zeros(length(posbinIdx), 2);

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
    
    train_posbinIdx = posbinIdx(trainSet, :); % x and y coordiantes as columns
    
    train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_posbinIdx, numofStates, xposcenters, yposcenters);
    
    
%     train_shuffle_posbinIdx = train_posbinIdx(randperm(size(trainSet, 1), size(trainSet, 1)), :);
    train_shuffle_posbinIdx = train_posbinIdx(randperm(length(trainSet)), :); % the occupancy will stay the same
    
    shuffle_train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_shuffle_posbinIdx, numofStates, xposcenters, yposcenters);
    
    
    % Decode the states in RUN data 
    
    test_timeBins = runBinnedfiring(:,testSet);
    
    obsProb = poisson_prob(test_timeBins, lambda,1);
      
    prior = 1/numofStates * ones(numofStates,1); % a uniform prior
    
    
    %(1) Decoding track positions based on actual lsPFs
    
    decodedPositions_actual(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'actual'); 
    
%     probabilityOverPosition(:,:,testSet) = rr;
    %--- shuffled time bin orders
    
%     decodedPositions_shuffledOrder(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledOrder'); 
    
    
    %--- shuffled gamma: it supposed to take away all the order and
    %coactivity information regarding the states
    
%     decodedPositions_shuffledlsPFs(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledlsPFs'); 
    

    %(2) Decoding track positions based on shuffled lsPFs (lsPFs calculated after shuffling the position indices)
    
    decodedPositions_shuffledlsPFs(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, shuffle_train_lsPFs, xposcenters,yposcenters, 'actual'); 

end


posbinIdx(find(posbinIdx(:,2) > length(yposcenters)), 2) = length(yposcenters);
posbinIdx(find(posbinIdx(:,1) > length(xposcenters)), 1) = length(xposcenters);


actualPositions = [yposcenters(posbinIdx(:,2))' xposcenters(posbinIdx(:,1))'];


PosDecError = sqrt((decodedPositions_actual(:, 1) - actualPositions(:, 1)).^2  + (decodedPositions_actual(:, 2) - actualPositions(:, 2)).^2);
% PosDecError(isnan(PosDecError)) = [];


% PosDecError_shuffledOrder = sqrt((decodedPositions_shuffledOrder(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledOrder(:, 2) - xposcenters(posbinIdx(:,1))').^2);

PosDecError_shuffledlsPFs = sqrt((decodedPositions_shuffledlsPFs(:, 1) - yposcenters(posbinIdx(:,2))').^2  + (decodedPositions_shuffledlsPFs(:, 2) - xposcenters(posbinIdx(:,1))').^2);    
% PosDecError_shuffledlsPFs(isnan(PosDecError_shuffledlsPFs)) = [];



end


function decodedPositions = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, shuffleMode)



noTestBins = size(obsProb, 2);
% numofStates = size(transmat, 1);
shuffledOrder = randperm(noTestBins, noTestBins);

    
if strcmp(shuffleMode, 'shuffledOrder')
    [~, ~, gamma(:, shuffledOrder)] = fwdback(prior, transmat, obsProb(:, shuffledOrder));
else
    [~, ~, gamma] = fwdback(prior, transmat, obsProb);
end


decodedPositions = zeros(noTestBins, 2);
for jj = 1: noTestBins
    
        
%     if strcmp(shuffleMode, 'shuffledGamma')
%         stateProbability = gamma(randperm(numofStates, numofStates), jj);
%     else
%         stateProbability = gamma(:, jj);
%     end
    
    stateProbability = gamma(:, jj);
    
    probabilityOverPosition = sum(repmat(permute(stateProbability, [3 2 1]), length(yposcenters), length(xposcenters)) .* train_lsPFs, 3); % averaging the states' lsPFs weighted by the probability of the states given in the current time bin

    decodedPositions(jj, 1) = yposcenters(ceil(sum((1:length(yposcenters))' .* sum(probabilityOverPosition, 2))));
    decodedPositions(jj, 2) = xposcenters(ceil(sum((1:length(xposcenters)) .* sum(probabilityOverPosition, 1))));
end


end


function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-dpdf','-r0')

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

function [A_B_corresp_pval,  sortAstates] = stateCorrespondence(testData, lastPBEoffirstPeriod, lambdaA, transmatA, lambdaB, transmatB, uniformTranprobs, seqContentCorrespondence)


if ~isempty(lastPBEoffirstPeriod) % for generating the pooled time swap data, I generate the shuffled data separately for each period and them consider them together
    
    testData1 = testData(1:lastPBEoffirstPeriod, :);
    testData2 = testData(lastPBEoffirstPeriod+1 : end, :);
end


noShuffles = 100;

% Decoding the states in the test period given model A and B (trained on the periods A and B, resp.)

% noTimebins = size(testData, 2);

% Decoded states given model A
numofStatesA = size(lambdaA, 2);

if ~uniformTranprobs
    currTransmatA = transmatA;
else
    currTransmatA = 1/numofStatesA * ones(numofStatesA);% to make the transition probabilities irrelevant 
end


priorA = 1/numofStatesA * ones(numofStatesA,1);


nPBEs = size(testData, 1);

gammagivenA = cell(size(testData));
for pbe = 1: nPBEs
    obsProb =  poisson_prob(testData{pbe}, lambdaA, 1);
    [~, ~, gammagivenA{pbe}] = fwdback(priorA, currTransmatA, obsProb);
end

cnctgammagivenA = cell2mat(gammagivenA');

noTimebins = size(cnctgammagivenA, 2);


% Decoded states given model B
numofStatesB = size(lambdaB, 2);

if ~uniformTranprobs
    currTransmatB = transmatB;
else
    currTransmatB = 1/numofStatesB * ones(numofStatesB);
end


priorB = 1/numofStatesB * ones(numofStatesB,1);


gammagivenB = cell(size(testData));
for pbe = 1: nPBEs
    obsProb =  poisson_prob(testData{pbe}, lambdaB, 1);
    [~, ~, gammagivenB{pbe}] = fwdback(priorB, currTransmatB, obsProb);
end

cnctgammagivenB = cell2mat(gammagivenB');




cnctgammagivenA_chance = zeros(numofStatesA, noTimebins, noShuffles);
cnctgammagivenB_chance = zeros(numofStatesB, noTimebins, noShuffles);

switch seqContentCorrespondence
    case 0 
        
        for ii = 1: noShuffles
            cnctgammagivenA_chance(:,:, ii) = cnctgammagivenA(:, randperm(noTimebins));
            cnctgammagivenB_chance(:,:, ii) = cnctgammagivenB(:, randperm(noTimebins));
        end
        
    case 1
        
        % obs probs given A
        obsProb_A1 = cell(size(testData1));
        for pbe = 1: size(testData1, 1) % we calculate the observation probabilities only once and for the shuffle data we shuffle the obs probs instead of the time-binned firing rates
            obsProb_A1{pbe} =  poisson_prob(testData1{pbe}, lambdaA, 1);
        end

        obsProb_A2 = cell(size(testData2));
        for pbe = 1:size(testData2, 1)
            obsProb_A2{pbe} = poisson_prob(testData2{pbe}, lambdaA, 1);
        end

        % obs probs given B
        obsProb_B1 = cell(size(testData1));
        for pbe = 1: size(testData1, 1) 
            obsProb_B1{pbe} =  poisson_prob(testData1{pbe}, lambdaB, 1);
        end

        obsProb_B2 = cell(size(testData2));
        for pbe = 1:size(testData2, 1)
            obsProb_B2{pbe} = poisson_prob(testData2{pbe}, lambdaB, 1);
        end



        for ii = 1: noShuffles  % it's like pooled time swap
            ii
            [obsProb_A1_ts, shuffleOrder1] = genTimeSwap(obsProb_A1);
            [obsProb_A2_ts, shuffleOrder2] = genTimeSwap(obsProb_A2);
            obsProb_A_ts = [obsProb_A1_ts; obsProb_A2_ts];

            obsProb_B1_ts = genTimeSwap(obsProb_B1, shuffleOrder1);
            obsProb_B2_ts = genTimeSwap(obsProb_B2, shuffleOrder2);
            obsProb_B_ts = [obsProb_B1_ts; obsProb_B2_ts];


            % given model A
            nPBEs = size(obsProb_A_ts, 1);

            gammagivenA_chance = cell(size(obsProb_A_ts));
            for pbe = 1: nPBEs
                [~, ~, gammagivenA_chance{pbe}] = fwdback(priorA, currTransmatA, obsProb_A_ts{pbe});
            end
            cnctgammagivenA_chance(:,:, ii) = cell2mat(gammagivenA_chance');


            % given model B
            nPBEs = size(obsProb_B_ts, 1);

            gammagivenB_chance = cell(size(obsProb_B_ts));
            for pbe = 1: nPBEs
                [~, ~, gammagivenB_chance{pbe}] = fwdback(priorB, currTransmatB, obsProb_B_ts{pbe});
            end
            cnctgammagivenB_chance(:,:, ii) = cell2mat(gammagivenB_chance');
        end
        
    case 2
        
        shuffledTmatA = genshuffles(currTransmatA, noShuffles, []);
        shuffledTmatB = genshuffles(currTransmatB, noShuffles, []);
        for ii = 1: noShuffles  % transition matrix shuffle
            ii
            % gamma for shuffled data given model A
            gammagivenA_chance = cell(size(testData));
            for pbe = 1: nPBEs
                obsProb =  poisson_prob(testData{pbe}, lambdaA, 1);
                [~, ~, gammagivenA_chance{pbe}] = fwdback(priorA, shuffledTmatA(:,:,randi(noShuffles)), obsProb);
            end
            cnctgammagivenA_chance(:,:, ii) = cell2mat(gammagivenA_chance');

            % gamma for shuffled data given model B
            gammagivenB_chance = cell(size(testData));
            for pbe = 1: nPBEs
                obsProb =  poisson_prob(testData{pbe}, lambdaB, 1);
                [~, ~, gammagivenB_chance{pbe}] = fwdback(priorB, shuffledTmatB(:,:,randi(noShuffles)), obsProb);
            end
            cnctgammagivenB_chance(:,:, ii) = cell2mat(gammagivenB_chance');

        end
end
    

A_B_corresp = (cnctgammagivenA * cnctgammagivenB')/noTimebins; % A states in rows and B states in columns

A_B_corresp_chance = zeros(numofStatesA, numofStatesB, noShuffles);

for ii = 1:noShuffles
    A_B_corresp_chance(:,:, ii) = (cnctgammagivenA_chance(:,:, ii) * cnctgammagivenB_chance(:,:, ii)')/noTimebins;
end


A_B_corresp_pval = zeros(numofStatesA, numofStatesB);
for istate = 1: numofStatesA
    for jstate = 1:numofStatesB
        
        [~, pval] = ttest2(A_B_corresp(istate, jstate), squeeze(A_B_corresp_chance(istate, jstate, :)), 'Tail', 'right', 'alpha', 0.01);
        A_B_corresp_pval(istate, jstate) = -log10(pval);
        
    end
end


% sorting the states in one period based on the maximum corresponding
% states from anoter period

[~, maxBstates] = max(A_B_corresp_pval, [], 2);

[~, sortAstates] = sort(maxBstates);

end





function plotPanel(stateCorrespMat, rowSortInd, model, model2comp, ccodeRange, currTitle)


% 
% for t = 1:5
%     [i,j] = find(stateCorrespMat == max(stateCorrespMat(:)));
%     stateCorrespMat(i,j) = mean(stateCorrespMat(:));
% end

if ~isempty(ccodeRange)
    
    imagesc(stateCorrespMat(rowSortInd, :), ccodeRange); colormap('jet'); % 
else
    imagesc(stateCorrespMat(rowSortInd, :)); colormap('jet'); % 
end


set(gca, 'YDir', 'normal', 'ytick', 1:length(rowSortInd), 'yticklabel', rowSortInd, 'fontsize', 10)
set(gca, 'XTick', [], 'YTick', [])

ylabel([model ' states'] , 'fontsize', 10)
xlabel([model2comp ' states'], 'fontsize', 10)

title(currTitle, 'fontsize', 10, 'fontweight', 'normal')
% colorbar

end



function plotErrorCdf(data, shuffledGamma, ylabelneeded)
% plotErrorCdf(data, shuffledOrder, shuffledGamma, ylabelneeded)

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
[h3, bins] = hist(shuffledGamma, 100);
h3c = cumsum(h3)/sum(h3);

shuffleMed = median(shuffledGamma);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(shuffleMed, 0.1, sprintf('%.1f cm', shuffleMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')




xlim([0 max([data; shuffledGamma])])
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


function out = Lambdasimilarity(x,y)

out = diag(x' * y ./ (sum(x)'*sum(y)));

end


function [PopVectSim, sortAstates] = PopVecSimilarity(x, y)

PopVectSim = x' * y ./ (sum(x)'*sum(y));


[~, maxBstates] = max(PopVectSim, [], 2);
[~, sortAstates] = sort(maxBstates);

% PopVectSim = PopVectSim(sortAstates, :);

end


function plotHistCorrespondence(data, chance, firstPeriod, secondPeriod, currTitle)


bins = linspace(min([data; chance]), max([data; chance]), 50);
h = hist(data, bins)/length(data);


chance = chance(randi(length(chance), 5000, 1));
h_chance = hist(chance, bins)/length(chance);
pval = ranksum(data, chance, 'tail', 'right');


set(gca, 'fontsize', 10, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.005, 0.01])

hold on

bar(bins, h_chance, 'FaceColor', 'none', 'EdgeColor', 'k') 
bar(bins, h, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5)


h_chance_sm = conv(h_chance, gausswindow(2,5), 'same');
h_sm = conv(h, gausswindow(2,5), 'same');
plot(bins, h_chance_sm, 'color', [0.7 0.7 0.7], 'linewidth', 3)
plot(bins, h_sm, 'k', 'linewidth', 3)

ylim([-0.01 0.2])

ylim1=get(gca,'ylim');
xlim1=get(gca,'xlim');

if pval < 0.001
    text(mean(xlim1),mean(ylim1), 'p < 0.001', 'fontsize', 10)
else
    text(mean(xlim1),mean(ylim1), sprintf('p = %.2f', pval), 'fontsize', 10)
end

xlabel('corres. w/in third period', 'fontsize', 10)
title({currTitle ; '' ; [firstPeriod ' and ' secondPeriod]}, 'fontsize', 10, 'fontweight', 'normal')

end

function [outPBEs, shuffleOrder] = genTimeSwap(inPBEs, varargin)

noPBEs = size(inPBEs, 1);

cnctPBEs = cell2mat(inPBEs');

if nargin > 1
    shuffleOrder = varargin{1};
else
    shuffleOrder = randperm(size(cnctPBEs, 2));
end

Tswap_cnctPBEs = cnctPBEs(:, shuffleOrder);

outPBEs = cell(size(inPBEs));
for pbe = 1: noPBEs
    
    nTimeBins = size(inPBEs{pbe}, 2);
    
    outPBEs{pbe} = Tswap_cnctPBEs(:, 1:nTimeBins);
    Tswap_cnctPBEs(:, 1:nTimeBins) = [];
    
end    

end
