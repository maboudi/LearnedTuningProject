% Train the model on a period and test on another period (for example, pre-sleep and run periods)

clear; clc;

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

animal = input('\n Enter the animal name (in quotes) \n');
animal(1) = upper(animal(1));


sessionNumber = input('\n Enter the session number \n');
% sessionNumber = 2;

sessionName = [animal 'Maze' num2str(sessionNumber)];
sessionName2 = [animal '-maze' num2str(sessionNumber)];

VarList = {'spikes','behavior','position','speed','basics','ripple'};

for var = 1 : length(VarList)
    load([currDir '/pooled/wake-' VarList{var} '.mat'])
end

spikes = eval(['spikes.' sessionName]);

behavior = eval(['behavior.' sessionName]);
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

fileinfo = struct('name', sessionName2, 'animal', animal, 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', Fs, 'nCh', nCh, 'lfpSampleRate', lfpSampleRate, 'pix2cm', 0.28); 
fileinfo.pix2cm = 0.28; % in cm

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

breakLap = input('\n Enter the lap number for truncating the data\n(If no partitioning is needed just press enter) \n');  

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
activeTimes = bvrTimeList(bvrState == 4,:);
noActEpochs = length(activeTimes);


%% PBEs

    %% detemining population burst periods




time_resolution = 0.001; % in second
threshZ = 3; % sdf with 3 std deviation above the mean


% pre-run sleep period

fileinfo.tbegin = behavior.time(1,1) ; %% this is temporary assigning. Later will change it when analyzing run/post-sleep data
fileinfo.tend = behavior.time(1,2);


exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4

pbEvents_pre = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);


eventsBinnedfiringPRE = finalBinningResult(pbEvents_pre, mySpikes, fileinfo);
[poissonEventsBinnedfiringPRE, timeSwapEventsBinnedfiringPRE, cycleEventsBinnedfiringPRE] = genSurrogates(eventsBinnedfiringPRE);




% run period

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);

exclude = bvrTimeList(ismember(bvrState, [1 2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

pbEvents_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);


eventsBinnedfiringRUN = finalBinningResult(pbEvents_run, mySpikes, fileinfo);
[poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, cycleEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);


% post-run sleep period

fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);


exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4

pbEvents_POST = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);


eventsBinnedfiringPOST = finalBinningResult(pbEvents_POST, mySpikes, fileinfo);
[poissonEventsBinnedfiringPOST, timeSwapEventsBinnedfiringPOST, cycleEventsBinnedfiringPOST] = genSurrogates(eventsBinnedfiringPOST);






%% HMM sequence detection

%%% In this step we test the extent to which an HMM trained on the
%%% candidate events is able to detect sequence events among non-sequence
%%% events

%%% we do this by k-fold cross-validation: using k-1 folds for training and
%%% and the remaining fold for sequence detection 


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



noShuffle = 1000;





%%% train PRE___ test RUN

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;

trainPeriod = 'PRE';
testPeriod = 'RUN';
curr_title = ['Train on ' trainPeriod '_ test on ' testPeriod];

HMMPercentile_pre_run = modelCongruence(eventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates

HMMPercentile_pre_run_p = modelCongruence(poissonEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMPercentile_pre_run_ts = modelCongruence(timeSwapEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMPercentile_pre_run_t = modelCongruence(cycleEventsBinnedfiringRUN, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');



figure;
set(gcf, 'position', [40 480 2527 686])

subplot(1,3,1)
congruencePlot(HMMPercentile_pre_run, HMMPercentile_pre_run_p, HMMPercentile_pre_run_ts, HMMPercentile_pre_run_t, 1, curr_title)


%%% train RUN___ test POST

curr_transmat = transmatRUN;
curr_lambda = lambdaRUN;

trainPeriod = 'RUN';
testPeriod = 'POST'; 
curr_title = ['Train on ' trainPeriod '_ test on ' testPeriod];

HMMPercentile_run_post = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates


HMMPercentile_run_post_p = modelCongruence(poissonEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMPercentile_run_post_ts = modelCongruence(timeSwapEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMPercentile_run_post_t = modelCongruence(cycleEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');

subplot(1,3,2)
congruencePlot(HMMPercentile_run_post, HMMPercentile_run_post_p, HMMPercentile_run_post_ts, HMMPercentile_run_post_t, 0, curr_title)



%%% train PRE___ test POST

curr_transmat = transmatPRE;
curr_lambda = lambdaPRE;


trainPeriod = 'PRE';
testPeriod = 'POST'; 
curr_title = ['Train on ' trainPeriod '_ test on ' testPeriod];

HMMPercentile_pre_post = modelCongruence(eventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'actual');

% _____________surrogates

HMMPercentile_pre_post_p = modelCongruence(poissonEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'poisson');
HMMPercentile_pre_post_ts = modelCongruence(timeSwapEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'timeswap');
HMMPercentile_pre_post_t = modelCongruence(cycleEventsBinnedfiringPOST, curr_transmat, curr_lambda, noShuffle, fileinfo, trainPeriod, testPeriod, 'temporal');


subplot(1,3,3)
congruencePlot(HMMPercentile_pre_post, HMMPercentile_pre_post_p, HMMPercentile_pre_post_ts, HMMPercentile_pre_post_t, 0, curr_title)


currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/crossTest/'];
mkdir(FileBase)



print(gcf, [FileBase fileinfo.name '_crossTest'], '-dpng')
print(gcf, [FileBase fileinfo.name '_crossTest'], '-dpdf')

save([FileBase fileinfo.name '_crossTest.mat'], 'HMMPercentile_pre_run', 'HMMPercentile_pre_run_p', 'HMMPercentile_pre_run_ts', 'HMMPercentile_pre_run_t', ...
    'HMMPercentile_run_post', 'HMMPercentile_run_post_p', 'HMMPercentile_run_post_ts', 'HMMPercentile_run_post_t', ...
    'HMMPercentile_pre_post', 'HMMPercentile_pre_post_p', 'HMMPercentile_pre_post_ts', 'HMMPercentile_pre_post_t')


function eventsBinnedfiring = finalBinningResult(pbEvents2, mySpikes, fileinfo)

    %% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events


%%% Calculating the number of firing pyramidal units within each event

qclus = [1 2 3 9]; % pyramidal and interneurons

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

eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


end


function [poissonEventsBinnedfiring, timeSwapEventsBinnedfiring, cycleEventsBinnedfiring] = genSurrogates(eventsBinnedfiring)


%% Generating shuffle surrogate PBE datasets

% generate poisson simulated dataset

binDur = 0.02;
noEvents = size(eventsBinnedfiring, 1);

poissonEventsBinnedfiring = poissonSpikeTrain(eventsBinnedfiring, binDur); 
% cmprPoisson2Raw(poissonEventsBinnedfiring, eventsBinnedfiring, longORshort, fileinfo, binDur); %%% compare the poisson simulated data with the actual data distribution

% generate time swap dataset (coherent shuffle within each event, keeping the cofirings the same)

timeSwapEventsBinnedfiring = timeswap(eventsBinnedfiring, binDur);


% time bins cycle shuffle (independently for each unit, incoherent shuffle)

cycleEventsBinnedfiring = cell(noEvents, 1);
for evt = 1 : noEvents
    
    currEvt = eventsBinnedfiring(evt, :);
    
    randShift = randi(size(currEvt{2}, 2)-1, 1, size(currEvt{2}, 1)); % generate shift amount for each row
    
    
    % Because we are doing row cycle shuffle, we need to transpose the data
    % matrix and after doing column cycle shuffles transpose it again
    
    cycleEventsBinnedfiring{evt, 2} = column_cycle_shuffle(currEvt{2}', randShift)';
    cycleEventsBinnedfiring{evt, 1} = column_cycle_shuffle(currEvt{1}', randShift*20)'; % note the shift here
    
end

end



function HMMPercentile = modelCongruence(testData, transmat, lambda, noShuffle, fileinfo, trainPeriod, testPeriod, surrogateData)



noEvents = size(testData, 1);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, noShuffle); %% log likelihood of shuffle data

HMMPercentile = zeros(noEvents, 1);


for evt = 1 : noEvents
    
    currEvent = testData{evt, 2};
%     noTimebins = size(currEvent, 2);

    B = poisson_prob(currEvent, lambda,1);
    prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the actual data
    [~, ~, ~,  dataLL(evt), ~] = fwdback(prior, transmat, B);



    %%%% null distribution of observation likelihhod
    %%%% shuffling the transmat matrix

    shuffled_transmat = zeros(size(transmat));

    for sn = 1 : noShuffle

        for s = 1 : numofStates

            %%% redistribute the transition probabilities within each 
            %%% row with an exception that the self-transitions (matrix diagonal elements) are not changed

            randgen = randperm(numofStates)';
            randgen(randgen == s) = [];

            shuffleInd = [randgen(1:s-1); s; randgen(s : end)]; 

            shuffled_transmat(:,s) = transmat(s, shuffleInd); 

        end             

        [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, shuffled_transmat, B); %%% B and Prior are the same as in the case of raw model
    end




    for method = 1%:2


        HMMPercentile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/noShuffle;

    end
    
    if mod(evt, 20) == 0
        fprintf(1, [trainPeriod '-' testPeriod '-' surrogateData '__event %d, HMMPercentile = %f\n'], evt, HMMPercentile(evt, 1));
    end
    
    
end



currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/congruence/'];
mkdir(FileBase)


save([FileBase '/Train_' trainPeriod '___Test_' testPeriod '__' surrogateData '.mat'], 'HMMPercentile', 'dataLL', 'nullLL')

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


p1 = plot(bins, cumhist, 'color', [50 150 50]/255, 'linewidth', 2);
line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [50 150 50]/255,'LineStyle', '--', 'linewidth', 1)

p2 = plot(bins, cumhist_ts, 'color', [30 144 200]/255, 'linewidth', 2);
line([0 bins(threshBin)], [crosspnt_ts crosspnt_ts], 'color', [30 144 200]/255,'LineStyle', '--', 'linewidth', 1)

p3 = plot(bins, cumhist_p, 'color', [180 150 110]/255, 'linewidth', 2);
line([0 bins(threshBin)], [crosspnt_p crosspnt_p], 'color', [180 150 110]/255,'LineStyle', '--', 'linewidth', 1)

p4 = plot(bins, cumhist_tc, 'color', [150 0 0]/255, 'linewidth', 2);
line([0 bins(threshBin)], [crosspnt_tc crosspnt_tc], 'color', [150 0 0]/255,'LineStyle', '--', 'linewidth', 1)

noEvents = length(HMMprctile);
line([bins(threshBin) bins(threshBin)], [0 noEvents], 'color', 'k','LineStyle', '--', 'linewidth', 1)


if okLegend
    legend([p1 p2 p3 p4],{'actual','time swap','poisson','temporal shuffle'})
end

set(gca, 'fontsize', 14)
xlabel('percentile score', 'fontsize', 20)

if okLegend
ylabel('cumulative number of PBEs', 'fontsize', 20)
end

title(curr_title, 'fontsize', 20)


end
