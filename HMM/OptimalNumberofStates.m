clear; clc; close all

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[2 3], [1 2 3], 1};

for rr = 1:length(rats) % loop for all the rats and sessions
    
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

% exclude all the periods which we are not intersetd for now
exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_pre = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPRE, secondaryPBEs_pre, acceptedprmEvts_pre] = finalBinningResult(primaryPBEs_pre, mySpikes, fileinfo);



% RUN

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [1 2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);



% POST

fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);


exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_post = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringPOST, secondaryPBEs_post, acceptedprmEvts_post] = finalBinningResult(primaryPBEs_post, mySpikes, fileinfo);




%% Determining the optimal number of states




% %% Within each period
% 
% 
searchSet = 5:5:60;
iter = 10;
searchSet = repelem(searchSet, iter);

setSize = length(searchSet);
% 
% numofFolds = 5;
% 
% % PRE cross validation
% 
% numofPBEs = size(eventsBinnedfiringPRE, 1);
% 
% likelihoodPRE = zeros(numofPBEs, setSize);
% 
% for ii = 1: length(searchSet) 
% 
%      numofStates = searchSet(ii);
%     [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPRE, numofStates, numofFolds, fileinfo, []); 
%     
%     likelihoodPRE(:, ii) = calculateLikelihoods(eventsBinnedfiringPRE, testEvts, CVtransmats, CVlambdas);
%     
% end
% 
% 
% 
% % RUN cross-validation
% 
% 
% numofPBEs = size(eventsBinnedfiringRUN, 1);
% 
% likelihoodRUN = zeros(numofPBEs, setSize);
% 
% for ii = 1: length(searchSet) 
%     
%     numofStates = searchSet(ii);
%     [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN, numofStates, numofFolds, fileinfo, []); 
%     
%     likelihoodRUN(:, ii) = calculateLikelihoods(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas);
%     
% end
% 
% 
% 
% % POST cross-validation
% 
% 
% numofPBEs = size(eventsBinnedfiringPOST, 1);
% 
% likelihoodPOST = zeros(numofPBEs, setSize);
% 
% for ii = 1: length(searchSet) 
% 
%     numofStates = searchSet(ii);
%     [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringPOST, numofStates, numofFolds, fileinfo, []); 
%     
%     likelihoodPOST(:, ii) = calculateLikelihoods(eventsBinnedfiringPOST, testEvts, CVtransmats, CVlambdas);
%     
% end
% 
% 



%% Between the periods


noActiveUnits = size(eventsBinnedfiringPRE{1,2}, 1);

% prior0 = normalise(rand(numofStates,1));
% transmat0 = mk_stochastic(rand(numofStates,numofStates));
% lambda0 = rand(noActiveUnits, numofStates);


% train on PRE, test on RUN or POST

likelihoodPRE_RUN = zeros(size(eventsBinnedfiringRUN, 1), setSize);
likelihoodPRE_POST = zeros(size(eventsBinnedfiringPOST, 1), setSize);


fprintf('\n%s _ train on PRE, test on RUN & POST started...\n', sessionName)

for ii = 1: length(searchSet)
    
    numofStates = searchSet(ii);
    
    if mod(ii, 10) == 1
       fprintf('\nNumber of states = %d .', numofStates)
    else
        fprintf(' .')
    end
    
   
    prior0 = normalise(rand(numofStates,1));
    transmat0 = mk_stochastic(rand(numofStates,numofStates));
    lambda0 = rand(noActiveUnits, numofStates);

    [~, ~, transmatPRE, lambdaPRE] = phmm_em(eventsBinnedfiringPRE(:,2), prior0, transmat0, lambda0, 'max_iter', 100, 'verbose', 0);
    
   likelihoodPRE_RUN(:, ii) = modelCongruence(eventsBinnedfiringRUN, transmatPRE, lambdaPRE);
   
   likelihoodPRE_POST(:, ii) = modelCongruence(eventsBinnedfiringPOST, transmatPRE, lambdaPRE);
    
end


% train on RUN, test on POST

numofPBEs = size(eventsBinnedfiringPOST, 1);

likelihoodRUN_POST = zeros(numofPBEs, setSize);

fprintf('\n%s _ train on RUN, test on POST started...', sessionName)

for ii = 1: length(searchSet)
    
    numofStates = searchSet(ii);
    
    if mod(ii, 10) == 1
       fprintf('\nNumber of states = %d .', numofStates)
    else
        fprintf(' .')
    end
    
    prior0 = normalise(rand(numofStates,1));
    transmat0 = mk_stochastic(rand(numofStates,numofStates));
    lambda0 = rand(noActiveUnits, numofStates);
   [~, ~, transmatRUN, lambdaRUN] = phmm_em(eventsBinnedfiringRUN(:,2), prior0, transmat0, lambda0, 'max_iter', 100, 'verbose', 0);
   
   likelihoodRUN_POST(:, ii) = modelCongruence(eventsBinnedfiringPOST, transmatRUN, lambdaRUN);
 
end

% 
% % train on PRE, test on POST
% 
% numofPBEs = size(eventsBinnedfiringPOST, 1);
% 
% likelihoodPRE_POST = zeros(numofPBEs, setSize);
% 
% fprintf('\n%s _ train on PRE, test on POST started...', sessionName)
% 
% for ii = 1: length(searchSet)
%     
%     numofStates = searchSet(ii);
%     
%     if mod(ii, 10) == 1
%        fprintf('\nNumber of states = %d .', numofStates)
%     else
%         fprintf(' .')
%     end
%     
%     prior0 = normalise(rand(numofStates,1));
%     transmat0 = mk_stochastic(rand(numofStates,numofStates));
%     lambda0 = rand(noActiveUnits, numofStates);
%     [~, ~, transmatPRE, lambdaPRE] = phmm_em(eventsBinnedfiringPRE(:,2), prior0, transmat0, lambda0, 'max_iter', 100, 'verbose', 0);
% 
%     likelihoodPRE_POST(:, ii) = modelCongruence(eventsBinnedfiringPOST, transmatPRE, lambdaPRE);
% 
%     
% end
% 




figure;
set(gcf, 'position', [801 697 1648 624])

% cross-validations (within-periods)
% 
% subplot(2,3,1)
% likelihoodVSNumofStates(searchSet, likelihoodPRE, 'PRE (CV)')
% 
% subplot(2,3,2)
% likelihoodVSNumofStates(searchSet, likelihoodRUN, 'RUN (CV)')
% 
% subplot(2,3,3)
% likelihoodVSNumofStates(searchSet, likelihoodPOST, 'POST (CV)')

% congruence of PBEs with other periods' models

subplot(1,3,1)
likelihoodPRE_RUN = reshape(likelihoodPRE_RUN, length(likelihoodPRE_RUN)*iter, setSize/iter);
likelihoodVSNumofStates(5:5:60, likelihoodPRE_RUN, 'Train on PRE and test on RUN')

subplot(1,3,2)
likelihoodRUN_POST = reshape(likelihoodRUN_POST, length(likelihoodRUN_POST)*iter, setSize/iter);
likelihoodVSNumofStates(5:5:60, likelihoodRUN_POST, 'Train on RUN and test on POST')

subplot(1,3,3)
likelihoodPRE_POST = reshape(likelihoodPRE_POST, length(likelihoodPRE_POST)*iter, setSize/iter);
likelihoodVSNumofStates(5:5:60, likelihoodPRE_POST, 'Train on PRE and test on POST')



currDir = pwd;
FileBase = [currDir '/' fileinfo.name  '/HMM/optimalnumofStates/'];
mkdir(FileBase)


print(gcf, [FileBase fileinfo.name '_likelihood_vs_noStates'], '-dpng')
savepdf(gcf, [FileBase fileinfo.name '_likelihood_vs_noStates'])


save([FileBase fileinfo.name '_likelihood_vs_noStates.mat'], ...
    'likelihoodPRE_RUN', 'likelihoodRUN_POST', 'likelihoodPRE_POST')
%     likelihoodPRE, likelihoodRUN, likelihoodPOST, ...
    
 

end
end





%% Functions


function [eventsBinnedfiring, secondaryPBEs, idx_acceptedEvents] = finalBinningResult(pbEvents2, mySpikes, fileinfo)

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
secondaryPBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:4)]; 


eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


end



function likelihood = calculateLikelihoods(data, testEvts, transmat, lambda)



noEvents = size(data, 1);
numofFolds = length(testEvts);
numofStates = size(transmat, 1);


likelihood = zeros(1, noEvents); %% log likelihood of the raw data


for fold = 1 : numofFolds
    
    curr_transmat = transmat(:, :, fold);
    curr_lambda = lambda(:, :, fold);
    
    
    for evt = testEvts{fold}
        
        currEvent = data{evt, 2};

        B = poisson_prob(currEvent, curr_lambda,1);
        prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, ~,  likelihood(evt), ~] = fwdback(prior, curr_transmat, B);
        
    end
end

end





function likelihood = modelCongruence(testData, transmat, lambda)


noEvents = size(testData, 1);
numofStates = size(transmat, 1);


likelihood = zeros(1, noEvents); %% log likelihood of the raw data

for evt = 1 : noEvents
    
    currEvent = testData{evt, 2};
%     noTimebins = size(currEvent, 2);

    B = poisson_prob(currEvent, lambda,1);
    prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states

    %%% observation likelihood for the actual data
    [~, ~, ~,  likelihood(evt), ~] = fwdback(prior, transmat, B);
end

end



function likelihoodVSNumofStates(searchSet, likelihood, currTitle)


noPBEs = size(likelihood, 1);

means = mean(likelihood);
% means = means-means(1);

% likelihood = likelihood/means(1);
errs = std(likelihood)/noPBEs;


h = errorbar(searchSet, means, errs, 'color', 'k', 'linewidth', 2);
set(gca, 'fontsize', 14, 'linewidth', 2, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'box', 'off', 'Layer', 'Top')

xlabel('Number of States', 'fontsize', 14)
ylabel('Likelihood of PBEs', 'fontsize', 14)

axis square

title(currTitle, 'fontsize', 14)


end



function savepdf(gcf, pdfname)

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf, pdfname,'-depsc','-r0')

end
