lear; clc; close all

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {1, 1};

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



% RUN

fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [1 2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);


[poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, temporalEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);



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



% RUN
pbePeriod = 'RUN';
activeUnits = 1:size(eventsBinnedfiringPRE{1,2}, 1);

% Train an HMM on PBEs
[transmatRUN, lambdaRUN] = trainHMM(eventsBinnedfiringRUN(:,2), activeUnits, numofStates);


lsPFs = plot_lsPF(runBinnedfiring, transmatRUN, lambdaRUN, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, behavior, pbePeriod);


%% HMM sequence detection

% In this step we test the extent to which an HMM trained on the
% candidate events is able to detect sequence events among non-sequence
% events

% we do this by k-fold cross-validation: using k-1 folds for training and
% and the remaining fold for sequence detection %% this is temporary assigning. Later will change it when analyzing run/post-sleep data


numofFolds = 5; %% 5-fold cross-validation
numofStates = 40;
noShuffle = 2000;



% RUN cross-validation

pbePeriod = 'RUN';

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiringRUN, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctile_run = HMMcongruence_crossvalid(eventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
HMMprctile_run_p = HMMcongruence_crossvalid(poissonEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
HMMprctile_run_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
HMMprctile_run_t = HMMcongruence_crossvalid(temporalEventsBinnedfiringRUN, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);



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


function [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behavior, binDur, posBinSize, fileinfo)


%% binning the run data


qclus = [1 2 3 9]; % pyramidal and interneurons

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
