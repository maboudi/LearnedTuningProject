clear; clc;

currDir = '/home/kouroshmaboudi/Documents/HMM project';
cd(currDir)


%% loading session data

animal = input('\n Enter the animal name (in quotes) \n');
animal(1) = upper(animal(1));


% sessionNumber = input('\n Enter the session number \n');
sessionNumber = 1;

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


% %%% truncate the spikes if partitioning is desired

% runSpikeTimes = mySpikes.t(find(mySpikes.lap > 0));
% fileinfo.tbegin = runSpikeTimes(1);
% fileinfo.tend = runSpikeTimes(end);


% if ~isempty(breakLap)
%     
%     breakpoint = mySpikes.t(find(mySpikes.lap == breakLap, 1, 'first'));
%     mySpikes = truncateSpike(mySpikes, firstORsecond, breakpoint, fileinfo.tbegin, fileinfo.tend, Fs);
%     
%     if firstORsecond == 2
%        mySpikes.lap = mySpikes.lap - breakLap;
%     end
%     
% end

%% PBEs

    %% detemining population burst periods


% pre-run sleep period

% fileinfo.tbegin = behavior.time(1,1) ; %% this is temporary assigning. Later will change it when analyzing run/post-sleep data
% fileinfo.tend = behavior.time(1,2);


% run period

% fileinfo.tbegin = behavior.time(2,1);
% fileinfo.tend = behavior.time(2,2);


% post-run sleep period

fileinfo.tbegin = behavior.time(3,1); 
fileinfo.tend = behavior.time(3,2);


% changing the zero time of spike timestamps 

% mySpikes.t = mySpikes.t - fileinfo.tbegin;
% mySpikes.t = mySpikes.t * Fs/1e6;


% Theta periods to be excluded (the timing information should be in second)

exclude = bvrTimeList(ismember(bvrState, [2 3 4]),:); % nrem=1, rem=2, quiet=3, wake=4
% exclude = [exclude; [8.52e10 fileinfo.tend]]; % just in case of the run PBEs
% exclude = (exclude - fileinfo.tbegin)/1e6;


time_resolution = 0.001; % in second
threshZ = 3; % sdf with 3 std deviation above the mean


[pbEvents2, spikeDensity] = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);


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
eventsBinnedfiring2 = timeBinning(pbEvents2, mySpikes, qclus, binDur, Fs);


% Remove the flanking zero-firing periods(silent bins) from the beginnig
% and end of each event. Since the interneurons were involved in calculating
% the PBEs' boundaries and here we are binning just pyramidal units' spike
% trains, some silent bins are expected. 

[eventsBinnedfiring, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring2, pbEvents2(:, 1), binDur, Fs);


% Qualify the events based on the minumum active units and length criteria
% Considering only the event with number of active units firing above 10 percent of the total or 5 (whichever is greater) and
% duration of at least 4 20ms-time bins.


activeUnits = unique(mySpikes.unit);
% idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4 & eventLen <= 25)); % maximum length 500 ms (???)
idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4 & eventLen <= 25)); % maximum length 500 ms (???)

eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


% Final population burst events

PBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:end)];
% PBEs = PBEs(1:1500, :);

noEvents = size(PBEs, 1)


% Generate output files

currDir = pwd;
Folderbase = [currDir '/' fileinfo.name  '/Population Burst Events'];
mkdir(Folderbase)

save([Folderbase '/' 'PBEvents.mat'], 'PBEs')

Filename = [Folderbase '/' fileinfo.name '.pbe.evt'];
MakeEvtFile(PBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, Fs, 1)


% Plot spike density function and the the boundaries of PBEs

figure;

hold on

plot((0:(length(spikeDensity)-1))*time_resolution, spikeDensity, 'linewidth', 1, 'color','k')

yL = get(gca,'YLim');

PBEs_sec = (PBEs(:,1:2)-fileinfo.tbegin)*1/Fs;

for epoch = 1 : size(PBEs_sec,1)
    
    p = patch([PBEs_sec(epoch,1) PBEs_sec(epoch,2) PBEs_sec(epoch,2) PBEs_sec(epoch,1)],[yL(1) yL(1) yL(2) yL(2)], ...
       'r','edgecolor', 'none');
    set(p,'FaceAlpha', 0.5)
    
end

plot((1:length(spikeDensity))*time_resolution, ones(1,length(spikeDensity))*threshZ, '--r', 'linewidth', 2)

hold off

ylabel('Firing rate(Hz)', 'fontsize', 20)
xlabel('time(sec)', 'fontsize', 20)
set(gca,'fontsize', 16);

saveas(gcf, [Folderbase '/spikeDensityFunction.fig'])



%% Generating shuffle surrogate PBE datasets

% generate poisson simulated dataset

binDur = 0.02;

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


%% Analyzing the hidden states regarding sparsity, and state place fields
% 
% FileBase = [currDir '/' fileinfo.name '/HMM/sparsity'];
% mkdir(FileBase)
% 
% 
% numofStates = 40; % we have fixed the number of states to 30 here
% numIter = 100; % because different initializations lead to differnt results we do the training several times and calculate the distribution
% 
% 
% % actual data
% dataType = 'actual'; 
% modelSparsity(eventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)
% 
% 
% %%%% train HMMs on surrogate datasets %%%
% 
% % (1)poisson simulated spike trains
% dataType = 'poisson';
% modelSparsity(poissonEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)
% 
% % (2)incoherent time bin shuffle
% dataType = 'timeCycle';
% modelSparsity(cycleEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)
% 
% % (3)time swap
% dataType = 'timeswap';
% modelSparsity(timeSwapEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)



%% lsPF and decoding of position using PBE models



%% binning run data


binDur = 0.25; % in sec 

% active periods during run
% ActivePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < 8.52e10, :);
ActivePeriods = bvrTimeList(bvrState == 4 & bvrTimeList(:,1) >= behavior.time(2,1) & bvrTimeList(:,2) < behavior.time(2,2), :);

% ActivePeriods = ActivePeriods(50:end, :);

runBinnedfiring2 = timeBinning(ActivePeriods , mySpikes, qclus, binDur, Fs); % the theta periods


% Calculate the track positions corresponding to each time bin

periods = ActivePeriods* 1/lfpSampleRate;
runBinPos = cell(size(periods, 1), 1);

for ii = 1 : size(periods, 1)
    
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/lfpSampleRate, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2);
    runBinPos{ii}(end, :) = [];
end
    

% runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields

posBinSize = 0.5; % in cm

xpos = fileinfo.xyt(:,1);
ypos = fileinfo.xyt(:,2);


noXPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
noYPosBins = floor((max(ypos) - min(ypos))/posBinSize);

xposBins = min(xpos): posBinSize: max(xpos); % center of the position bins
xposcenters = xposBins + posBinSize/2;
xposcenters(end) = [];


yposBins = min(ypos): posBinSize: max(ypos); 
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

numofStates = 40;
gamma_avg = HMM_StatesPlaceField_2D(eventsBinnedfiring, runBinnedfiring, posbinIdx, [], numofStates, xposBins, yposBins);

gamma_avg2 = permute(gamma_avg, [2 1 3]); % (1:70,1:80,:)



figure; 

runidx = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2));
for ii = 1: numofStates
    subplot(8,5, ii)
    plot(fileinfo.xyt(runidx, 1), fileinfo.xyt(runidx, 2), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
    hold on
    h1 = imagesc(xposcenters, yposcenters, gamma_avg2(:,:, ii));
    set(h1, 'AlphaData', 0.7)

    
    
    colormap('hot')
    set(gca,'YDir','normal')
    set(gca, 'XTick', [], 'YTick', [], 'box', 'off')
    xlim([min(xpos) max(xpos)])
    ylim([min(ypos) yposcenters(70)])

end


print(gcf, 'TedMaze1_postPBEs_lsPF', '-dpng')


%% Cross-validation analysis
% 
% k = 5; % number of folds
% % noEvents = 1000;
% 
% numofStates = 40;
% 
% foldSize = floor(noEvents/k);
% 
% likelihood_actual = zeros(noEvents, 1);
% likelihood_timeSwap = zeros(noEvents, 1);
% likelihood_poisson = zeros(noEvents, 1);
% likelihood_timeCycle = zeros(noEvents, 1);
% 
% for fold = 1 : k
%     
%     testSet = (fold - 1)*foldSize + 1 : fold*foldSize;
%     trainSet =  setdiff(1:noEvents, testSet);
%     
%         
%     [transmat, lambda] = trainHMM(eventsBinnedfiring(trainSet, 2), activeUnits, numofStates);
%     
%     
%     
%     % calculating the likelihoods of the test set events
%     
%     prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
% 
%     for evt = 1: length(testSet)
%         
%         % actual data
%         B = poisson_prob(eventsBinnedfiring{testSet(evt), 2}, lambda,1);
%         [~, ~, ~,  likelihood_actual(testSet(evt)), ~] = fwdback(prior, transmat, B);
%         
%         % time swap
%         B = poisson_prob(timeSwapEventsBinnedfiring{testSet(evt), 2}, lambda,1);
%         [~, ~, ~,  likelihood_timeSwap(testSet(evt)), ~] = fwdback(prior, transmat, B);
%         
%         % poisson
%         B = poisson_prob(poissonEventsBinnedfiring{testSet(evt), 2}, lambda,1);
%         [~, ~, ~,  likelihood_poisson(testSet(evt)), ~] = fwdback(prior, transmat, B);
%         
%         % timeCycle (incoherent)
%         B = poisson_prob(cycleEventsBinnedfiring{testSet(evt), 2}, lambda,1);
%         [~, ~, ~,  likelihood_timeCycle(testSet(evt)), ~] = fwdback(prior, transmat, B);
%     end
% end
%         
% figure; load('actual.mat', 'lambda_unit_sparsity')
% actualGini = lambda_unit_sparsity;
% save('post_cv.mat', 'likelihood_actual','likelihood_poisson','likelihood_timeSwap','likelihood_timeCycle')

%% HMM sequence detection

%%% In this step we test the extent to which an HMM trained on the
%%% candidate events is able to detect sequence events among non-sequence
%%% events

%%% we do this by k-fold cross-validation: using k-1 folds for training and
%%% and the remaining fold for sequence detection 


% numofFolds = 5; %% 5-fold cross-validation
% numofStates = 40;
% 
% [CVtransmat, CVlambda, testEvts] = trainCVmodels(eventsBinnedfiring, numofStates, numofFolds, fileinfo);
% 
% 
% 
% noShuffle = 1000;
% 
% dataType = 'actual';
% HMMprctile = measureHMMcongruence(eventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo);
% 
% 
% dataType = 'poisson';
% HMMprctile_p = measureHMMcongruence(poissonEventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo);
% 
% 
% dataType = 'timeswap';
% HMMprctile_ts = measureHMMcongruence(timeSwapEventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo);
% 
% 
% dataType = 'temporal';
% HMMprctile_tc = measureHMMcongruence(cycleEventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo);
% 
% 
% 
% figure;
% hold on
% 
% [h, bins] = hist(HMMprctile, 100);
% cumhist = cumsum(h);
% threshBin = find(bins > 0.95, 1, 'first');
% crosspnt = cumhist(threshBin);
% 
% h_p = hist(HMMprctile_p, bins);
% cumhist_p = cumsum(h_p);
% crosspnt_p = cumhist_p(threshBin);
% 
% 
% h_ts = hist(HMMprctile_ts, bins);
% cumhist_ts = cumsum(h_ts);
% crosspnt_ts = cumhist_ts(threshBin);
% 
% 
% h_tc = hist(HMMprctile_tc, bins);
% cumhist_tc = cumsum(h_tc);
% crosspnt_tc = cumhist_tc(threshBin);
% 
% 
% colors = [[50 150 50]/255; [150 0 0]/255; [30 144 200]/255; [180 150 110]/255];
% 
% 
% 
% p1 = plot(bins, cumhist, 'color', [50 150 50]/255, 'linewidth', 2);
% line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [50 150 50]/255,'LineStyle', '--', 'linewidth', 1)
% 
% p2 = plot(bins, cumhist_ts, 'color', [30 144 200]/255, 'linewidth', 2);
% line([0 bins(threshBin)], [crosspnt_ts crosspnt_ts], 'color', [30 144 200]/255,'LineStyle', '--', 'linewidth', 1)
% 
% p3 = plot(bins, cumhist_p, 'color', [180 150 110]/255, 'linewidth', 2);
% line([0 bins(threshBin)], [crosspnt_p crosspnt_p], 'color', [180 150 110]/255,'LineStyle', '--', 'linewidth', 1)
% 
% p4 = plot(bins, cumhist_tc, 'color', [150 0 0]/255, 'linewidth', 2);
% line([0 bins(threshBin)], [crosspnt_tc crosspnt_tc], 'color', [150 0 0]/255,'LineStyle', '--', 'linewidth', 1)
% 
% 
% line([bins(threshBin) bins(threshBin)], [0 noEvents], 'color', 'k','LineStyle', '--', 'linewidth', 1)
% 
% legend([p1,p2,p3,p4],'actual','time swap','poisson','temporal shuffle')
% 
% set(gca, 'fontsize', 14)
% xlabel('percentile score', 'fontsize', 20)
% ylabel('cumulative number of PBEs', 'fontsize', 20)
% 
% title('TedMaze1 - Quiet Wake PBEs', 'fontsize', 20)
% 
% print(gcf, 'TedMaze1_QuietWakePBEs2', '-dpng')
%     
