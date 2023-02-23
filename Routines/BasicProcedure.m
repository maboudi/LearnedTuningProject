function BasicProcedure(DatasetNum2, longORshort2, datasetDir, animalname)


DatasetNum = str2num(DatasetNum2);
longORshort = str2num(longORshort2);

los = {'long';'short'};

cd(datasetDir)

load('IIdata') % loading info data related to all sessions in the current folder


%% current session info

fileinfo = IIdata(DatasetNum);
fileinfo.animal = animalname;

try
FileBase = [datasetDir '/' fileinfo.name '/' fileinfo.name];
catch
    return
end

[Fs, lfpSampleRate, nCh] = getParameters([FileBase '.xml']);

% Add the parameters read from .xml file to fileinfo

fileinfo.Fs = Fs;
fileinfo.lfpSampleRate = lfpSampleRate;
fileinfo.nCh = nCh; 


% Load the spike data
load([FileBase '.spikeII.mat']);

%%% truncate the data in order to separate the datasets before and after
%%% track shortening (long and short)

firstLapSpikes = spike.t(find(spike.lap == 1));

[timediff, idx] = max([0; diff(firstLapSpikes)]);

breakpoint = firstLapSpikes(idx); %%% the time of first spike in the first lap after shortening the track

if timediff > max(diff(spike.t)) % to make sure breaking not happen between successive spike times in case of not shortened track

    spike = truncateSpike(spike, longORshort, breakpoint, [], [], Fs, fileinfo.pix2cm);

else
    
    breakpoint = spike.t(end);
    % there is no need anymore to truncate the data
    if longORshort == 2 % quit when there no shortening happen
	fprintf('There is no short part!\n')
    return
    end
end


%%% Basic calcultion to signify the units 

spike.unit = zeros(length(spike.t), 1); %%% defining a new field which signifies a unique unit for each cluster

currUnit = 0;
for shank = 1 : max(spike.shank)
    
    shankClusters = unique(spike.cluster(spike.shank == shank)); %%% clusters of the current shank
    for ii = 1 : length(shankClusters)
        
        currUnit = currUnit + 1;
        spkInd = find(spike.shank == shank & spike.cluster == shankClusters(ii)); 
        spike.unit(spkInd) = currUnit;

    end
    
end


%for Etienne

allIntClusters = [];

shanks = unique(spike.shank(find(spike.qclu == 5)));

for s = 1 : length(shanks)
    
    intClusters = unique(spike.cluster(find(spike.shank == shanks(s) & spike.qclu == 5)));
    
    if size(intClusters, 2) > size(intClusters, 1)
        intClusters = intClusters';
    end
        
    
    allIntClusters = [allIntClusters; [shanks(s)*ones(length(intClusters), 1) intClusters]];
    
end


save(['/Volumes/Recording_Data/Datasets/DibaBuzsaki Data' '/' fileinfo.name '-' fileinfo.animal '-' los{longORshort} '-' 'interneuronClusters.mat'], 'allIntClusters') % The name signifies that the events are from processing sdf file without further binning and filtering, etc.


%% defining the tuning of the neurons (putitive pyramidal units) 

speedThresh = 10; %% only the spikes above a 
posBinSize = 2.5; %% the size of each position bin in cm

qclus = [1 2 4 8 9]; %% pyramidal units

PCtuning(spike, qclus, fileinfo, longORshort, speedThresh, posBinSize); 


%%% loading information from Run periods; active units and their tunings

[activeUnits, tuningRL, tuningLR, RLtemplate, LRtemplate, posbincenters]= loadPFinfo(fileinfo, longORshort);

%%% smoothing the tunings 

sigma = 2;
halfwidth = 5;

tuningRL = smooth(tuningRL, sigma, halfwidth); %% RL or leftward tuning
tuningLR = smooth(tuningLR, sigma, halfwidth); %% LR or rightward tuning


%% detemining population burst periods

% Theta periods to be excluded

thetaPeriods = load([FileBase '.theta.1']);

if longORshort == 1
    idx = find(thetaPeriods(:,2)*Fs/lfpSampleRate < breakpoint);
else
    idx = find(thetaPeriods(:,1)*Fs/lfpSampleRate > breakpoint);
end

thetaPeriods = thetaPeriods(idx, :);
thetaPeriods = thetaPeriods/lfpSampleRate; % convert the periods to second

% exclude the other half of the data (if considering the long segment then exclude the short part from calculations)


if longORshort == 1
    
    sessionpart2exclude = [breakpoint/Fs (fileinfo.tend -fileinfo.tbegin)/1e6];
    
else
    sessionpart2exclude = [0.001 breakpoint/Fs];
     
end

exclude = [thetaPeriods; sessionpart2exclude];


time_resolution = 0.001; % in second
threshZ = 3; % sdf with 3 std deviation above the mean
qclus = setdiff(unique(spike.qclu), 5); % all but interneurons

[pbEvents2, spikeDensity] = PBPeriods(spike, fileinfo, qclus,  longORshort, time_resolution, threshZ, exclude);


length(pbEvents2)
%% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events

qclus = [1 2 4 5 8 9]; % 5 is the interneuron

%%% Calculating the number of firing pyramidal units within each event

noFiringUnits = zeros(size(pbEvents2, 1), 1);
for evt = 1 : size(pbEvents2, 1)
    
    evtSpkInd = find(spike.t > pbEvents2(evt, 1) & spike.t < pbEvents2(evt, 2) & ismember(spike.qclu, qclus));
    noFiringUnits(evt) = length(unique(spike.unit(evtSpkInd)));
end

%%% binning the spikes within each event

binDur = 0.02; % 20 ms bins (beside 1 ms binning for visualizing the rasters) 
eventsBinnedfiring2 = timeBinning(pbEvents2, spike, qclus, binDur, Fs);


% Remove the flanking zero-firing periods(silent bins) from the beginnig
% and end of each event. Since the interneurons were involved in calculating
% the PBEs' boundaries and here we are binning just pyramidal units' spike
% trains, some silent bins are expected. 

[eventsBinnedfiring, eventLen, eventBeg, eventEnd] = removeSideSilentBins(eventsBinnedfiring2, pbEvents2(:, 1), binDur, Fs);


% Qualify the events based on the minumum active units and length criteria
% Considering only the event with number of active units firing above 10 percent of the total or 5 (whichever is greater) and
% duration of at least 4 20ms-time bins.

idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4 & eventLen <= 25)); % maximum length 500 ms (???)

eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


% Final population burst events

PBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:end)];

noEvents = size(PBEs, 1);


% Generate output files

currDir = pwd;
Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/Population Burst Events'];
mkdir(Folderbase)

save([Folderbase '/' 'PBEvents.mat'], 'PBEs')

Filename = [Folderbase '/' fileinfo.name '.pbe.evt'];
MakeEvtFile(PBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, Fs, 1)


% Plot spike density function and the the boundaries of PBEs

figure('Visible','off');

hold on

plot((1:length(spikeDensity))*time_resolution, spikeDensity, 'linewidth', 1, 'color','k')

yL = get(gca,'YLim');

PBEs_ms = PBEs(:,1:2)*1000/Fs;

for epoch = 1 : size(PBEs_ms,1)
    
    p = patch([PBEs_ms(epoch,1) PBEs_ms(epoch,2) PBEs_ms(epoch,2) PBEs_ms(epoch,1)]*time_resolution,[-yL(2) -yL(2) yL(2) yL(2)], ...
       'r','edgecolor', 'none');
    set(p,'FaceAlpha', 0.5)
    
end

plot((1:length(spikeDensity))*time_resolution, ones(1,length(spikeDensity))*threshZ, '--r', 'linewidth', 2)

hold off

ylabel('Firing rate(Hz)', 'fontsize', 20)
xlabel('time(sec)', 'fontsize', 20)
set(gca,'fontsize', 16);

saveas(gcf, [Folderbase '/spikeDensityFunction.fig'])


%%% Measure the extent of overlap between Ripple events and PBEs %%%

% Calculate the ripple periods

% Here we intended to detect only ripples that have on average relatively large
% amplitudes over the all four shanks resided in CA1. Alternatively, we can
% include also ripple events (with usually smaller amplitude) picked by
% just one or a subset of shanks. 

best_channels = BigRippleChannels(fileinfo,1);
fileinfo.bestch = best_channels;

rippPowThresh = 3; 
rippleEvents = RippleDetect(fileinfo, rippPowThresh); % in lfpSampleRate


% Save the results

Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/RippleEvents'];
mkdir(Folderbase)

save([Folderbase '/' 'Ripples_thresh' num2str(rippPowThresh) '.mat'], 'rippleEvents')

Filename = [Folderbase '/' fileinfo.name '.rip.evt'];
MakeEvtFile(rippleEvents(:,1:3), Filename,{'str','peak','end'},lfpSampleRate); % make evt file for neuroscope browsing

% PBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:end)];

noEvents = size(PBEs, 1);


% Generate output files

currDir = pwd;
Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/Population Burst Events'];
mkdir(Folderbase)

save([Folderbase '/' 'PBEvents.mat'], 'PBEs')

Filename = [Folderbase '/' fileinfo.name '.pbe.evt'];
MakeEvtFile(PBEs(:, 1:3), Filename, {'beg', 'end', 'peak'}, Fs, 1)


% Plot spike density function and the the boundaries of PBEs

figure('Visible','off');

hold on

plot((1:length(spikeDensity))*time_resolution, spikeDensity, 'linewidth', 1, 'color','k')

yL = get(gca,'YLim');

PBEs_ms = PBEs(:,1:2)*1000/Fs;

for epoch = 1 : size(PBEs_ms,1)
    
    p = patch([PBEs_ms(epoch,1) PBEs_ms(epoch,2) PBEs_ms(epoch,2) PBEs_ms(epoch,1)]*time_resolution,[-yL(2) -yL(2) yL(2) yL(2)], ...
       'r','edgecolor', 'none');
    set(p,'FaceAlpha', 0.5)
    
end

plot((1:length(spikeDensity))*time_resolution, ones(1,length(spikeDensity))*threshZ, '--r', 'linewidth', 2)

hold off

ylabel('Firing rate(Hz)', 'fontsize', 20)
xlabel('time(sec)', 'fontsize', 20)
set(gca,'fontsize', 16);

saveas(gcf, [Folderbase '/spikeDensityFunction.fig'])


% Find the overlap 

rplOverlap = zeros(noEvents, 1);
noRipples = size(rippleEvents, 1);

Ripples = rippleEvents/lfpSampleRate; % in sec
PBEvents = PBEs/Fs; % in sec

for evt = 1 : noEvents
    currEvt = PBEvents(evt,:);
    for rip = 1:noRipples
        
        if (currEvt(1) >= Ripples(rip, 1) && currEvt(1) <= Ripples(rip, 2)) || ...
                (currEvt(2) >= Ripples(rip, 1) && currEvt(2) <= Ripples(rip, 2))...
                || (currEvt(1) <= Ripples(rip, 1) && currEvt(2) >= Ripples(rip, 2))
            
            rplOverlap(evt) = 1;
            break
        end
    end
end

PopBurst_Ripple_overlap = length(find(rplOverlap))/noEvents; %#ok<NASGU>

Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/Population Burst Events'];
save([Folderbase '/' 'PopBurst_Ripple_overlap.mat'], 'PopBurst_Ripple_overlap')


%%% Visulaize the PBEs and ripple events in regard to track positions %%%

% Position

tbegin = fileinfo.tbegin;
tend = fileinfo.tend;
xyt = fileinfo.xyt;
pix2cm = fileinfo.pix2cm;

withinRange = find(xyt(:,3) >= tbegin & xyt(:,3) <= tend); % selecting only the position samples occured during neural recording
xyt = xyt(withinRange, :);

timepnts = (xyt(:,3) - tbegin)/10^6; % u seconds to seconds

xpos = xyt(:,1) * pix2cm; % Converting the position data (from pixel number to centimeter)
ypos = xyt(:,2) * pix2cm;

% Speed

diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];
difftt = [1; diff(timepnts)];

% Calculation of speed in just x direction

speed = abs(diffx)./difftt; %% speed before smoothing
% speed = sqrt(diffx.^2 + diffy.^2)./difftt; 

sigma = 15; %%% smoothing the speed, the positon sampling rate is around 30 Hz, so duration of the sigma is about 500 ms
halfwidth = 3 * sigma;

smoothwin = gausswindow(sigma, halfwidth);
speed = conv(speed, smoothwin, 'same');

% PBEs' and ripples's position and speed on the track

popBurstPos = interp1(timepnts, xpos, PBEvents(:,3));
ripplePos = interp1(timepnts, xpos, Ripples(:, 2)); % peak of the ripples

popBurstSpeed = interp1(timepnts, speed, PBEvents(:,3));


% Plot the events (PBEs and ripples) in relation to the track positions

figure('Visible','off');
subplot(5,1,1:3) %%% spike times as a function of position

spikeSpeed = interp1(timepnts, speed, spike.t/Fs);
spikePos = interp1(timepnts, xpos, spike.t/Fs); 

plot(spikePos, spike.t/Fs, '.','color', 'k','markersize', 0.5)
hold on
plot(spikePos(spikeSpeed < 5), spike.t((spikeSpeed < 5))/Fs, '.','color', [0 191 255]/255,'markersize', 0.5)

ylabel('Time(Sec)', 'fontsize', 16)
xrange = xlim;

h1 = subplot(5,1,4);
[n1, pb] = hist(ripplePos, 50);
bar(pb, n1, 'FaceColor', 'k', 'EdgeColor', 'none')
set(h1,'xlim', xrange)
ylabel('#Ripples', 'fontsize', 12)

h2 = subplot(5,1,5);
[n2, pb]= hist(popBurstPos, 50);
bar(pb, n2, 'FaceColor', 'k', 'EdgeColor', 'none')
set(h2, 'xlim',xrange)
xlabel('Position(cm)', 'fontsize', 16)
ylabel('#PopBurst', 'fontsize', 12)
set(gca,'fontsize', 16);

axis([h1 h2],[xrange 0 max([n1 n2])])

Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/Population Burst Events'];
saveas(gcf, [Folderbase '/' 'populationBusrtVSripplePosition.fig'])


% Lets find out what is speed distribution of the ripple-overlapping and -nonoveralpping PBEs 
% We may find non-overlapping events with speed above 5 cm/s but would be
% negligibel if there is just a few of them

[~, bins] = hist(popBurstSpeed, 50);

h1 = hist(popBurstSpeed(find(rplOverlap)), bins);
h2 = hist(popBurstSpeed(find(~rplOverlap)), bins);


figure;
hold on

bar(bins, h1, 'facecolor', 'none', 'edgecolor', 'g') % ovelapping
bar(bins, h2, 'facecolor', 'none', 'edgecolor', 'r') % non-overlapping

xlabel('Speed(cm/s)', 'fontsize', 20)
ylabel('Number of events', 'fontsize', 20)
title('Speed histogram of PBEs (overlap vs non-overlap)', 'fontsize', 20)
set(gca, 'fontsize', 16)

saveas(gcf, [Folderbase '/' 'popBurstsSpeed.fig'])


%% charachterizing firings of the neurons within each event, and how the neurons particpated in the events

unitsParticipation2(eventsBinnedfiring, activeUnits, RLtemplate, LRtemplate, fileinfo, longORshort,  binDur);


%% binning run data


binDur = 0.25; % in sec 

runBinnedfiring2 = timeBinning(thetaPeriods .* Fs/lfpSampleRate, spike, qclus, binDur, Fs); % the theta periods

runBinnedfiring{1,1} = cell2mat(runBinnedfiring2(:,1)'); % concatenating all the run binned spikes
runBinnedfiring{1,2} = cell2mat(runBinnedfiring2(:,2)'); % the same as the first column but with milisecond bins


% Calculate the track positions corresponding to each time bin

runBinPos = cell(size(thetaPeriods, 1), 1);
periods = thetaPeriods ./ lfpSampleRate;
for jj = 1 : size(periods, 1)
    
    runBinPos{jj} = interp1(timepnts, xpos, (periods(jj, 1):binDur:periods(jj, 2))+binDur/2);
    runBinPos{jj}(end) = [];
end
    
runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields

noPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
positionBins = min(xpos): posBinSize :max(xpos); % center of the position bins

[~, posbinIdx] = histc(runBinPos, positionBins);



%%% Calculate the track's profile of occupancy time  %%%%
% (we can add this to the PCtuning function?)


ii_ok = zeros(length(xpos), 1);

for jj = 1 : size(thetaPeriods, 1)
    
    excludePer = thetaPeriods(jj,:) ./ lfpSampleRate;
    
    ii_ok(find(timepnts > excludePer(1) & timepnts < excludePer(2))) = 1;
end

xpos_ok = xpos(find(ii_ok));


nSamples_perPosbin = histc(xpos_ok, positionBins);
nSamples_perPosbin(end) = [];
posbin_Occ = nSamples_perPosbin*mean(diff(timepnts));

figure; 
x0=0;
y0=0;
width=400;
height=150;
set(gcf,'units','points','position',[x0,y0,width,height])

plot(posbincenters, posbin_Occ, 'k', 'linewidth', 2)

set(gca, 'color', 'none', 'box', 'off')

difpos = posbincenters(2)- posbincenters(1);
xlim([posbincenters(1)-difpos/2 posbincenters(end)+difpos/2])

xlabel('Position(cm)', 'fontsize', 10)
ylabel('Occupancy(sec)', 'fontsize', 10)


pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])

Folderbase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/PlaceFields'];
saveas(gcf, [Folderbase '/occupancy_time.fig'])
    


%% Generating shuffle datasets

% generate poisson simulated dataset

binDur = 0.02;

poissonEventsBinnedfiring = poissonSpikeTrain(eventsBinnedfiring, binDur); 
cmprPoisson2Raw(poissonEventsBinnedfiring, eventsBinnedfiring, longORshort, fileinfo, binDur); %%% compare the poisson simulated data with the actual data distribution

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



% Spike identity shuffle

spikeIDEventsBinnedfiring = cell(noEvents, 2);

for evt = 1 : noEvents
    
    currEvt = eventsBinnedfiring(evt, :);
    noUnits = size(currEvt{1}, 1);
    noTimebins = size(currEvt{1}, 2);
    
    shuffle = zeros(size(currEvt{1})); % 1ms time bins
    for ii = 1 : noTimebins
        randIdx = randperm(noUnits);
        shuffle(:, ii) = currEvt{1}(randIdx, ii); % 1 ms bin
    end  
    
    spikeIDEventsBinnedfiring{evt,1} = shuffle;
    spikeIDEventsBinnedfiring{evt,2} = permute(sum(reshape(shuffle, [noUnits, 20, size(shuffle,2)/20]), 2), [1,3,2]);

end


% Unit idenitity shuffle
% unit ID shuffles done independently for each event

unitIDEventsBinnedfiring = cell(noEvents, 2);

for evt = 1 : noEvents
    randIDs = randperm(noUnits);
    
    for binning = 1:2
        unitIDEventsBinnedfiring{evt, binning} = eventsBinnedfiring{evt, binning}(randIDs, :);
    end
end


%% Replay(sequence) Detection (Using Bayesian Decoding, and Hidden Markov Model)
% 
% 
% %% Bayesian Decoding 
% 
% dataType = 'actual';
% [BDScoreZ, mcPval, PosteriorProbMatrix, begPosition, endPosition] = BDreplayDetect(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur, longORshort, fileinfo, dataType);
% 
% % The last dimension of the output cell arrays and matrices
% % signifies different run directions. Note that for decoding the positions
% % in each time bin, we can use either of the tunings corresponding to each
% % direction, i.e RL and LR. We can use the the probability distributions themselves or
% % calculate a third distribution by summing across them. 
% 
% 
% %% replay detection in surrogate datasets (Poisson simulated data)
% 
% 
% %%% Bayesian decoding replay detection for poisson simulated data
% 
% dataType = 'poisson';
% [BDScoreZ_p, mcPval_p, PosteriorProbMatrix_p, begPosition_p, endPosition_p] = BDreplayDetect(poissonEventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur, longORshort, fileinfo, dataType);
% 
% 
% %% Time swap
% 
% 
% %%% Bayesian decoding replay detection for time swap dataset
% 
% dataType = 'TimeSwap';
% [BDScoreZ_ts, mcPval_ts, PosteriorProbMatrix_ts, begPosition_ts, endPosition_ts] = BDreplayDetect(timeSwapEventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur, longORshort, fileinfo, dataType);
% 
% 

%% Hidden Markov modeling (HMM)

% 
% % What is the optimal number of states 
% 
% 
% % (1)cross-validation
% 
% currDir = pwd;
% FileBase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/HMM/Optimal number of States/cross_validation' ];
% mkdir(FileBase)
% 
% Kfold = 5;
% maxNumStates = 60;
% stepSize = 5;
% HMM_numberStates(eventsBinnedfiring, activeUnits, Kfold, maxNumStates, stepSize, FileBase);
% 
% 
% % (2)states' (virtual) place fields
% 
% FileBase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/HMM/Optimal number of States/States_PlaceField' ];
% mkdir(FileBase)
% 
% dims = 10:10:60; % the set of number of states
% 
% SpatialInfo_uniOccProb = cell(length(dims), 1);
% SpatialInfo_actOccProb = cell(length(dims), 1);
% 
% for ii = 1:length(dims)
%     
%     numofStates = dims(ii);
%     
%     [SpatialInfo_uniOccProb{ii}, SpatialInfo_actOccProb{ii}] = ...
%         HMM_StatesPlaceField(eventsBinnedfiring, runBinnedfiring, ...
%         posbinIdx, posbin_Occ, activeUnits, numofStates, posbincenters, Filebase);
% end
% 
% % plot boxplots each for a specified number of states
% 
% group = [];
% for ii = 1:lenght(dims)
%     group = [group; ii*ones(dims(ii), 1)];
% end
% 

% figure('visible', 'off');
% 
% x0=0;
% y0=0;
% width=300;
% height=300;
% set(gcf,'units','points','position',[x0,y0,width,height])

% plotbox(group, cell2mat(SpatialInfo_uniOccProb), dims) % Spatial Information(Uni)
% xlabel('Number of states', 'fontsize', 14)
% ylabel('Spatial Information(Uni)', 'fontsize', 14)
% set(gca, 'fontsize', 14)

% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
% print(gcf,[FileBase '/' 'Spatial Information(Uni)'],'-dpdf','-r0')


% figure('visible', 'off');
% set(gcf,'units','points','position',[x0,y0,width,height])

% plotbox(group, cell2mat(SpatialInfo_actOccProb), dims) % Spatial Information(Act)
% xlabel('Number of states', 'fontsize', 14)
% ylabel('Spatial Information(Act)', 'fontsize', 14)
% set(gca, 'fontsize', 14)

% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
% print(gcf,[FileBase '/' 'Spatial Information(Act)'],'-dpdf','-r0')



%% Analyzing the hidden states regarding sparsity, and state place fields

FileBase = [datasetDir '/' fileinfo.name '/data part' num2str(longORshort) '/HMM/sparsity'];
mkdir(FileBase)


numofStates = 30; % we have fixed the number of states to 30 here
numIter = 50; % because different initializations lead to differnt results we do the training several times and calculate the distribution


% actual data
dataType = 'actual'; 
modelSparsity(eventsBinnedfiring(1:1000, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)


%%%% train HMMs on surrogate datasets %%%

% (1)poisson simulated spike trains
dataType = 'poisson';
modelSparsity(poissonEventsBinnedfiring(1:1000, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)

% % (2)incoherent time bin shuffle
% dataType = 'timeCycle';
% modelSparsity(cycleEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)

% (3)time swap
dataType = 'timeswap';
modelSparsity(timeSwapEventsBinnedfiring(1:1000, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)

% % (4)spike identity shuffle
% dataType = 'SpikeID';
% modelSparsity(spikeIDEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)
% 
% % (5)unit identity shuffle
% dataType = 'UnitID';
% modelSparsity(unitIDEventsBinnedfiring(:, 2), activeUnits, numofStates, [FileBase '/' dataType], numIter)


% States place fields for actual data and position-shuffled data

FileBase = [datasetDir '/' fileinfo.name '/data part' num2str(longORshort) '/HMM/state_placeFields'];
mkdir(FileBase)

numofStates = 30;

dir1 = [FileBase '/position_intact'];
dir2 = [FileBase '/position_shuffled'];

mkdir(dir1)
mkdir(dir2)

[actual_SIuni, actual_SIact] = ...
        HMM_StatesPlaceField(eventsBinnedfiring, runBinnedfiring, ...
        posbinIdx, posbin_Occ, activeUnits, numofStates, posbincenters, dir1);
  
posbinIdx_shuffle = posbinIdx(randperm(length(posbinIdx)));

[shuffle_SIuni, shuffle_SIact] = ...
        HMM_StatesPlaceField(eventsBinnedfiring, runBinnedfiring, ...
        posbinIdx_shuffle, posbin_Occ, activeUnits, numofStates, posbincenters, dir2);
    

% plot boxplots each for a specified number of states

group = [];
for ii = 1:2
    group = [group; ii*ones(30, 1)];
end



figure('visible', 'off');

x0=0;
y0=0;
width=300;
height=300;
set(gcf,'units','points','position',[x0,y0,width,height])

plotbox(group, [actual_SIuni; shuffle_SIuni] , [1 2]) % Spatial Information(Uni)

xlabel('Number of states', 'fontsize', 14)
ylabel('Spatial Information(Uni)', 'fontsize', 14)
set(gca, 'fontsize', 14)

pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(gcf,[FileBase '/' 'Spatial Information(Uni)'],'-dpdf','-r0')





figure('visible', 'off');
set(gcf,'units','points','position',[x0,y0,width,height])

plotbox(group, [actual_SIact; shuffle_SIact] , [1 2]) % Spatial Information(Act)
xlabel('Number of states', 'fontsize', 14)
ylabel('Spatial Information(Act)', 'fontsize', 14)
set(gca, 'fontsize', 14)

pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)])
print(gcf,[FileBase '/' 'Spatial Information(Act)'],'-dpdf','-r0')


% 
% %% HMM: training scale effect
% % Determining the asymptotic amount of data needed to approximate an
% % optimal model.
% % Using cross validation, we divide the data to two halves, for training
% % and test. Then we train models each time based on a different amount of
% % training data and then evaluate the performance by measuring
% % log Likelihood in the validation set.
% 
% if noEvents > 300
% 
%     HMMScaling(eventsBinnedfiring, activeUnits, fileinfo, longORshort);
%     
% end
% 
% 
% %% HMM sequence detection
% 
% %%% In this step we test the extent to which an HMM trained on the
% %%% candidate events is able to detect sequence events among non-sequence
% %%% events
% 
% %%% we do this by k-fold cross-validation: using k-1 folds for training and
% %%% and the remaining fold for sequence detection 
% 
% 
% numofFolds = 5; %% 5-fold cross-validation
% 
% [CVtransmat, CVlambda, testEvts] = trainCVmodels(eventsBinnedfiring, tuningRL, tuningLR, positionBins, numofStates, numofFolds, binDur, fileinfo, longORshort);
% 
% 
% 
% noShuffle = 1000;
% 
% dataType = 'actual';
% [HMMScoreZ, binsCongSignificance, maxCongruenceLength] = measureHMMcongruence(eventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo, longORshort);
% 
% 
% dataType = 'poisson';
% [HMMScoreZ_p, binsCongSignificance_p, maxCongruenceLength_p] = measureHMMcongruence(poissonEventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo, longORshort);
% 
% 
% dataType = 'timeswap';
% [HMMScoreZ_ts, binsCongSignificance_ts, maxCongruenceLength_ts] = measureHMMcongruence(timeSwapEventsBinnedfiring, testEvts, CVtransmat, CVlambda, noShuffle, dataType, fileinfo, longORshort);
% 
% 
% 
% 
% 
% 
% %% Visualization of the events and comparing between BD and HMM in detection of sequences
% 
% FileBase = [datasetDir '/' fileinfo.name '/data part' num2str(longORshort) '/methodsComparison'];
% 
% 
% dataType = 'actual';
% [confMat, confMatIndices] = calculateConfMat(HMMScoreZ, BDScoreZ(:,1,3), noEvents, FileBase, dataType); %% confusion Matrix
% 
% DetectionResultsVisualization(eventsBinnedfiring, confMatIndices, PosteriorProbMatrix, begPosition, endPosition, ...
%         binsCongSignificance, RLtemplate, LRtemplate, binDur, FileBase, dataType); %% visualization of the events
% 
% 
% 
% dataType = 'Poisson';
% [confMat_p, confMatIndices_p] = calculateConfMat(HMMScoreZ_p, BDScoreZ_p(:,1,3), noEvents, FileBase, dataType);
% 
% DetectionResultsVisualization(poissonEventsBinnedfiring, confMatIndices_p, PosteriorProbMatrix_p, begPosition_p, endPosition_p, ...
%         binsCongSignificance_p, RLtemplate, LRtemplate, binDur, FileBase, dataType); 
% 
% 
% 
% dataType = 'timeSwap';
% [confMat_ts, confMatIndices_ts] = calculateConfMat(HMMScoreZ_ts, BDScoreZ_ts(:,1,3), noEvents, FileBase, dataType);
% 
% DetectionResultsVisualization(timeSwapEventsBinnedfiring, confMatIndices_ts, PosteriorProbMatrix_ts, begPosition_ts, endPosition_ts, ...
%         binsCongSignificance_ts, RLtemplate, LRtemplate, binDur, FileBase, dataType); 
%     
%     
% 
% end
% 
% 
