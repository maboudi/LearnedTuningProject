clear; clc

cd('/home/kouroshmaboudi/Documents/LaurelData/ObjectLocationMemory/Charlie')

% Dataset directory
currDir = pwd;

%% session info

animal = 'charlie';
sessionName = 'PlacePreference-Charlie-20150624-01-nlx';

fileinfo = struct('name', sessionName, 'animal', animal, 'xyt', [], 'tbegin', [], 'tend', ...
                    [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'CA1thetach', [], 'CA', [1 1 1 1], 'pos2cm', 240); 


% (1)  Sleep 1: ~20 minute baseline sleep in homecage.
% (2)  Exposure 1: exposure to empty box (~6 min).
% (3)  Homecage 1: return to home cage (~3 min).
% (4)  Exposure 2: exposure to box with objects (~6 min).
% (5)  Homecage 2: return to homecage (~3 min). 
% (6)  Exposure 3: exposure to box with objects (~6 min).
% (7)  Homecage 3: return to homecage (~3 min). 
% (8)  Exposure 4: exposure to box with objects (~6 min). 
% (9)  Sleep 2: extended sleep (~4 hours). 
% (10) Memory recall test: exposure to box with displaced object (~6 min).
% (11) rest of the reocring for a short period (~20 min)


FileBase = [currDir '/' fileinfo.name '/' fileinfo.name];

% Loading the mat fiels
% varList = {'basics'; 'speed'; 'spikes'; 'hiAmpSpindles'};
varList = {'-basics'; '-speed'; '-spikes'; '.SleepState.states'};

for var = 1:length(varList)
    load([FileBase  varList{var} '.mat'])
end


% fileinfo.Fs = basics.SampleRate;

fileinfo.Fs = 1e6; % we don't use the original sample rate, because the time stamps were converted to micro sec. 
fileinfo.lfpSampleRate = basics.lfpSampleRate; 
fileinfo.nCh = basics.nChannels;

fileinfo.xyt = [speed.x*fileinfo.pos2cm speed.y*fileinfo.pos2cm speed.t]; 


% periods
periods = basics.period.time;


%% processing spikes

nUnits = length(spikes);

for unit = 1:length(spikes)
    
    spikes(unit).rate = length(spikes(unit).time)/((periods(end, 2) - periods(1,1))/1e6);

end



spike = struct('t', [], 'unit', [], 'shank', [], 'qclu', [], 'x', [], 'y', [], 'speed', [], 'rate', []);

for unit = 1: nUnits
    
    if spikes(unit).quality == 0
        continue
    end
    
    unitSpikes = spikes(unit).time;
    
    spike.t = [spike.t; unitSpikes];
    
    spike.unit = [spike.unit; repelem(unit, length(unitSpikes))'];
    spike.shank = [spike.shank; repelem(spikes(unit).id(1), length(unitSpikes))'];
    spike.qclu = [spike.qclu; repelem(spikes(unit).quality, length(unitSpikes))'];
    spike.rate = [spike.rate; repelem(spikes(unit).rate, length(unitSpikes))'];
    
    spike.x = [spike.x; spikes(unit).x*fileinfo.pos2cm];
    spike.y = [spike.y; spikes(unit).y*fileinfo.pos2cm];
    spike.speed = [spike.speed; spikes(unit).speed];
    
    
end


% sorting the spike time stamps

[~, spikeSortIdx] = sort(spike.t, 'ascend');
spike = structfun(@(x)(x(spikeSortIdx)), spike,'UniformOutput',false);


% spike.qclu(spike.rate > 2 & ismember(spike.qclu, [1 3])) = 2; % remove the high firing pyramidal units


%% PBEs

    %% detemining population burst periods

    
% exclude all the periods with speed above 5 cm/s

velocity = speed.v;

% smoothing the speed; the sampling rate is 100 Hz (after interpolation).

sigma = 100; % 1 sec
halfwidth = 3*sigma;
win = gausswindow(sigma, halfwidth); % smoothing kernel

velocity = conv(velocity, win, 'same');  

% exclude theta, or in case of Charlie (day 1) periods with high amplitude
% spindles

exclude = []; 

% exclude = [exclude HiAmpSpindles/fileinfo.lfpSampleRate*fileinfo.Fs]; % high volatge spindles


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the meanposbinIdx(:,1)



% memory test period

pers = 9;
combinedPBEs = [];

for ii = 1:length(pers)
    
fileinfo.tbegin = periods(pers(ii),1); 
fileinfo.tend = periods(pers(ii),2);


[primaryPBEs, sdat] = PBPeriods(spike, fileinfo, [], time_resolution, threshZ, exclude);
eventVelocity = interp1(speed.t, velocity, primaryPBEs(:,3));
primaryPBEs2 = primaryPBEs(eventVelocity < 5, :); % the exclusion criterion here is based on speed, ideally should be based on the theta
                                                  % some high attentive states could be accompanied with high amplitude theta oscillations 
                                            

combinedPBEs = [combinedPBEs; primaryPBEs2];  

end
                                            
                                                  
noPrimEvents = size(combinedPBEs, 1);



% Binning and assessing the duration and participation ratio for each
% event; excluding the events with duration less than 4 time bins(80 ms)
% and participation of less than 10% of active units


% remove all of the non-pyraidal units

keepIdx = find(ismember(spike.qclu, [1 3]));
spike = structfun(@(x)(x(keepIdx)), spike,'UniformOutput',false);

allActiveUnits = unique(spike.unit);
[~, sortIdx] = sort(allActiveUnits, 'ascend');


temp = spike.unit;

for ii = 1: length(sortIdx)
    temp(find(spike.unit == allActiveUnits(ii))) = sortIdx(ii);
end

spike.unit = temp;


[eventsBinnedfiring, secondaryPBEs, acceptedprmEvts] = finalBinningResult(combinedPBEs, spike, fileinfo);



savePBEs(combinedPBEs, secondaryPBEs, acceptedprmEvts, 'test', fileinfo, periods(1,1))

% [poissonEventsBinnedfiring, timeSwapEventsBinnedfiring, temporalEventsBinnedfiring, pooledtimeswapEventsBinnedfiring] = genSurrogates(eventsBinnedfiring);

noSecEvents = size(secondaryPBEs, 1);


%% population burst events during the exposure peridos


% memory test period

pers = [4 6 8];
excludePBEs = [];

for ii = 1:length(pers)
    
fileinfo.tbegin = periods(pers(ii),1); 
fileinfo.tend = periods(pers(ii),2);


[primaryPBEs, sdat] = PBPeriods(spike, fileinfo, [], time_resolution, threshZ, exclude);
eventVelocity = interp1(speed.t, velocity, primaryPBEs(:,3));
primaryPBEs2 = primaryPBEs(eventVelocity < 5, :); % the exclusion criterion here is based on speed, ideally should be based on the theta
                                                  % some high attentive states could be accompanied with high amplitude theta oscillations 
                                            

excludePBEs = [excludePBEs; primaryPBEs2];  

end


%% Ripple detection


% Here we intended to detect only ripples that have on average relatively large
% amplitudes over the all four shanks resided in CA1. Alternatively, we can
% include also ripple events (with usually smaller amplitude) picked by
% just one or a subset of shanks. 

best_channels = BigRippleChannels(fileinfo,1);
fileinfo.bestch = best_channels;

rippPowThresh = 3; 
rippleEvents = RippleDetect(fileinfo, rippPowThresh); % in lfpSampleRate (beg peak end normalized amplitude)

Filename = [FileBase '.spl.evt'];
MakeEvtFile(rippleEvents(:, 1:3), Filename, {'beg', 'peak', 'end'}, fileinfo.lfpSampleRate, 1)


eventVelocity = interp1(speed.t, velocity, (rippleEvents(:,2)/fileinfo.lfpSampleRate *1e6) + periods(1,1));
rippleEvents2 = rippleEvents(eventVelocity < 5, :); % refining the ripples events; excluding the ones that happen during mobile state


%% Detected PBEs overlap with the ripples


rplOverlap = zeros(noSecEvents, 1);

noRipples = size(rippleEvents2, 1);

Ripples = rippleEvents2;
Ripples(:, 1:3) = (rippleEvents2(:, 1:3)/fileinfo.lfpSampleRate) + periods(1,1)/fileinfo.Fs; % in sec
PBEvents = secondaryPBEs/fileinfo.Fs; % in sec

for evt = 1 : noSecEvents
    currEvt = PBEvents(evt,:);
    for rip = 1:noRipples
        
        if (currEvt(1) >= Ripples(rip, 1) && currEvt(1) <= Ripples(rip, 3)) || ...
                (currEvt(2) >= Ripples(rip, 1) && currEvt(2) <= Ripples(rip, 3))...
                || (currEvt(1) <= Ripples(rip, 1) && currEvt(2) >= Ripples(rip, 3))
            
            rplOverlap(evt) = 1;
            break
        end
    end
end

PopBurst_Ripple_overlap = length(find(rplOverlap))/length(PBEvents)



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

pbePeriod = 'extended sleep';

FileBase = [currDir '/' fileinfo.name  '/HMM/lsPFs/'];
mkdir(FileBase)


binDur = 0.1; % bin duration of the behavioral data (I need to test with different bin durations)
numofStates = 80;
posBinSize = 1; 


% % active run periods
% 
% velocityIdx = velocity > 10;
% 
% tran2HiVelPnts = find([0; diff(velocityIdx)] == 1);
% tran2LoVelPnts = find([0; diff(velocityIdx)] == -1);
% 
% 
% if length(tran2HiVelPnts) == length(tran2LoVelPnts)+1
%     
%     tran2LoVelPnts = [tran2LoVelPnts; find(velocityIdx, 1, 'last')];
%     
% elseif length(tran2LoVelPnts) == length(tran2HiVelPnts)+1
% 
%     tran2HiVelPnts = [find(velocityIdx, 1, 'first'); tran2HiVelPnts];
%     
% elseif tran2LoVelPnts(1) < tran2HiVelPnts(1)
% 
%     tran2LoVelPnts = [tran2LoVelPnts; find(velocityIdx, 1, 'last')];
%     tran2HiVelPnts = [find(velocityIdx, 1, 'first'); tran2HiVelPnts];
%     
% end
% 
% activePeriods = [tran2HiVelPnts tran2LoVelPnts];
% activePeriods = activePeriods/100*1e6 + periods(1,1);
% 
% 
% 
% % Active periods that happen during the exposure of the animal to the
% % objects
% 
% 
% boxActivePeriods = activePeriods(activePeriods(:,1) > periods(4,1) & activePeriods(:,2) < periods(8,2), :);



runPeriods = periods([4 6 8], :);

SecRunPeriods = [];
for ii = 1: size(runPeriods, 1)
    
    excludes = excludePBEs(excludePBEs(:,1) > runPeriods(ii,1) & excludePBEs(:,2) < runPeriods(ii,2), :);
     
   
    if ~isempty(excludes)
        SecRunPeriods = [SecRunPeriods; runPeriods(ii,1) excludes(1,1)];

        for jj = 1: size(excludes, 1)-1
            SecRunPeriods = [SecRunPeriods; excludes(jj,2) excludes(jj+1, 1)];
        end

        SecRunPeriods = [SecRunPeriods; excludes(end, 2) runPeriods(ii,2)];
    else
        SecRunPeriods = [SecRunPeriods; runPeriods(ii,1) runPeriods(ii,2)];
    end
    
end



runDuration = (SecRunPeriods(:, 2) - SecRunPeriods(:, 1))/1e6;
% initialPosition = interp1(fileinfo.xyt(:, 3), fileinfo.xyt(:, 1), SecRunPeriods(:, 1));


boxThreshold = 70;
SecRunPeriods = SecRunPeriods(runDuration > binDur, :); %   



% Binnig the active run period
[runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(spike, SecRunPeriods, periods, binDur, posBinSize, fileinfo);



insideTheboxBins = find(posbinIdx(:, 1) > boxThreshold/posBinSize);
posbinIdx = posbinIdx(insideTheboxBins, :);


runBinnedfiring{1} = [];
runBinnedfiring{2} = runBinnedfiring{2}(:, insideTheboxBins);



% Train an HMM on PBEs

activeUnits = 1:size(eventsBinnedfiring{1,2}, 1);
[transmat, lambda] = trainHMM(eventsBinnedfiring(:,2), activeUnits, numofStates);

% 
% [lsPFs, max_prob] = plot_lsPF(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, periods, pbePeriod);
% 
% set(gcf, 'units', 'points')
% set(gcf, 'position', [341.8105   49.2632  730.6105  930.6947])
% print(gcf, [fileinfo.name '_lsPFs'],'-dpng','-r0')

% savepdf(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'])
% print(gcf, [FileBase fileinfo.name '_' pbePeriod '_lasPFs'], '-dpng')
% 
% figure; 
% [h, bins] = hist(max_prob, 20);
% bar(bins, h, 'EdgeColor', 'k', 'FaceColor', 'k')
% 
% set(gca, 'box', 'off', 'linewidth', 3, 'fontsize', 14)
% xlabel('activation probability', 'fontsize', 20)
% ylabel('number of states', 'fontsize', 20)
% 
% 
% set(gcf, 'units', 'points')
% pos = get(gcf, 'position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf, [fileinfo.name '_lsPFs'],'-dpng','-r0')
% 


[PosDecError, PosDecError_shuffledOrder, PosDecError_shuffledGamma] = positionDecodingError(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates);

figure;
set(gcf, 'Units', 'points', 'position', [539   375   560   420])


plotErrorCdf(PosDecError, PosDecError_shuffledOrder, PosDecError_shuffledGamma, 1)
% title('Decoding error', 'fontsize', 16)

set(gcf, 'units', 'points')
pos = get(gcf, 'position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf, [fileinfo.name '_decodingError'],'-dpdf','-r0')




%% HMM sequence detection
% 
% % In this step we test the extent to which an HMM trained on the
% % candidate events is able to detect sequence events among non-sequence
% % events
% % 
% % we do this by k-fold cross-validation: using k-1 folds for training and
% % and the remaining fold for sequence detection %% this is temporary assigning. Later will change it when analyzing run/post-sleep data
% 
% % first of all randomize the events
% 
% 
% 
% numofFolds = 5; %% 5-fold cross-validation
% numofStates = 40;
% noShuffle = 5000;
% 
% 
% acc = [];
% acc_pts = [];
% 
% for ii = 1:1
%     
%     ii
%     
% eventsBinnedfiring = eventsBinnedfiring(randperm(size(eventsBinnedfiring, 1)), :);
% 
% % test
% 
% pbePeriod = 'combined_exposures+interleavedrest';
% 
% [CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiring, numofStates, numofFolds, fileinfo, pbePeriod); 
% 
% HMMprctile2 = HMMcongruence_crossvalid(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
% % HMMprctile_p = HMMcongruence_crossvalid(poissonEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
% % HMMprctile_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
% % HMMprctile_t = HMMcongruence_crossvalid(temporalEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);
% 
% 
% % Pooled time swap with replacement and resampling the bins from all the
% % events not only the test set; we need to replace this surrogate with the
% % one with sampling with replacement and sampling from only test set
% % HMMprctile_pts2 = HMMcongruence_crossvalid(pooledtimeswapEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap_pooled', fileinfo, pbePeriod);
% 
% HMMprctile_pts2 = HMMcongruence_crossvalid_pts(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);
% % 
% % % HMMprctile_PBEshuffle = HMMcongruence_crossvalid_PBEshuffle(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual');
% % 
% 
% acc = [acc; HMMprctile2];
% acc_pts = [acc_pts; HMMprctile_pts2];
% end
% 
% 
% %% plot
% 
% figure;
% set(gcf, 'position', [874  430 393 255])
% 
% % cross-validations (within-periods)
% 
% % congruencePlot(HMMprctile, HMMprctile_p, HMMprctile_ts, HMMprctile_t, HMMprctile_pts, 1, 'Extended sleep (4 hr) Cross-Validation')
% 
% 
% % flip the curve like the one in the paper (something like 1 minus the curve in the top figure generated by the function congruencePlot)
% % congruencePlot2(HMMprctile2, HMMprctile_pts2, 1, 'combined exposures and rest Cross-Validation')
% congruencePlot2(acc, acc_pts, 1, 'combined exposures and rest Cross-Validation')
% 
% % 
% % figure;
% % idx = find(speed.t > fileinfo.tbegin & speed.t < fileinfo.tend);
% % plot(speed.x(idx), speed.y(idx), 'linewidth', 1)
% % 
% % hold on
% % 
% % for ii = 1:length(secondaryPBEs)
% %     plot(eventX(ii), eventY(ii), '.', 'markersize', HMMprctile_PBEshuffle(ii)*30)
% % end
% % 


%% Functions

function [eventsBinnedfiring, secondaryPBEs, idx_acceptedEvents] = finalBinningResult(pbEvents2, mySpikes, fileinfo)

    %% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events


%%% Calculating the number of firing pyramidal units within each event

qclus = [1 3]; % only pyramidal

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
idx_acceptedEvents = find(noFiringUnits >= max(4, floor(0.1 * length(activeUnits))) & (eventLen >= 4)); % maximum length 500 ms (???)
secondaryPBEs = [eventBeg(idx_acceptedEvents) eventEnd(idx_acceptedEvents) pbEvents2(idx_acceptedEvents, 3:4)]; 


eventsBinnedfiring = eventsBinnedfiring(idx_acceptedEvents, :); % considering only the qualified event for the analyses


end


function savePBEs(primaryPBEs, secondaryPBEs, acceptedEvtsIdx, analysisPeriod, fileinfo, t0) %#ok<INUSL>


currDir = pwd;
Folderbase = [currDir '/' fileinfo.name '/PopulationBurstEvents/' analysisPeriod];
mkdir(Folderbase)

save([Folderbase '/' fileinfo.name  '-' 'PBEs.mat'], 'primaryPBEs', 'secondaryPBEs', 'acceptedEvtsIdx') 

Filename = [Folderbase '/' fileinfo.name '.prm.evt'];
MakeEvtFile(primaryPBEs(:, 1:3)-t0, Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)

Filename = [Folderbase '/' fileinfo.name '.scn.evt'];
MakeEvtFile(secondaryPBEs(:, 1:3)-t0, Filename, {'beg', 'end', 'peak'}, fileinfo.Fs, 1)


end


function [poissonEventsBinnedfiring, timeSwapEventsBinnedfiring, temporalEventsBinnedfiring, pooledtimeswapEventsBinnedfiring] = genSurrogates(eventsBinnedfiring)


%% Generating shuffle surrogate PBE datasets

% generate poisson simulated dataset

binDur = 0.02;
noEvents = size(eventsBinnedfiring, 1);

poissonEventsBinnedfiring = poissonSpikeTrain(eventsBinnedfiring, binDur); 
% cmprPoisson2Raw(poissonEventsBinnedfiring, eventsBinnedfiring, longORshort, fileinfo, binDur); %%% compare the poisson simulated data with the actual data distribution

% generate time swap dataset (coherent shuffle within each event, keeping the cofirings the same)

timeSwapEventsBinnedfiring = timeswap(eventsBinnedfiring, binDur);


% time bins temporal shuffle (independently for each unit, incoherent shuffle)

temporalEventsBinnedfiring = cell(noEvents, 2);
for evt = 1 : noEvents
    
    currEvt = eventsBinnedfiring(evt, :);
    
    randShift = randi(size(currEvt{2}, 2)-1, 1, size(currEvt{2}, 1)); % generate shift amount for each row
    
    
    % Because we are doing row cycle shuffle, we need to transpose the data
    % matrix and after doing column cycle shuffles transpose it again
    
    temporalEventsBinnedfiring{evt, 2} = column_cycle_shuffle(currEvt{2}', randShift)';
    temporalEventsBinnedfiring{evt, 1} = column_cycle_shuffle(currEvt{1}', randShift*20)'; % note the shift here
    
end

pooledtimeswapEventsBinnedfiring = cell(noEvents, 2);
pooledEvents = cell2mat(eventsBinnedfiring(:, 2)');
noallBins = size(pooledEvents, 2);

for evt = 1:noEvents
    
    currEvt = eventsBinnedfiring{evt, 2};
    
    for bin = 1 : size(currEvt, 2)
        
        pooledtimeswapEventsBinnedfiring{evt, 2}(:, bin) = pooledEvents(:, randi(noallBins, 1));
        
    end


end

end



function congruencePlot(HMMprctile, HMMprctile_p, HMMprctile_ts, HMMprctile_tc, HMMprctile_pts,  okLegend, curr_title)

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

h_pts = hist(HMMprctile_pts, bins);
cumhist_pts = cumsum(h_pts);
crosspnt_pts = cumhist_pts(threshBin);

% h_pbe = hist(HMMprctile_PBEshuffle, bins);
% cumhist_pbe = cumsum(h_pbe);
% crosspnt_pbe = cumhist_pbe(threshBin);


p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9], 'linewidth', 4);

p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

p2 = plot(bins, cumhist_ts, 'color', [100 180 255]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_ts crosspnt_ts], 'color', [100 180 255]/255,'LineStyle', '--', 'linewidth', 2)

p3 = plot(bins, cumhist_p, 'color', [255 180 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_p crosspnt_p], 'color', [255 180 100]/255,'LineStyle', '--', 'linewidth', 2)

p4 = plot(bins, cumhist_tc, 'color', [255 100 100]/255, 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_tc crosspnt_tc], 'color', [255 100 100]/255,'LineStyle', '--', 'linewidth', 2)

p5 = plot(bins, cumhist_pts, 'color', 'k', 'linewidth', 4);
line([0 bins(threshBin)], [crosspnt_pts crosspnt_pts], 'color', 'k','LineStyle', '--', 'linewidth', 2)

% p6 = plot(bins, cumhist_pbe, 'color', 'k', 'LineStyle', '-.', 'linewidth', 4);
% line([0 bins(threshBin)], [crosspnt_pbe crosspnt_pbe], 'color', 'k','LineStyle', '--', 'linewidth', 2)

noEvents = length(HMMprctile);
line([bins(threshBin) bins(threshBin)], [0 noEvents], 'color', 'k','LineStyle', '--', 'linewidth', 2)


if okLegend
    legend([p1 p2 p3 p4 p5 p0],{'actual','time swap','poisson','temporal shuffle', 'pooled ts', 'pbe shuffle', 'chance line'}, 'Location', 'northwest')
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


function HMMprctile = HMMcongruence_crossvalid_PBEshuffle(data, testEvts, transmat, lambda, noShuffle, dataType)



noEvents = size(data, 1);
numofFolds = length(testEvts);
numofStates = size(transmat, 1);



dataLL = zeros(1, noEvents); %% log likelihood of the raw data
nullLL = zeros(noEvents, noShuffle); %% log likelihood of shuffle data

HMMprctile = zeros(noEvents, 1);

pooledEvents = cell2mat(data(:, 2)');
noallBins = size(pooledEvents, 2);


for fold = 1 : numofFolds
    
    curr_transmat = transmat(:, :, fold);
    curr_lambda = lambda(:, :, fold);


    
    for evt = testEvts{fold}
        
        currEvent = data{evt, 2};
        noBins = size(currEvent, 2);

        B = poisson_prob(currEvent, curr_lambda,1);
        prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, ~,  dataLL(evt), ~] = fwdback(prior, curr_transmat, B);
        

        %%%% null distribution of observation likelihhod
        
        for sn = 1 : noShuffle
            
            shuffledEvent = pooledEvents(:, randi(noallBins, [1 , noBins]));
            B = poisson_prob(shuffledEvent, curr_lambda,1);

            [~, ~, ~, nullLL(evt,sn), ~] = fwdback(prior, curr_transmat, B);
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
            HMMprctile(evt, method) = length(find(nullLL(evt,:,method) < dataLL(evt)))/noShuffle;

        end
        
        
        if mod(evt, 50) == 0
            fprintf(1, 'cross-validation_%s_event %d, HMMPercentile = %f\n', dataType, evt, HMMprctile(evt, 1));
        end
        
      
    end

end
% 
% save([FileBase '/cvScores_' period '_' dataType '.mat'], 'HMMPercentile', 'dataLL', 'nullLL')
% 
% % save([FileBase '/' dataType '.mat'], 'HMMScoreZ', 'dataLL', 'nullLL', 'partialLikeli', 'bins_nullDistPercentiles', 'binsCongSignificance', 'maxCongruenceLength', 'totalCongruenceLength')
% 

end


function HMMPercentile = HMMcongruence_crossvalid_pts(data, testEvts, transmat, lambda, noShuffle, dataType, fileinfo, period)


currDir = pwd;
% FileBase = [currDir '/' fileinfo.name '/data part' num2str(whichPart) '/HMM/sequenceDetection/measure congruence/'];
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

        shuffle_transmat(s,:) = transmat(s, shuffleInd); % shuffling within the columns
    end
    
    % if doing column-wise shuffle
    
%     shuffleTransmats(:, :, sn) = shuffle_transmat./repmat(sum(shuffle_transmat, 2), [1, numofStates]);
       shuffleTransmats(:, :, sn) = shuffle_transmat;
    
end

end


function [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData(mySpikes, activePeriods, behPeriods, binDur, posBinSize, fileinfo)


%% binning the run data


qclus = [1 3]; % pyramidal neurons

runBinnedfiring2 = timeBinning(activePeriods , mySpikes, qclus, binDur, fileinfo.Fs); 



% Calculate the track positions corresponding to each time bin

% periods = activePeriods * 1/fileinfo.lfpSampleRate;

periods = activePeriods;
runBinPos = cell(size(periods, 1), 1);

binDur = binDur * 1e6;

for ii = 1 : size(periods, 1)
    
%     runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/fileinfo.Fs, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3), fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin

    runBinPos{ii}(end, :) = [];
end
    

% runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields


% posBinSize = 0.01; % in cm
runidx = find(fileinfo.xyt(:, 3) > behPeriods(4,1) & fileinfo.xyt(:, 3) < behPeriods(8,2));

xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);


% noXPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
% noYPosBins = floor((max(ypos) - min(ypos))/posBinSize);

xposBins = min(xpos): posBinSize: max(xpos); % center of the position bins
xposBins(end) = max(xpos);

xposcenters = xposBins + posBinSize/2;
xposcenters(end) = [];


yposBins = min(ypos): posBinSize:max(ypos); 
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



function [runBinnedfiring, posbinIdx, xposcenters, yposcenters] = binRunData2(mySpikes, activePeriods, behPeriods, binDur, posBinSize, exclude, fileinfo)


%% binning the run data


qclus = [1 3]; % pyramidal neurons

runBinnedfiring2 = timeBinning(activePeriods , mySpikes, qclus, binDur, fileinfo.Fs); 


% exclude the PBEs





% Calculate the track positions corresponding to each time bin

% periods = activePeriods * 1/fileinfo.lfpSampleRate;

periods = activePeriods;
runBinPos = cell(size(periods, 1), 1);

binDur = binDur * 1e6;

for ii = 1 : size(periods, 1)
    
%     runBinPos{ii} = interp1(fileinfo.xyt(:, 3)*1/fileinfo.Fs, fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin
    runBinPos{ii} = interp1(fileinfo.xyt(:, 3), fileinfo.xyt(:,[1 2]), (periods(ii, 1):binDur:periods(ii, 2))+binDur/2); %% at the center of each time bin

    runBinPos{ii}(end, :) = [];
end
    

% runBinPos = cell2mat(runBinPos');


% Go from the continuous position to position bins

% Defining the position bins; the same we had in calculating the place
% fields

% posBinSize = 0.01; % in cm
runidx = find(fileinfo.xyt(:, 3) > behPeriods(4,1) & fileinfo.xyt(:, 3) < behPeriods(8,2));

xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);


% noXPosBins = floor((max(xpos) - min(xpos))/posBinSize); % posBinSize: the length of each position bin in cm
% noYPosBins = floor((max(ypos) - min(ypos))/posBinSize);

xposBins = min(xpos): posBinSize: max(xpos); % center of the position bins
xposBins(end) = max(xpos);

xposcenters = xposBins + posBinSize/2;
xposcenters(end) = [];


yposBins = min(ypos): posBinSize:max(ypos); 
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


function [lsPFs, max_prob] = plot_lsPF(runBinnedfiring, transmat, lambda, posbinIdx, xposcenters, yposcenters, numofStates, fileinfo, periods, pbePeriod)


% numofStates = 40;
lsPFs = HMM_StatesPlaceField_2D(runBinnedfiring{2}, transmat, lambda, posbinIdx, [], numofStates, xposcenters, yposcenters);



% gamma_avg = permute(gamma_avg, [2 1 3]); 
% clim = [min(lsPFs(:)) max(lsPFs(:))];


% plot the lsPFs

figure;

numStatesperRow = 6;
boxThreshold = 120; 

runidx = find(fileinfo.xyt(:, 3) > periods(4,1) & fileinfo.xyt(:, 3) < periods(8,2) & fileinfo.xyt(:, 1) > boxThreshold);
xpos = fileinfo.xyt(runidx,1);
ypos = fileinfo.xyt(runidx,2);

max_prob = zeros(numofStates, 1);

for ii = 1: numofStates
    
    subplot(ceil(numofStates/numStatesperRow), numStatesperRow, ii)
    
    plot(fileinfo.xyt(runidx, 1), fileinfo.xyt(runidx, 2), '.', 'markersize', 1, 'MarkerEdgeColor', [0.7 0.7 0.7])
    
    hold on
    h1 = imagesc(xposcenters, yposcenters, lsPFs(:,:, ii));
    set(h1, 'AlphaData', 0.7)
    max_prob(ii) = max(max(lsPFs(:,:, ii)));
    text(min(xposcenters), max(yposcenters)+10, sprintf('%.1f', max_prob(ii)*100))

    
    colormap('hot')
    set(gca,'YDir','normal')
    
    cornerstate = (ceil(numofStates/numStatesperRow)-1) * numStatesperRow + 1;
    
    if ii == cornerstate
        
%         set(gca, 'XTick', [min(xpos) max(xpos)], 'YTick', [min(ypos) max(ypos)], 'box', 'off')
        set(gca, 'box', 'off')
        xticks([min(xpos) max(xpos)])
        xticklabels({sprintf('%.1f',  min(xpos))  sprintf('%.1f',  max(xpos))})
        
        yticks([min(ypos) max(ypos)])
        yticklabels({sprintf('%.1f',  min(ypos))  sprintf('%.1f',  max(ypos))})
        
        xlabel('X (cm)', 'fontsize', 14)
        ylabel('Y (cm)', 'fontsize', 14)
        
    else
        
        set(gca, 'XTick', [], 'YTick', [], 'box', 'off')
    end

    xlim([min(xpos) max(xpos)])
    ylim([min(ypos) max(ypos)])

end

% suptitle([fileinfo.name '-' pbePeriod])
set(gcf, 'position', [341.8105   49.2632  730.6105  930.6947])

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
    
    
    train_shuffle_posbinIdx = train_posbinIdx(randperm(size(trainSet, 1), size(trainSet, 1)), :);
    shuffle_train_lsPFs = HMM_StatesPlaceField_2D(train_timeBins, transmat, lambda, train_shuffle_posbinIdx, [], numofStates, xposcenters, yposcenters);
    
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
    
%     decodedPositions_shuffledGamma(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, train_lsPFs, xposcenters,yposcenters, 'shuffledGamma'); 
    
    decodedPositions_shuffledGamma(testSet, :) = lsPFPositionDecode(obsProb, transmat, prior, shuffle_train_lsPFs, xposcenters,yposcenters, 'actual'); 

    

        
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
    
    % not sure abou the following line 
%     probabilityOverPosition = probabilityOverPosition ./ repmat(sum(sum(probabilityOverPosition, 1),2), [size(probabilityOverPosition, 1) size(probabilityOverPosition, 2)]);
    
    

    decodedPositions(jj, 1) = yposcenters(ceil(sum((1:length(yposcenters))' .* sum(probabilityOverPosition, 2))));
    decodedPositions(jj, 2) = xposcenters(ceil(sum((1:length(xposcenters)) .* sum(probabilityOverPosition, 1))));
end


end


function plotErrorCdf(data, shuffledOrder, shuffledGamma, ylabelneeded)

hold on

bins = linspace(min([data; shuffledGamma]), max([data; shuffledGamma]), 100);

% actual data
% [h1, bins] = hist(data, 100);
h1 = hist(data, bins);

h1c = cumsum(h1)/sum(h1);

dataMed = median(data);
 
curve1 = plot(bins, h1c, 'linewidth', 4, 'color', [150, 150, 255]/255);
line([dataMed dataMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
line([0 dataMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 150, 255]/255)
text(dataMed, 0.03, sprintf('%.2f', dataMed), 'color', [150, 150, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')

% shuffled order
[h2, bins] = hist(shuffledOrder, 100);
h2c = cumsum(h2)/sum(h2);

shuffleMed = median(shuffledOrder);

curve2 = plot(bins, h2c, 'linewidth', 4, 'color', [150, 200, 255]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [150, 200, 255]/255)
text(shuffleMed, 0.07, sprintf('%.2f', shuffleMed), 'color', [150, 200, 255]/255, 'fontsize', 12, 'FontWeight', 'Bold')


% shuffled Gamma
[h3, bins] = hist(shuffledGamma, 100);
h3c = cumsum(h3)/sum(h3);

shuffleMed = median(shuffledGamma);

curve3 = plot(bins, h3c, 'linewidth', 4, 'color', [160, 160, 160]/255);
line([shuffleMed shuffleMed],[0 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
line([0 shuffleMed],[0.5 0.5], 'linewidth', 2, 'LineStyle', '--', 'color', [160, 160, 160]/255)
text(shuffleMed, 0.1, sprintf('%.2f', shuffleMed), 'color', [160, 160, 160]/255, 'fontsize', 12, 'FontWeight', 'Bold')




xlim([0 max([data; shuffledOrder; shuffledGamma])])
set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01])


xlabel('decoding error(cm)', 'fontsize', 20)

if ylabelneeded
    ylabel('cumulative ratio', 'fontsize', 20)
    
%     legend('actual', 'shuffle state probability', 'shuffle position index', 'Location', 'northwest'); 
    legend([curve1, curve2, curve3], 'actual', 'shuffled run bins', 'shuffled lsPFs', 'Location', 'northwest'); 

    legend boxoff 
end

axis square
end

