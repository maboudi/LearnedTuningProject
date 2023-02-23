% clear; clc
% 
% cd('/home/kouroshmaboudi/Documents/LaurelData/ObjectLocationMemory/Charlie')
% 
% % Dataset directory
% currDir = pwd;
% 
% %% session info
% 
% animal = 'charlie';
% sessionName = 'PlacePreference-Charlie-20150624-01-nlx';
% 
% fileinfo = struct('name', sessionName, 'animal', animal, 'xyt', [], 'tbegin', [], 'tend', ...
%                     [], 'Fs', [], 'nCh', [], 'lfpSampleRate', [], 'CA1thetach', [], 'CA', [1 1 1 1]); 
% 
% 
% % (1)  Sleep 1: ~20 minute baseline sleep in homecage.
% % (2)  Exposure 1: exposure to empty box (~6 min).
% % (3)  Homecage 1: return to home cage (~3 min).
% % (4)  Exposure 2: exposure to box with objects (~6 min).
% % (5)  Homecage 2: return to homecage (~3 min). 
% % (6)  Exposure 3: exposure to box with objects (~6 min).
% % (7)  Homecage 3: return to homecage (~3 min). 
% % (8)  Exposure 4: exposure to box with objects (~6 min). 
% % (9)  Sleep 2: extended sleep (~4 hours). 
% % (10) Memory recall test: exposure to box with displaced object (~6 min).
% % (11) rest of the reocring for a short period (~20 min)
% 
% 
% FileBase = [currDir '/' fileinfo.name '/' fileinfo.name];
% 
% % Loading the mat fiels
% % varList = {'basics'; 'speed'; 'spikes'; 'hiAmpSpindles'};
% varList = {'-basics'; '-speed'; '-spikes'; '.SleepState.states'};
% 
% for var = 1:length(varList)
%     load([FileBase  varList{var} '.mat'])
% end
% 
% 
% % fileinfo.Fs = basics.SampleRate;
% 
% fileinfo.Fs = 1e6; % we don't use the original sample rate, because the time stamps were converted to micro sec. 
% fileinfo.lfpSampleRate = basics.lfpSampleRate; 
% fileinfo.nCh = basics.nChannels;
% 
% fileinfo.xyt = [speed.x; speed.y; speed.t]'; 
% 
% 
% % periods
% periods = basics.period.time;
% 
% 
% %% processing spikes
% 
% nUnits = length(spikes);
% 
% for unit = 1:length(spikes)
%     
%     spikes(unit).rate = length(spikes(unit).time)/((periods(end, 2) - periods(1,1))/1e6);
% 
% end
% 
% 
% 
% spike = struct('t', [], 'unit', [], 'shank', [], 'qclu', [], 'x', [], 'y', [], 'speed', [], 'rate', []);
% 
% for unit = 1: nUnits
%     
%     if spikes(unit).quality == 0
%         continue
%     end
%     
%     unitSpikes = spikes(unit).time;
%     
%     spike.t = [spike.t; unitSpikes];
%     
%     spike.unit = [spike.unit; repelem(unit, length(unitSpikes))'];
%     spike.shank = [spike.shank; repelem(spikes(unit).id(1), length(unitSpikes))'];
%     spike.qclu = [spike.qclu; repelem(spikes(unit).quality, length(unitSpikes))'];
%     spike.rate = [spike.rate; repelem(spikes(unit).rate, length(unitSpikes))'];
%     
%     spike.x = [spike.x; spikes(unit).x];
%     spike.y = [spike.y; spikes(unit).y];
%     spike.speed = [spike.speed; spikes(unit).speed];
%     
%     
% end
% 
% 
% % sorting the spike time stamps
% [spike.t, spikeSortIdx] = sort(spike.t, 'ascend');
% 
% spike.unit = spike.unit(spikeSortIdx);
% spike.shank = spike.shank(spikeSortIdx);
% spike.qclu = spike.qclu(spikeSortIdx);
% spike.rate = spike.rate(spikeSortIdx);
% 
% spike.x = spike.x(spikeSortIdx);
% spike.y = spike.y(spikeSortIdx);
% spike.speed = spike.speed(spikeSortIdx);
% 
% 
% spike.qclu(spike.rate > 2 & ismember(spike.qclu, [1 3])) = 2;
% 
% %% PBEs
% 
%     %% detemining population burst periods
% 
%     
% % exclude all the periods with speed above 5 cm/s
% 
% velocity = speed.v;
% 
% % smoothing the speed; the sampling rate is 100 Hz (after interpolation).
% 
% sigma = 100; % 1 sec
% halfwidth = 3*sigma;
% win = gausswindow(sigma, halfwidth); % smoothing kernel
% 
% velocity = conv(velocity, win, 'same');  
% 
% % exclude theta, or in case of Charlie (day 1) periods with high amplitude
% % spindles
% 
% exclude = []; 
% 
% % exclude = [exclude HiAmpSpindles/fileinfo.lfpSampleRate*fileinfo.Fs]; % high volatge spindles
%     
% 
% 
% 
% time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
% threshZ = 2; % sdf with 3 std deviation above the meanposbinIdx(:,1)
% 
% 
% 
% % memory test period
% % 
% % pers = [4 5 6 7 8 9 10 11];
% % combinedPBEs = [];
% % 
% % for ii = 1:length(pers)
% %     
% % fileinfo.tbegin = periods(pers(ii),1); 
% % fileinfo.tend = periods(pers(ii),2);
% % % fileinfo.tend = periods(9,1)+2*60*60*fileinfo.Fs; %% first two hours of sleep
% % 
% % 
% % [primaryPBEs, sdat] = PBPeriods(spike, fileinfo, [], time_resolution, threshZ, exclude);
% % eventVelocity = interp1(speed.t, velocity, primaryPBEs(:,3));
% % primaryPBEs2 = primaryPBEs(eventVelocity < 5, :); % the exclusion criterion here is based on speed, ideally should be based on the theta
% %                                                   % some high attentive
% %                                                   % states could be
% %                                                   % accompanied with high
% %                                                   % amplitude theta
% %                                                   % oscillations
% %                                                   
% % 
% % combinedPBEs = [combinedPBEs; primaryPBEs2];  
% % 
% % end
% %                                                   
% %                                              
% %                                                   
% %                                                   
% % noEvents = size(combinedPBEs, 1)
% 
% 
% % % exclude the interneurons and MUAs
% % 
% % INspikeidx = find(ismember(spike.qclue, [4 8]));
% % 
% % spike.t(INspikeidx) = [];
% % spike.unit(INspikeidx) = [];
% % spike.shank(INspikeidx) = [];
% % spike.qclu(INspikeidx) = [];
% % spike.x(INspikeidx) = [];
% % spike.y(INspikeidx) = [];
% % spike.speed(INspikeidx) = [];
% 
% 
% % Binning and assessing the duration and participation ratio for each
% % event; excluding the events with duration less than 4 time bins(80 ms)
% % and participation of less than 10% of active units
% 
% 
% % remove all of the non-pyraidal units
% remIdx = find(~ismember(spike.qclu, [1 3]));
% 
% spike.t(remIdx) = [];
% spike.qclu(remIdx) = [];
% spike.unit(remIdx) = [];
% spike.shank(remIdx) = [];
% spike.rate(remIdx) = [];
% 
% spike.x(remIdx) = [];
% spike.y(remIdx) = [];
% spike.speed(remIdx) = [];
% 
% % 
% % 
% % [eventsBinnedfiring, secondaryPBEs, acceptedprmEvts] = finalBinningResult(combinedPBEs, spike, fileinfo);
% % 
% % % eventX = interp1(speed.t, speed.x, secondaryPBEs(:,3));
% % % eventY = interp1(speed.t, speed.y, secondaryPBEs(:,3));
% % 
% % savePBEs(combinedPBEs, secondaryPBEs, acceptedprmEvts, 'test', fileinfo, periods(1,1))
% % 
% % [poissonEventsBinnedfiring, timeSwapEventsBinnedfiring, temporalEventsBinnedfiring, pooledtimeswapEventsBinnedfiring] = genSurrogates(eventsBinnedfiring);
% % 
% 
% 
% % Ripple detection
% 
% 
% % Here we intended to detect only ripples that have on average relatively large
% % amplitudes over the all four shanks resided in CA1. Alternatively, we can
% % include also ripple events (with usually smaller amplitude) picked by
% % just one or a subset of shanks. 
% 
% best_channels = BigRippleChannels(fileinfo,1);
% fileinfo.bestch = best_channels;
% 
% rippPowThresh = 3; 
% rippleEvents = RippleDetect(fileinfo, rippPowThresh); % in lfpSampleRate
% 
% Filename = [FileBase '.spl.evt'];
% MakeEvtFile(rippleEvents(:, 1:3), Filename, {'beg', 'peak', 'end'}, fileinfo.lfpSampleRate, 1)
% 
% eventVelocity = interp1(speed.t, velocity, (rippleEvents(:,2)/fileinfo.lfpSampleRate *1e6) + periods(1,1));
% rippleEvents2 = rippleEvents(eventVelocity < 5, :); % refining the ripples events; excluding the ones that happen during mobile state


% % calculating the overlap between the detected population burst events and
% % the ripples (percentage)
% 
% 
% rplOverlap = zeros(noEvents, 1);
% 
% noRipples = size(rippleEvents2, 1);
% 
% Ripples = rippleEvents2;
% Ripples(:, 1:3) = (rippleEvents2(:, 1:3)/fileinfo.lfpSampleRate) + periods(1,1)/fileinfo.Fs; % in sec
% PBEvents = combinedPBEs/fileinfo.Fs; % in sec
% 
% for evt = 1 : noEvents
%     currEvt = PBEvents(evt,:);
%     for rip = 1:noRipples
%         
%         if (currEvt(1) >= Ripples(rip, 1) && currEvt(1) <= Ripples(rip, 3)) || ...
%                 (currEvt(2) >= Ripples(rip, 1) && currEvt(2) <= Ripples(rip, 3))...
%                 || (currEvt(1) <= Ripples(rip, 1) && currEvt(2) >= Ripples(rip, 3))
%             
%             rplOverlap(evt) = 1;
%             break
%         end
%     end
% end
% 
% PopBurst_Ripple_overlap = length(find(rplOverlap))/length(PBEvents)
% 


% if we are going to use the ripple events instead of the population burst
% events for further process

periodRipples2 = rippleEvents2;
periodRipples2(:, 1:3) = periodRipples2(:, 1:3)/fileinfo.lfpSampleRate *1e6 + periods(1,1);

periodRipples = periodRipples2;
periodRipples(:,2) = periodRipples2(:,3);
periodRipples(:,3) = periodRipples(:,2);

periodRipples = periodRipples(find(periodRipples(:, 2) >= periods(9, 1) & periodRipples(:,2) <= periods(9, 1)+2.5*60*60*1e6), :);

% 
% 
% % Passing only the ripples that happens during wake (immobile wake as a
% % proxy for quiet wake) periods
% 
% wakePeriods = SleepState.ints.WAKEstate;
% wakePeriods = [(wakePeriods(:, 1)-1)*1e6+1 wakePeriods(:,2)*1e6]+periods(1,1);
% 
% wakePeriods(4:end, :) = [];
% 
% 
% qwIdx = zeros(size(periodRipples, 1), 1);
% for rip = 1: size(periodRipples, 1)
%     
%     for ii = 1: size(wakePeriods, 1)
%         
%         if periodRipples(rip, 1) >= wakePeriods(ii, 1) & periodRipples(rip, 2) < wakePeriods(ii, 2)
%             qwIdx(rip) = 1;
%             break
%         end
%     end
% end
% 
% 
% qwRippls = periodRipples(find(qwIdx), :);


[eventsBinnedfiring, secondaryrpl, acceptedprmEvts] = finalBinningResult(periodRipples, spike, fileinfo);
% [~, secondaryrpl, acceptedprmEvts] = finalBinningResult(rippleEvents2, spike, fileinfo);


Filename = [FileBase '.rpl.evt'];
MakeEvtFile(secondaryrpl(:, 1:3)-periods(1,1), Filename, {'beg', 'peak', 'end'}, fileinfo.Fs, 1)

%% HMM sequence detection

% In this step we test the extent to which an HMM trained on the
% candidate events is able to detect sequence events among non-sequence
% events
% 
% we do this by k-fold cross-validation: using k-1 folds for training and
% and the remaining fold for sequence detection %% this is temporary assigning. Later will change it when analyzing run/post-sleep data

% first of all randomize the events



numofFolds = 5; %% 5-fold cross-validation
numofStates = 30;
noShuffle = 1000;


acc = [];
acc_pts = [];

for ii = 1:10
    
    ii
    
eventsBinnedfiring = eventsBinnedfiring(randperm(size(eventsBinnedfiring, 1)), :);

% test

pbePeriod = 'combined_exposures+interleavedrest';

[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiring, numofStates, numofFolds, fileinfo, pbePeriod); 

HMMprctile2 = HMMcongruence_crossvalid(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual', fileinfo, pbePeriod);
% HMMprctile_p = HMMcongruence_crossvalid(poissonEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'poisson', fileinfo, pbePeriod);
% HMMprctile_ts = HMMcongruence_crossvalid(timeSwapEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap', fileinfo, pbePeriod);
% HMMprctile_t = HMMcongruence_crossvalid(temporalEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'temporal', fileinfo, pbePeriod);


% Pooled time swap with replacement and resampling the bins from all the
% events not only the test set; we need to replace this surrogate with the
% one with sampling with replacement and sampling from only test set
% HMMprctile_pts2 = HMMcongruence_crossvalid(pooledtimeswapEventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'timeswap_pooled', fileinfo, pbePeriod);

HMMprctile_pts2 = HMMcongruence_crossvalid_pts(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'pooledTimeSwap', fileinfo, pbePeriod);
% 
% % HMMprctile_PBEshuffle = HMMcongruence_crossvalid_PBEshuffle(eventsBinnedfiring, testEvts, CVtransmats, CVlambdas, noShuffle, 'actual');
% 

acc = [acc; HMMprctile2];
acc_pts = [acc_pts; HMMprctile_pts2];
end

figure;
set(gcf, 'position', [801 82 1648 1239])

% cross-validations (within-periods)

% congruencePlot(HMMprctile, HMMprctile_p, HMMprctile_ts, HMMprctile_t, HMMprctile_pts, 1, 'Extended sleep (4 hr) Cross-Validation')


% flip the curve like the one in the paper (something like 1 minus the curve in the top figure generated by the function congruencePlot)
% congruencePlot2(HMMprctile2, HMMprctile_pts2, 1, 'combined exposures and rest Cross-Validation')
congruencePlot2(acc, acc_pts, 1, 'combined exposures and rest Cross-Validation')


figure;
idx = find(speed.t > fileinfo.tbegin & speed.t < fileinfo.tend);
plot(speed.x(idx), speed.y(idx), 'linewidth', 1)

hold on

for ii = 1:length(secondaryPBEs)
    plot(eventX(ii), eventY(ii), '.', 'markersize', HMMprctile_PBEshuffle(ii)*30)
end



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
idx_acceptedEvents = find(noFiringUnits >= floor(0.1 * length(activeUnits)) & (eventLen >= 4)); % maximum length 500 ms (???)
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

sigma = 3;
halfwidth = 2*sigma;

win = gausswindow(sigma, halfwidth);

h = conv(h, win, 'same');



cumhist = 100-cumsum(h); %%  for the paper
% threshBin = find(bins > 99, 1, 'first');
% crosspnt = cumhist(threshBin);


h_pts = hist(HMMprctile_pts*100, bins);
h_pts = h_pts/sum(h_pts)*100;


h_pts = conv(h_pts, win, 'same');

cumhist_pts = 100-cumsum(h_pts); %% for the paper
% crosspnt_pts = cumhist_pts(threshBin);



% p0 = line([0 1], [0 length(HMMprctile)], 'color', [0.9 0.9 0.9],
% 'linewidth', 4); % the chance line-what would be the chance distribution
% is not established yet

% p1 = plot(bins, cumhist, 'color', [0 200 100]/255, 'linewidth', 4);
p1 = area(bins, h, 'FaceColor', [0 200 100]/255, 'FaceAlpha', 0.3, 'EdgeColor', [0 200 100]/255, 'linewidth', 3);


% line([0 bins(threshBin)], [crosspnt crosspnt], 'color', [0 200 100]/255,'LineStyle', '--', 'linewidth', 2)

% p2 = plot(bins, cumhist_pts, 'color', 'b', 'linewidth', 4);
p2 = area(bins, h_pts, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', [0.7 0.7 0.7], 'linewidth', 3);


% line([0 bins(threshBin)], [crosspnt_pts crosspnt_pts], 'color', 'b','LineStyle', '--', 'linewidth', 2)


% noEvents = length(HMMprctile);
% line([bins(threshBin) bins(threshBin)], [0 1], 'color', 'k','LineStyle', '--', 'linewidth', 2)


if okLegend
    legend([p1 p2],{'actual', 'pooled time swap'}, 'Location', 'northwest')
    legend boxoff 
end

set(gca, 'fontsize', 14, 'linewidth', 3, 'TickDir', 'out', 'TickLength',[0.02, 0.01], 'Layer', 'top')
xlabel('percentile threshold', 'fontsize', 16)


if okLegend
    ylabel('% sig. PBEs', 'fontsize', 16)
end

title(curr_title, 'fontsize', 14)

% ylim([0 100])
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