clear; clc

cd('/home/kouroshmaboudi/Documents/datasets/LaurelData/ObjectLocationMemory/Charlie')

% Dataset directory
currDir = pwd;

%% session info

animal = 'charlie';
sessionName = 'PlacePreference-Charlie-20150623-01-nlx';

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
% (11) rest of the recording for a short period (~20 min)


FileBase = [currDir '/' fileinfo.name '/' fileinfo.name];

% Loading the mat files
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


spike.qclu(spike.rate > 2 & ismember(spike.qclu, [1 3])) = 2; % remove the high firing pyramidal units




%% PBEs

%% detemining population burst periods

    
% exclude all the periods with speed above 5 cm/s

velocity = speed.v;

% smoothing the speed; the sampling rate is 100 Hz (after interpolation).

sigma     = 100; % 1 sec
halfwidth = 3 * sigma;
win       = gausswindow(sigma, halfwidth); % smoothing kernel

velocity = conv(velocity, win, 'same');  

% exclude theta, or in case of Charlie (day 1) periods with high amplitude
% spindles

exclude = []; 

% exclude = [exclude HiAmpSpindles/fileinfo.lfpSampleRate*fileinfo.Fs]; % high volatge spindles


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the meanposbinIdx(:,1)



% The extended sleep period following the last exposure to the task box

pers = 9;
combinedPBEs = [];

for ii = 1:length(pers)
    
fileinfo.tbegin = periods(pers(ii),1); 
fileinfo.tend   = periods(pers(ii),2);

velocityFilter = 0;
[primaryPBEs, sdat] = PBPeriods(spike, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);
eventVelocity       = interp1(speed.t, velocity, primaryPBEs(:,3));
primaryPBEs2        = primaryPBEs(eventVelocity < 5, :); % the exclusion criterion here is based on speed, ideally should be based on the theta
                                                         % some high attentive states could be accompanied with high amplitude theta oscillations 
                                            
combinedPBEs = [combinedPBEs; primaryPBEs2];  

end
                                            
                                                  
noPrimEvents = size(combinedPBEs, 1);



% Binning and assessing the duration and participation ratio for each
% event; excluding the events with duration less than 4 time bins(80 ms)
% and participation of less than 10% of active units


% remove all of the non-pyraidal units

keepIdx = find(ismember(spike.qclu, [1 3]));
spike   = structfun(@(x)(x(keepIdx)), spike, 'UniformOutput', false);

allActiveUnits = unique(spike.unit);
[~, sortIdx] = sort(allActiveUnits, 'ascend');


temp = spike.unit;

for ii = 1: length(sortIdx)
    temp(find(spike.unit == allActiveUnits(ii))) = sortIdx(ii);
end

spike.unit = temp;

included_qclus = [1 3];
[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(combinedPBEs, spike, included_qclus, fileinfo);


noSecEvents = size(secondaryPBEs, 1);



%% Ripple detection


% % Here we intended to detect only ripples that have on average relatively large
% % amplitudes over the all four shanks resided in CA1. Alternatively, we can
% % include also ripple events (with usually smaller amplitude) picked by
% % just one or a subset of shanks. 
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


%% Detected PBEs overlap with the ripples


% rplOverlap = zeros(noSecEvents, 1);
% 
% noRipples = size(rippleEvents2, 1);
% 
% Ripples = rippleEvents2;
% Ripples(:, 1:3) = (rippleEvents2(:, 1:3)/fileinfo.lfpSampleRate) + periods(1,1)/fileinfo.Fs; % in sec
% PBEvents = secondaryPBEs/fileinfo.Fs; % in sec
% 
% for evt = 1 : noSecEvents
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
% PopBurst_Ripple_overlap = length(find(rplOverlap))/length(PBEvents);



%% HMM sequence detection

% In this step we test the extent to which an HMM trained on the
% candidate events is able to detect sequence events among non-sequence
% events
% 
% we do this by k-fold cross-validation: using k-1 folds for training and
% and the remaining fold for sequence detection %% this is temporary assigning. Later will change it when analyzing run/post-sleep data

% first of all randomize the events


nFolds = 3; %% 5-fold cross-validation
nStates = 20;
nShuffles = 100;
sm = 1;

[HMMSeqScores.data(:, sm), HMMSeqScores.pts(:, sm), statesProbDists] = HMMCrossValidAnalysis(binnedPBEs,  nStates, nShuffles, nFolds, currDir, 'PRE', sm);


%% plot example PBEs

nPBEs2Plot = 100;

% [~, sortInd] = sort(HMMSeqScores.data, 'descend');
% selectPBEs = sortInd(1:nPBEs2Plot);

selectPBEs = randperm(size(binnedPBEs, 1), nPBEs2Plot);
% selectPBEs = examplePBEs;
figure;
set(gcf, 'position', [1 35 1440 2406])

firstRow  = [1:10 21:30 41:50 61:70 81:90 101:110 121:130 141:150 161:170 181:190];
secondRow = firstRow + 10;

gaussianWin = gausswindow(1,1);


for ii = 1: nPBEs2Plot
    
    currPBE = selectPBEs(ii);
    firings2 = zeros(size(binnedPBEs{currPBE, 1}));
    for jj = 1:size(binnedPBEs{currPBE, 1}, 2)
        firings2(:, jj) = conv(binnedPBEs{currPBE, 1}(:, jj), gaussianWin, 'same');
    end
    firings  = zeros(size(firings2));
    firings(firings2 > 0) = 1;
    
    h1 = subplot(20 ,10, firstRow(ii));
    imagesc(firings);
    colormap(h1, flipud(gray))
%     caxis([0 10])
%     title(sprintf('event %d, sc=%.2f', currPBE, HMMSeqScores.data(currPBE)))
    
    text(1,10, sprintf('event %d, sc=%.2f', currPBE, HMMSeqScores.data(currPBE)), 'fontsize', 6, 'color', 'b')

    h2 = subplot(20, 10, secondRow(ii));
    imagesc(statesProbDists{currPBE})
    colormap(h2, flipud(bone))
    
end


%% plot example PBEs showing transition between states within them

examplePBEs = [ 384 2222 1959 2148 3121 2367 779 1993]; %99 2678 2240 1134 2893
examplePBEs = sort(examplePBEs, 'ascend'); 

cnctPBEs      = cell2mat(binnedPBEs(examplePBEs, 1)');
cnctStateProb = cell2mat(statesProbDists(examplePBEs)');
cnctStateProb = repelem(cnctStateProb, 1, 20);

figure; 



a1 = subplot('Position',[.1 .3 .8 .3]);
imagesc(cnctStateProb); colormap(flipud(bone))

summedLen= 0; 
for ii = 1:numel(examplePBEs)
    loc = summedLen + size(statesProbDists{examplePBEs(ii)},2)*20;
    line([loc loc], [0 size(statesProbDists{examplePBEs(ii)},1)], 'linestyle', '--', 'linewidth', 2, 'color', 'b')
    summedLen = summedLen + size(statesProbDists{examplePBEs(ii)},2)*20;
end

ylim([0 30])

a2 = subplot('Position',[.1 .6 .8 .3]);

for unit = 1 : size(cnctPBEs,1)

    spikes = find(cnctPBEs(unit,:));

    if ~isempty(spikes)
        for spk = 1 : length(spikes)
            plot([spikes(spk), spikes(spk)], [unit-1.25, unit+1.25], 'color', 'k', 'linewidth', 2) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
            hold on
        end
    end

end

summedLen= 0; 
for ii = 1:numel(examplePBEs)
    loc = summedLen+size(binnedPBEs{examplePBEs(ii), 1},2);
    line([loc loc], [0 size(binnedPBEs{examplePBEs(ii), 1},1)], 'linestyle', '--', 'linewidth', 2, 'color', 'b')
    summedLen = summedLen + size(binnedPBEs{examplePBEs(ii), 1},2);
end


set(a2,'YTick',[],'XTick',[]);

xlim([0 size(cnctPBEs, 2)])
ylim([0 30])

linkaxes([a1 a2], 'x')
%% Functions


function [HMMprctile, HMMprctile_pts, statesProbDists] = HMMCrossValidAnalysis(eventsBinnedfiring, nStates, nShuffles, nFolds, directory, pbePeriod, sm)


subDirectory = fullfile(directory, pbePeriod);
mkdir(subDirectory)


% randomizing the temporal order of the PBEs
rndgen                     = randperm(size(eventsBinnedfiring, 1)); 
eventsBinnedfiring_rnd     = eventsBinnedfiring(rndgen, :); 



% train cross-validation models build on 80% of PBEs
[CVtransmats, CVlambdas, testEvts] = trainCVmodels(eventsBinnedfiring_rnd, nStates, nFolds); 

% calculating congruence of remaining PBEs with the cross-validation models
[HMMprctile, statesProbDists]      = HMMcongruence_crossvalid(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'actual', subDirectory, sm);


% back to the original order
HMMprctile(rndgen)      = HMMprctile; 
statesProbDists(rndgen) = statesProbDists;


% calculating congruence of shuffled PBEs (pooled time swapped) with the
% cross-validation models (the same models trained on actual PBEs)
HMMprctile_pts    = HMMcongruence_crossvalid_pts(eventsBinnedfiring_rnd, testEvts, CVtransmats, CVlambdas, nShuffles, 'pooledTimeSwap', subDirectory, sm); % pts: pooled time swap

% back to the original order
HMMprctile_pts    = HMMprctile_pts(rndgen);


end