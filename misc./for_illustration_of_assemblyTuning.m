function basicProcedure_GrosmarkOriginal(sessionNumber)


% addpath(genpath('/nfs/turbo/umms-kdiba/Kourosh/ReplayPreplayAnalyses'))

% currDir = '/nfs/turbo/umms-kdiba/Kourosh/Grosmark_originalClusters';

currDir = '/home/kouroshmaboudi/Documents/HMM_project/Grosmark_originalClusters';
cd(currDir)


%% loading session data

% VarList = {'Spikes','Epochs','Position'};

% note that the spikes within the .mat file are only from the sorted units,
% for MUA I need to load the clu and res files again to consider all the
% spikes

dirresults = dir(currDir);

noSessions = 8;
sessionNames = cell(noSessions, 1);


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100

mazeLimits = [-inf inf -inf 13 ;... % Achilles_10252013
              -inf inf -inf inf;... % Achilles_11012013
              -inf inf -140 inf;... % Buddy_06272013
              -inf inf -inf inf;... % Cicero_09012014
              -inf inf -inf inf;... % Cicero_09102014
              -inf inf -inf inf;... % Cicero_09172014
              -inf inf -inf inf;... % Gatsby_08022013
              -20  inf -inf 110];   % Gatsby_08282013

mazeShapes = {'linear', 'circular', 'linear', 'linear', 'circular', 'linear', 'linear', 'circular'};

    
    
sessionName = dirresults(sessionNumber + 2).name;

load(fullfile(currDir, 'NoveltySessInfoMatFiles', [sessionName '_sessInfo.mat']))


spikes = sessInfo.Spikes;
behavior = sessInfo.Epochs;
position = sessInfo.Position;


behavior.time = [behavior.PREEpoch; behavior.MazeEpoch; behavior.POSTEpoch];

speed = calspeed(position); % animal's velocity





%% sessioninfo

fileBase = fullfile(currDir, sessionName, sessionName);

fileinfo = struct('name', sessionName, 'animal', sessionName(1:strfind(sessionName, '_')-1), 'xyt', [], 'tbegin', [], 'tend', [], 'Fs', [], 'lfpSampleRate', [], 'nCh', []); 


mainDir = fullfile(currDir, ['analysesResults_' datestr(now, 'dd-mmm-yyyy')], fileinfo.name);
mkdir(mainDir)

% Load the recording configurations from the xml file

Par = LoadXml(fileBase);
% fileinfo.lfpSampleRate = Par.lfpSampleRate;

fileinfo.nCh = Par.nChannels;
fileinfo.Fs = 1; % Since the spike timestamps are already in seconds we don't need to use the original sampling rate of recording
fileinfo.lfpSampleRate = 1;



% RUN high power theta periods

thetaPeriods = importEVTfile(fileBase); % make sure that we know details about procedures and parameters which were used for the detection of these periods 


% Position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

tpos = position.TimeStamps';


% We need to exclude long periods during which the positional data is
% missing for PBE detection (if we are not using theta)


temp   = [0; diff(isnan(xpos))]; 
starts = find(temp == 1);
ends   = find(temp == -1);

if ends(1) < starts(1)
    starts = [1; starts];
end

if starts(end) > ends(end)
    ends = [ends; length(xpos)];
end

nanPeriods = [tpos(starts) tpos(ends)];
nanDurs    = diff(nanPeriods')'; 
nanPeriods = nanPeriods(nanDurs > 5, :);


fileinfo.xyt = [xpos ypos position.TimeStamps']; 

figure; plot(xpos, ypos, '.', 'markersize', 3)



%% Based on the lookup table exclude the positions outside the track boundaries

fileinfo.xyt(find(ypos > mazeLimits(sessionNumber, 4) | ypos < mazeLimits(sessionNumber, 3) | xpos > mazeLimits(sessionNumber, 2) | xpos < mazeLimits(sessionNumber, 1)), 1:2) = NaN;


linearPos = linearizePosition(fileinfo, behavior, mazeShapes{sessionNumber});

fileinfo.xyt2(:, 1) = linearPos;
fileinfo.xyt2(:, 2) = fileinfo.xyt(:, 3);


% % refining the boundaries of RUN period
% 
% firstTimePnt = find(~isnan(linearPos), 1, 'first');
% 
% fileinfo.xyt2(:, 1) = linearPos(firstTimePnt: end); 
% fileinfo.xyt2(:, 2) = fileinfo.xyt(firstTimePnt: end, 3);
% 
% 
% % The following line of code was added to remove the Nans introduced
% % because of excluding the part before the animal was put on the track
% 
% behavior.time(2,1) = behavior.time(2,1) + firstTimePnt * diff(fileinfo.xyt2(1:2,2));


% % updated
if ismember(sessionNumber, [2 5 8]) 
    direction = 'uni'; 
    [laps, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
else
    direction = 'bi';
    [lapsStruct, turningPeriods] = calculateLapTimings(fileinfo, speed, direction, mainDir); 
    
    if length(lapsStruct.RL) > length(lapsStruct.LR)
       lapsStruct.RL(1,:) = [];
%        behavior.MazeEpoch(1,:) = lapsStruct.LR(1,:);
    end
    
    totNumLaps = size(lapsStruct.RL, 1) + size(lapsStruct.LR, 1);
    laps = zeros(totNumLaps, 2);
    laps(1:2:totNumLaps, :)  = lapsStruct.LR;
    laps(2:2:totNumLaps, :)  = lapsStruct.RL;
end

laps(:, 3) = 1:size(laps, 1); 

% %


fileinfo.xyt2(:, 3) = zeros(size(fileinfo.xyt2(:, 1))); % labeling the postion samples with the calculated laps (if not part of any lap the label is zero)

for ii = 1: length(laps)

   idx =  find(fileinfo.xyt2(:, 2) > laps(ii, 1) & fileinfo.xyt2(:, 2) < laps(ii, 2));
   fileinfo.xyt2(idx, 3) = laps(ii, 3);
          
end

runSpeedThresh = 10; % cm/s

%% Behavioral states

% nrem = 1, Drowsy = 3.5; rem = 2, Intermediate = 1.5, wake = 4


% Based on the description of the dataset on CRCNS, the Intermediate could
% be considered a different type of NREM, characterized by high spindle (12-20 Hz) power and low movement/EMG (but considerbale change in respect to NREM).
% So, we can include them as NREM in our analyses. 
% Drowsy states happens in transition from wake to NREM, REM to NREM or within the NREM periods,
% charachterized by low overall spectral power. Can we consider them as NREM???


bvrTimeList = [behavior.Wake; behavior.Drowsy; behavior.NREM; behavior.Intermediate; behavior.REM];
bvrState = [4   *   ones(length(behavior.Wake), 1); ... % including both the active and quiet wake periods
            3 *   ones(length(behavior.Drowsy), 1); ... % considering Drowsy as a part of NREM, for now 
            1   *   ones(length(behavior.NREM), 1); ... 
            1   *   ones(length(behavior.Intermediate), 1); ... % considering Intermediate as a part of NREM, for now 
            2   *   ones(length(behavior.REM), 1)];



%% making a new spike data just to conform with the format of Kamran's data


qual2consider = 'all';

spikeStruct = spikeBehaviorAnalysis3(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo);



% We might use multiunit (any detected spike) as well for defining the population burst events
% REMOVE CLUSTER #0 as it might contain artifacts 

numberofShanks = length(dir([fileBase '.res.*']));

MUA.t = [];
for shank = 1:numberofShanks
    
    resContent  = load([fileBase '.res.' num2str(shank)]);
    
    cluContent  = load([fileBase '.clu.' num2str(shank)]);
    clusterIDs  = cluContent(2:end);
    
    spikes2Add = resContent(clusterIDs ~= 0);
    
    MUA.t = [MUA.t; spikes2Add]; % the first is the total number of clusters detected for the shank
    
end


MUA.t = sort(MUA.t, 'ascend');
% The originla sampling frequency is 20000

MUA.t = MUA.t./20000;



%% Spatial tuning of the units


close all

subfolder = fullfile(mainDir, 'PlaceFields');
mkdir(subfolder)

% decide whether to exclude the turning periods from calculations


%%% 2D spatial tuning

% This section needs to be updated with changes similar to what I did for
% spatialTuning_1D
% 
% if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular
%     spatialTunings2D = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'uni', 1, runSpeedThresh, fileinfo.Fs, subfolder);
% else
%     spatialTunings2D_LR = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'LR', 1, runSpeedThresh, fileinfo.Fs, subfolder);
%     spatialTunings2D_RL = spatialTuning_2D(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'RL', 1, runSpeedThresh, fileinfo.Fs, subfolder);
% end


%%% 1D spatial tuning: using linearized position

% the place fields are calculated separately for the left and right
% direction


if ismember(sessionNumber, [2 5 8]) % sessions 2, 5, and 8 are circular

    [spatialTunings, PF_sorted, runTemplate, spatialInfo, conslapsRatio, diffWithAvg] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, thetaPeriods, [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings, [], spatialInfo, [], conslapsRatio, [], behavior, [], fileinfo, subfolder);
else

    [spatialTunings_LR, PF_sorted_LR, runTemplate_LR,  spatialInfo_LR, conslapsRatio_LR, diffWithAvg_LR] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'LR', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
    [spatialTunings_RL, PF_sorted_RL, runTemplate_RL,  spatialInfo_RL, conslapsRatio_RL, diffWithAvg_RL] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'RL', 2, runSpeedThresh, [], fileinfo.Fs, subfolder);
%     firingLapbyLap(spikeStruct, spatialTunings_LR, spatialTunings_RL, spatialInfo_LR, spatialInfo_RL, conslapsRatio_LR, conslapsRatio_RL, behavior, [], fileinfo, subfolder);

    [spatialTunings_biDir, PF_sorted_biDir, runTemplate_biDir,  spatialInfo_biDir, conslapsRatio_biDir, diffWithAvg_biDir, posbincenters] = spatialTuning_1D_tempModifications(spikeStruct, [1 2 3], fileinfo, behavior, [], [], speed, 'uni', 2, runSpeedThresh, [], fileinfo.Fs, subfolder, ccl);
end

save(fullfile(subfolder, 'biDirectional.mat'), 'spatialTunings_biDir', 'PF_sorted_biDir', 'runTemplate_biDir', 'spatialInfo_biDir', 'conslapsRatio_biDir', 'diffWithAvg_biDir')
% close all




%% Ripple Detection

% After visualizing deteceted events using different metrics (MUA and
% ripple ampplitude) on Neuroscope, I realized that ripple is more reliable


best_channels = BigRippleChannels(fileBase, fileinfo);
fileinfo.RippleChannels = best_channels(1);
%%
threshSD = 3;
[ripplePeriods, rippleLFP, ripplePower] = RippleDetect(fileBase, fileinfo, threshSD); % the results in seconds

%%

% doing some small experiments ..

poi = [24259 24269];
qclus = [1 2 3];

binned_poi = finalBinningResult(poi, spikeStruct, qclus, fileinfo); 

fileinfo.lfpSampleRate = 1;

fileinfo.tbegin = poi(1); 
fileinfo.tend   = poi(2);

time_resolution = 0.001;
threshZ = 3;

PBEs = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, [], 0);
binnedPBEs = finalBinningResult(PBEs, spikeStruct, qclus, fileinfo);

currRipple = rippleLFP(poi(1)*1250:poi(2)*1250);
PBEs = (PBEs(:, 1:3) - poi(1))*1000;


figure;

ax(1) = subplot(5,1,1);
hold on
plot(currRipple, 'k')

for ii = 1:size(PBEs, 1)
    
   patch([PBEs(ii, 1) PBEs(ii, 2) PBEs(ii, 2) PBEs(ii, 1)]*1.25, [-1000 -1000 1000 1000], [0, 113, 178]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.1) 
    
end


xlim([0 length(currRipple)])

set(ax(1),'XColor', 'none','YColor','none')



ax(2) = subplot(5,1,[2:5]);
hold on
currEvent  = binned_poi{1, 1};

spikesLR = currEvent(runTemplate_biDir, :); %%% for raster of place cells in each direction

cl = colormap('lines');
cl = cl(randi(64, size(spikesLR,1),1), :);
for unit = 1 : size(spikesLR,1)
                
    spikes = find(spikesLR(unit,:));

    if ~isempty(spikes)
        for spk = 1 : length(spikes)
            plot([spikes(spk), spikes(spk)], [1.5*unit-1.25, 1.5*unit+1.25], 'color', cl(unit, :), 'linewidth', 1.5) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
            hold on
        end
    end

end

xlim([0 size(spikesLR, 2)])
ylim([0 1.5*length(runTemplate_LR)+3])

% linkaxes(ax, 'x')
set(ax(2),'XColor', 'none','YColor','none')


for ii = 1:size(PBEs, 1)
    
   patch([PBEs(ii, 1) PBEs(ii, 2) PBEs(ii, 2) PBEs(ii, 1)], [0 0 1.5*length(runTemplate_LR)+3 1.5*length(runTemplate_LR)+3], [0, 113, 178]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.1) 
    
end

line([1 1000], [0 0], 'linewidth', 3, 'color', 'k')

xl=get(gca,'XLim');
yl=get(gca,'YLim');
ht = text(0.5*xl(1)+0.5*xl(2),0.6*yl(1)+0.4*yl(2),'My text');
set(ht,'Rotation',90)
set(ht,'FontSize',18)

ccl = cl;

%% 


spatialTunings_biDir(~spatialTunings_biDir) = 1e-4;

currEvent = binnedPBEs{ii, 1};


selectedUnits = find(ismember(1:size(binnedPBEs{3,2}), runTemplate_biDir(1:end)));
temp = find(sum(currEvent(selectedUnits, :), 2));
selectedUnits = selectedUnits(temp);

spatialTunings_biDir2 = spatialTunings_biDir(selectedUnits, :);


figure;

len = zeros(size(PBEs,1), 1);
ii = 3;

subplot(2,2,1)

hold on

len(ii) = size(currEvent, 2);

spikesLR = currEvent(selectedUnits, :); %%% for raster of place cells in each direction

for unit = 1 : size(spikesLR,1)
                
    spikes = find(spikesLR(unit,:));

    if ~isempty(spikes)
        for spk = 1 : length(spikes)
            plot([spikes(spk), spikes(spk)], [3*unit-1.25, 3*unit+1.25], 'color', 'k', 'linewidth', 1.5) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
            hold on
        end
    end

end

xlim([0 size(spikesLR, 2)])
ylim([0 3*length(selectedUnits)+3])
set(gca, 'box', 'on')
xticks([])
yticks([])
xticklabels([])
yticklabels([])

yl = ylim;

for jj = 1:floor(len(ii)/20)
   
    line([jj*20 jj*20], [yl(1) yl(2)], 'linestyle',':', 'linewidth', 1, 'color', 'k')
    
end

title('Spike trains', 'fontsize', 14)
ylabel('Units', 'fontsize', 14)
xlabel('time', 'fontsize', 14)

subplot(2,2,2)

spatialTunings_biDir3 = spatialTunings_biDir2 ./ repmat(max(spatialTunings_biDir2, [], 2), [1 size(spatialTunings_biDir2, 2)]);

hold on

for jj = 1:size(spatialTunings_biDir3, 1)

    fill([posbincenters fliplr(posbincenters)], [jj+0.9*spatialTunings_biDir3(jj, :) fliplr(jj*ones(size(spatialTunings_biDir3(jj, :))))], 'k','LineStyle','none')

    plot(posbincenters, jj+0.9*spatialTunings_biDir3(jj, :),'color', 'k','linewidth', 0.5);

    alpha(0.5)
end

ylim([0.5 43])
xlim([0 posbincenters(end)])
set(gca,'YColor','none')

title('Spatial tunings', 'fontsize', 14)
xlabel('position(cm)', 'fontsize', 14)


subplot(2,2,3)


currEvent_20  = binnedPBEs{3, 2};
currEvent_20  = currEvent_20(selectedUnits, :);

nPosBins = size(spatialTunings_biDir2, 2);

binDur = 0.02;

postProb = baysDecoder(currEvent_20, spatialTunings_biDir2, binDur);  
postProb = postProb./repmat(sum(postProb, 1), [nPosBins, 1]); 


currPPM   = postProb;
maxProb   = max(currPPM(:));
minProb   = min(currPPM(:));

Clim_low  = minProb;
Clim_high = minProb + 0.75 *(maxProb - minProb);


imagesc(1:size(currPPM, 2), 2*(1:size(currPPM, 1)), currPPM, [Clim_low Clim_high])
colormap('hot')
set(gca,'YDir','normal')

% xticks([])
% yticks([])
% xticklabels([])
% yticklabels([])

ylabel('Position(cm)', 'fontsize', 14)
xlabel('time(bins)', 'fontsize', 14)

oldPos = get(gca, 'position');
set(gca, 'position', [oldPos(1) oldPos(2) oldPos(3) 0.7*oldPos(4)])



%% 


postProb = baysDecoder(currEvent_20, spatialTunings_biDir2, binDur);  
postProb = postProb./repmat(sum(postProb, 1), [nPosBins, 1]); 


[nPositionBins, nTimeBins] = size(postProb);

fitLineDictionary = generateLineDictionary(nPositionBins, nTimeBins);


[replayScore, weightedCorr,~, ~, begPosition, endPosition] = replayevaluation_polar2(postProb, fitLineDictionary);



tuningshuffles = tuningsUnitIDshuffle_normalize(spatialTunings_biDir2, [], 1:size(spatialTunings_biDir2, 1), 1000, 'same-peak');


replayScore_shuffle = zeros(1000, 2);
weightedCorr_shuffle = zeros(1000, 2);

% time swap
for ii = 1:1000
   
    
    idx = randperm(size(postProb, 2));
    postProb_shuffle = postProb(:, idx);
    
    [replayScore_shuffle(ii,1), temp] = replayevaluation_polar2(postProb_shuffle, fitLineDictionary);
     weightedCorr_shuffle(ii,1) = abs(temp);

end


% unit ID

for ii = 1:1000
        
    postProb_shuffle = baysDecoder(currEvent_20, tuningshuffles(:,:,ii), binDur);  
    postProb_shuffle = postProb_shuffle./repmat(sum(postProb_shuffle, 1), [nPosBins, 1]); 

    
    [replayScore_shuffle(ii,2), temp] = replayevaluation_polar2(postProb_shuffle, fitLineDictionary);
    weightedCorr_shuffle(ii,2) = abs(temp);
end






bins_1 = linspace(min(replayScore_shuffle(:))-0.1, max(replayScore_shuffle(:))+0.1, 30);

[h1] = hist(replayScore_shuffle(:,1), bins_1); h1 = h1/sum(h1);
[h3] = hist(replayScore_shuffle(:,2), bins_1); h3 = h3/sum(h3);


bins_2 = linspace(min(weightedCorr_shuffle(:))-0.1, max(weightedCorr_shuffle(:))+0.1, 30);

[h2] = hist(weightedCorr_shuffle(:,1), bins_2); h2 = h2/sum(h2);
[h4] = hist(weightedCorr_shuffle(:,2), bins_2); h4 = h4/sum(h4);



figure;

subplot(2,1,1)
hold on

set(gca, 'fontsize', 12)
aa = plot(bins_1, h1, 'color', [204 0 0]/255, 'linewidth', 2);
bb = plot(bins_1, h3, 'color', [204 0 0]/255, 'linewidth', 2);
bb.Color(4) = 0.5;


line([replayScore replayScore], [0 0.2], 'color', [0.3 0.3 0.3], 'linewidth', 3)

legend([aa bb], 'time swap', 'unit-ID shuffle', 'fontsize', 12, 'location', 'best')
xlabel('Radon integral', 'fontsize', 14)
ylabel('Probability', 'fontsize', 14)


subplot(2,1,2)
hold on

set(gca, 'fontsize', 12)

aa = plot(bins_2, h2, 'color', [0 153 0]/255, 'linewidth', 2);
bb = plot(bins_2, h4, 'color', [0 153 0]/255, 'linewidth', 2);
bb.Color(4) = 0.5;

line([abs(weightedCorr) abs(weightedCorr)], [0 0.1], 'color', [0.3 0.3 0.3], 'linewidth', 3)
legend([aa bb], 'time swap', 'unit-ID shuffle', 'fontsize', 12 , 'location', 'best')
xlabel('Abs weighted correlation', 'fontsize', 14)
ylabel('Probability', 'fontsize', 14)





%%
figure;

hold on


currPPM   = postProb;
maxProb   = max(currPPM(:));
minProb   = min(currPPM(:));

Clim_low  = minProb;
Clim_high = minProb + 0.75 *(maxProb - minProb);


imagesc(1:size(currPPM, 2), 2*(1:size(currPPM, 1)), currPPM, [Clim_low Clim_high])

line([begPosition(1) endPosition(1)], 2*[begPosition(2) endPosition(2)], 'color', 'w', 'linewidth', 2)
line([begPosition(1) endPosition(1)], 2*[begPosition(2)+7 endPosition(2)+7], 'color', 'w', 'linestyle', ':', 'linewidth', 2)
line([begPosition(1) endPosition(1)], 2*[begPosition(2)-7 endPosition(2)-7], 'color', 'w', 'linestyle', ':', 'linewidth', 2)
colormap('hot')
set(gca,'YDir','normal')

xlim([0.5 nTimeBins+0.5])
ylim([0 2*nPositionBins])

% xticks([])
% yticks([])
% xticklabels([])
% yticklabels([])

ylabel('Position(cm)', 'fontsize', 14)
xlabel('time(bins)', 'fontsize', 14)



% 
% 
% if ii > 1
% ref = get(ax(1), 'position');
% oldPos = get(ax(ii), 'position');
% 
% set(ax(ii), 'position', [oldPos(1:2) ref(3)*len(ii)/len(1) oldPos(4)])
% 
% end

% end



%%

cnt = 1;
for ii = 1:size(PBEs,1)
   
   idx = find(binnedPBEs{ii, 2}(61, :));
   
   if ~isempty(idx)
   for jj = 1:length(idx) 
      
       indivBin{cnt} = binnedPBEs{ii, 1}(:, (idx(jj)-1)*20+1 : idx(jj)*20);
       indivBin_20{cnt} = binnedPBEs{ii, 2}(:, idx(jj));
       cnt = cnt+1;
   end
   end
end


figure;

for ii = 1:length(indivBin)
   
    ax(ii) = subplot(1, length(indivBin), ii);
    
    hold on
    currEvent  = indivBin{ii};
    len(ii) = size(currEvent, 2);

    spikesLR = currEvent(runTemplate_biDir, :); %%% for raster of place cells in each direction

    for unit = 1 : size(spikesLR,1)

        spikes = find(spikesLR(unit,:));

        if ~isempty(spikes)
            for spk = 1 : length(spikes)
                plot([spikes(spk), spikes(spk)], [1.5*unit-1.25, 1.5*unit+1.25], 'color', ccl(unit, :), 'linewidth', 1.5) % colorset(ceil(tempate_LR(unit)*128/noUnits), :)
                hold on
            end
        end

    end
    set(gca, 'box', 'on')
    xticks([])
    yticks([])
    xticklabels([])
    yticklabels([])
    xlim([-1 21])
    
end

linkaxes(ax, 'xy')


%% decode the positions in the example bins

nUnits = size(binnedPBEs{1,1},1);

currUnit = 61;
otherUnits = setdiff(1:nUnits, currUnit);

spatialTunings_LR(~spatialTunings_LR) = 1e-4;
spatialTunings_RL(~spatialTunings_RL) = 1e-4; 

nPosBins = size(spatialTunings_RL, 2);

binDur = 0.02;

decoding = zeros(nPosBins, length(indivBin));
unitsFiring = zeros(length(indivBin), 1);
for ii = 1:length(indivBin)
    
    currBin = indivBin_20{ii};

    postprobRL = baysDecoder(currBin(otherUnits), spatialTunings_LR(otherUnits, :), binDur);  
    postprobLR = baysDecoder(currBin(otherUnits), spatialTunings_RL(otherUnits, :), binDur);  


    %%% marginalized over directions %%%

    posteriorProbMatrix = postprobRL + postprobLR;
    decoding(:,ii) = posteriorProbMatrix./repmat(sum(posteriorProbMatrix, 1), [nPosBins, 1]); 
    
    unitsFiring(ii) = currBin(currUnit);
    
end





currPPM   = decoding;
maxProb   = max(currPPM(:));
minProb   = min(currPPM(:));

Clim_low  = minProb;
Clim_high = minProb + 0.75 *(maxProb - minProb);



figure;


for ii = 1:size(decoding, 2)
    
   ax(ii) = subplot(1, length(indivBin), ii);

    imagesc(decoding(:, ii), [Clim_low Clim_high])
    set(gca,'YDir','normal')
    colormap(ax(ii), 'hot')
    
    xticks([])
    yticks([])
    xticklabels([])
    yticklabels([])
 
end




%%

assemblyTunings = calAssemblyTunings(decoding, unitsFiring');


currPPM   = decoding;
maxProb   = max(currPPM(:));
minProb   = min(currPPM(:));

Clim_low  = minProb;
Clim_high = minProb + 0.75 *(maxProb - minProb);



figure;

imagesc(assemblyTunings, [Clim_low Clim_high])
set(gca,'YDir','normal')
colormap('hot')

xticks([])
yticks([])
xticklabels([])
yticklabels([])



figure;
imagesc(spatialTunings_biDir(61, :))
colormap('jet')

xticks([])
yticks([])
xticklabels([])
yticklabels([])














%% PBEs


%% detemining population burst periods


time_resolution = 0.001; % in secondprint(gcf, [FileBase fileinfo.name '_congruence'], '-dpng')
threshZ = 3; % sdf with 3 std deviation above the mean


% baseFile = fullfile(unsortedDataDir, sessionName, sessionName);
% 
% thetaPeriods = load([baseFile '.theta.1']);
% thetaPeriods = floor(thetaPeriods./fileinfo.lfpSampleRate);


fileinfo.lfpSampleRate = 1;

fileinfo.tbegin = behavior.time(1,1); 
fileinfo.tend = behavior.time(3,2);


subfolder = fullfile(mainDir, 'PBEs', 'wholeSession');
mkdir(subfolder)


exclude = bvrTimeList(ismember(bvrState, [2 4]), :); % nrem=1, rem=2, wake=4


velocityFilter = 1;
[primaryPBEs, sdat] = PBPeriods(spikeStruct, fileinfo, [], time_resolution, threshZ, exclude, velocityFilter);


qclus = [1 2 3];

[binnedPBEs, secondaryPBEs, nFiringUnits, PBElength] = finalBinningResult(primaryPBEs, spikeStruct, qclus, fileinfo); 
PBErippleIdx = ifContainRipples(secondaryPBEs, ripplePeriods);

secondaryPBEs(:, 5) = PBErippleIdx;

nPBEs = size(binnedPBEs, 1);


secondaryPBEs = [secondaryPBEs zeros(nPBEs, 4)];

for ii= 1:nPBEs
    
    pbeCenter = secondaryPBEs(ii, 3);
    
    boutInd        = find(bvrTimeList(:,1) < pbeCenter & bvrTimeList(:,2) > pbeCenter, 1, 'first');
    boutBrainState = bvrState(boutInd);
    
    secondaryPBEs(ii, 5 + boutBrainState) = 1;
    
end


baseStruct = struct('data', [], 'p', [], 'ts', [], 'pts', []);


PREbinnedPBEs       = baseStruct;
PREidx              = find(secondaryPBEs(:, 1) > behavior.time(1,1) & secondaryPBEs(:, 2) < behavior.time(1,2));
PREbinnedPBEs.data  = binnedPBEs(PREidx, :);
secondaryPBEs_PRE   = secondaryPBEs(PREidx, :);
PREbinnedPBEs       = genSurrogates(PREbinnedPBEs);


RUNbinnedPBEs       = baseStruct;
RUNidx              = find(secondaryPBEs(:, 1) > behavior.time(2,1) & secondaryPBEs(:, 2) < behavior.time(2,2));
RUNbinnedPBEs.data  = binnedPBEs(RUNidx, :);
secondaryPBEs_RUN   = secondaryPBEs(RUNidx, :);
RUNbinnedPBEs       = genSurrogates(RUNbinnedPBEs);


POSTbinnedPBEs       = baseStruct;
POSTidx              = find(secondaryPBEs(:, 1) > behavior.time(3,1) & secondaryPBEs(:, 2) < behavior.time(3,2));
POSTbinnedPBEs.data  = binnedPBEs(POSTidx, :);
secondaryPBEs_POST   = secondaryPBEs(POSTidx, :);
POSTbinnedPBEs       = genSurrogates(POSTbinnedPBEs);


BayesianReplayDetection_GrosmarkReclu_nov2020


end


function speed = calspeed(position)


cov2cmfac = 100; % the position are in m in the dataset, hence positions must be multiplied with 100


% 2D position data

xpos = position.TwoDLocation(:,1)* cov2cmfac;
ypos = position.TwoDLocation(:,2)* cov2cmfac;

timepnts = position.TimeStamps';


diffx = [0; abs(diff(xpos))];
diffy = [0; abs(diff(ypos))];

difftt = [1; diff(timepnts)];


% calculation of speed in just x direction (we then can use the speed or theta for filtering the events)

velocity = sqrt(diffx.^2 + diffy.^2)./difftt; % 2D velocity 

velocity(isnan(velocity)) = interp1(timepnts(~isnan(velocity)), velocity(~isnan(velocity)), timepnts(isnan(velocity)));

% smoothing the speed

sigma = 20; %%% smoothing the speed, the positon sampling rate is around 40 Hz (0.0256 sec sampling period), so duration of the sigma is about 25.6 ms times sigma (25.6*20 ~ 512 ms)
halfwidth = 3 * sigma;


smoothwin = gausswindow(sigma, halfwidth); % smoothing kernel
velocity = conv(velocity, smoothwin, 'same'); 


speed.v = velocity;
speed.t = timepnts;


end



function [assemblyTunings, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
%     [assemblyTunings, assemblyTunings_ci, p_of_x] = calAssemblyTunings(posteriorProbMatrix, unitFirings)
    
    % spike-weighted average of posteriors divided by average of the
    % posteriors
    
    nTimeBins = numel(unitFirings);
    nPosBins  = size(posteriorProbMatrix, 1);


    weightedSummation = sum(repmat(unitFirings, [nPosBins 1]) .* posteriorProbMatrix, 2)/nTimeBins;
    p_of_x = sum(posteriorProbMatrix, 2)/nTimeBins;
    
    
    assemblyTunings = weightedSummation; 
    
    
%     sigma2 = sum(repmat(unitFirings, [nPosBins 1]) .* ((posteriorProbMatrix ./ repmat(p_of_x, [1 nTimeBins])) - assemblyTunings).^2, 2)/ nTimeBins;
%     sigma = sqrt(sigma2);
%     
%     assemblyTunings_ci = assemblyTunings - sigma/sqrt(nTimeBins) * 1.645;
    
    
end
