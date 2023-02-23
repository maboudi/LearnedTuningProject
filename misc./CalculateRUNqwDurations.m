clear; clc; close all

addpath(genpath('/home/kouroshmaboudi/Documents/HMM_project/workingCodes'))

currDir = '/home/kouroshmaboudi/Documents/HMM_project';
cd(currDir)


%% loading session data

rats = {'Roy','Ted', 'Kevin'};
allsessionNumbers = {[1 2 3], [1 2 3], 1};


for rr = 2%:length(rats) % loop for all the rats and sessions

    
rat = rats{rr};

sessionNumbers = allsessionNumbers{rr};


for sessionNumber = 1%sessionNumbers

    
    
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


runQuietWakingPeriod = behavior.list(find(behavior.list(:,2) > behavior.time(2,1) & behavior.list(:,1) < behavior.time(2,2) & ismember(behavior.list(:,3), [1,3])), 1:2);

QWdur = runQuietWakingPeriod(:,2) - runQuietWakingPeriod(:,1);


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


% %%% see if we want to partition the data to halves (here based on the time gaps between the laps)
% 
% gapbtwLaps = laps(2:end, 1) - laps(1:end-1, 2);
% 
% figure; 
% plot(2:length(laps), gapbtwLaps)

% breakLap = input('\n Enter the lap number for truncating the data\n(If no partitioning is needed just press enter) \n');  


breakLap = [];

if ~isempty(breakLap)
    firstORsecond = input('\nselect a part of data to process (one or two)\n');
else
    firstORsecond = [];
end



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


%% RUN


fileinfo.tbegin = behavior.time(2,1);
fileinfo.tend = behavior.time(2,2);


exclude = bvrTimeList(ismember(bvrState, [2 4]),:); % nrem=1, rem=2, quiet=3, wake=4

primaryPBEs_run = PBPeriods(mySpikes, fileinfo, [], time_resolution, threshZ, exclude);

[eventsBinnedfiringRUN, secondaryPBEs_run, acceptedprmEvts_run] = finalBinningResult(primaryPBEs_run, mySpikes, fileinfo);

% savePBEs(primaryPBEs_run, secondaryPBEs_run, acceptedprmEvts_run, 'RUN', fileinfo)

% [poissonEventsBinnedfiringRUN, timeSwapEventsBinnedfiringRUN, temporalEventsBinnedfiringRUN] = genSurrogates(eventsBinnedfiringRUN);



%% categorizing the PBEs based on the quiet waking period duration

qwind = zeros(size(secondaryPBEs_run, 1), 1);

for pbe = 1:size(secondaryPBEs_run, 1) 
    
    qwind(pbe) = find(secondaryPBEs_run(pbe, 3) > runQuietWakingPeriod(:,1)  &  secondaryPBEs_run(pbe, 3) < runQuietWakingPeriod(:,2)); 
    
end

PBE_surrQWdur = QWdur(qwind)/1e6;


% loading the sequence scores (RUN PBEs  given the PRE model)

load(['/home/kouroshmaboudi/Documents/HMM_project/' sessionName2 '/HMM/sequenceDetection/congruence(wi_and_bw)/' sessionName2 '_congruence.mat'], 'HMMprctile_PRE_RUN')


% correlation between the sequence score and the surrounding quiet waking
% period


[scoreQWcorr, scoreQWPval] = corr(PBE_surrQWdur, HMMprctile_PRE_RUN);



% fit a line to scatter
coeffs = polyfit(PBE_surrQWdur, HMMprctile_PRE_RUN, 1);
% Get fitted values
fittedX = linspace(min(PBE_surrQWdur), max(PBE_surrQWdur), 200);
fittedY = polyval(coeffs, fittedX);

fileBase = [currDir '/' fileinfo.name '/HMM/sequenceDetection/' ];


figure 
set(gcf, 'position', [-1809 458 653 514])

plot(PBE_surrQWdur, HMMprctile_PRE_RUN, '.k', 'markersize', 10)

hold on
% plot the fitted line
plot(fittedX, fittedY, 'r-', 'LineWidth', 3);

xlabel('surround QW duration(sec)', 'fontsize', 16)
ylabel('Sequence score given PRE of RUN PBEs', 'fontsize', 16)
title(fileinfo.name, 'fontsize', 16)

save([fileBase 'sequenceScorevsQWdur.mat'], 'PBE_surrQWdur', 'scoreQWcorr', 'scoreQWPval')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[fileBase 'sequenceScorevsQWdur'],'-dpdf','-r0')
print(gcf,[fileBase 'sequenceScorevsQWdur'],'-dpng','-r0')




end

end





% functions

function [eventsBinnedfiring, secondaryPBEs, idx_acceptedEvents] = finalBinningResult(pbEvents2, mySpikes, fileinfo)

%% binnig the spikes within an event and qualifying the event based on number of active units and length

%%% Pre-processing the population burst events


%%% Calculating the number of firing pyramidal units within each event

qclus = [1 2 3]; % pyramidal and interneurons

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
