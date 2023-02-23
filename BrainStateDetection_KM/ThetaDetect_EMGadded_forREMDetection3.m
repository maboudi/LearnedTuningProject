function thetaPeriods = ThetaDetect_EMGadded_forREMDetection3(basePath)


% theta detection separately performed on each epoch. The reason was that
% the theta frequency band and the frequency gap between thata and its
% first harmonic is different between the sleep and maze epochs. 


if ~exist('basePath','var')
    basePath = uigetdir(cd,...
        'Which recording(s) would you like to state score?');
    if isequal(basePath,0)
        return
    end  
end


[~, sessionName] = fileparts(basePath);

outputDir = fullfile('~/Documents/NCMLproject/StateDetectionResults', sessionName);
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
else
    load(fullfile(outputDir, [sessionName, '.SleepScoreMetrics.mat']), 'SleepScoreMetrics')
    thetaChannel = SleepScoreMetrics.THchanID;
end


Par = LoadXml(fullfile(basePath, [sessionName '.xml']));
lfpSampleRate = Par.lfpSampleRate;
nCh = Par.nChannels;


% epochs
epochNames = {'pre'; 'run'; 'post'; 'remaze'};


% for Bapun sessions add the follwoing steps
 
load(fullfile(basePath, [sessionName '.epochs.mat']))

pre  = double(pre);
post = double(post);

if exist('maze', 'var')
    maze = double(maze);
elseif exist('maze1', 'var')
    maze = double(maze1);
end


epochs.pre = pre;

epochs.run(1) = pre(2);
epochs.run(2) = maze(2);

epochs.post(1) = maze(2);
epochs.post(2) = post(2);

epochs.remaze(1) = post(2);
epochs.remaze(2) = re_maze(2);



% for Grosmark's sessions
% 
% load(fullfile('/home/kouroshmaboudi/Documents/NCMLproject/assemblyTuning_finalResults/', sessionName, 'spikes', [sessionName '.spikes.mat']), 'fileInfo') % for Grosmark sessions
% 
% behavior = fileInfo.behavior;
% 
% epochs.pre  = behavior.PREEpoch;
% epochs.run  = behavior.MazeEpoch;
% epochs.post = behavior.POSTEpoch;


% % theta and outside-theta frequency bands

% Rat S
theta_f_in.pre = [5 10];
theta_f_out.pre = [1 4 10 13];

theta_f_in.run = [5 10];
theta_f_out.run = [1 5 10 12];

theta_f_in.post = [4 8];
theta_f_out.post = [1 4 10 13];

theta_f_in.remaze = [5 10];
theta_f_out.remaze = [1 4 10 12];



% loading emg and sw calculated using Buzcode

emg      = SleepScoreMetrics.EMG;
slowWave = SleepScoreMetrics.broadbandSlowWave;
timePnts = SleepScoreMetrics.t_clus;


% fileName = fullfile('~/Documents/NCMLproject/assemblyTuning_finalResults', sessionName, 'spikes', [sessionName '.brainState.mat']);
% load(fileName, 'brainState')
% 
% emg      = brainState.emg;
% slowWave = brainState.slowWave;
% timePnts = brainState.sw_emg_timePnts;



%% loading the eeg file (I should add a procedure to calculate the best theta channel)

Eeg = readmulti(fullfile(basePath, [sessionName '.eeg']), nCh, thetaChannel+1); % added +1 becasue channel numbers in Buzcode starts from 0

Window = 1; %sec - for spectrogram calculation
FreqRange = [1 100]; % Hz - range of freq. for the spectrogram
SpecWindow = 2^round(log2(Window*lfpSampleRate));% choose window length as power of two
overlapRatio = 0.75;

nFFT = SpecWindow*4;
weeg = WhitenSignal(Eeg,lfpSampleRate*2000,1);

[Pxx,f,t] = mtcsglong(weeg,nFFT, lfpSampleRate, SpecWindow, overlapRatio*SpecWindow, 2,'linear',[],FreqRange);

% t starts from zero; the start time should be adjusted to the center of
% the first window
t = t + SpecWindow/2/lfpSampleRate;



%% detecting theta states by applying HHm on theta ratio


nStates = 2; % parameters for HMM

for iepoch = 1:numel(epochNames)
    
    currEpoch = epochNames{iepoch};
    period    = epochs.(currEpoch); % period during which the theta detection is performed
            
    in_period_idx.(currEpoch) = find(t >= period(1) & t < period(2)); 
    
    if numel(theta_f_in.(currEpoch)) == 2 
        thfin = find(f > theta_f_in.(currEpoch)(1) & f < theta_f_in.(currEpoch)(2));    
    elseif numel(theta_f_in.(currEpoch)) == 4 % in case we are including the second harmonic as well
        thfin = find((f > theta_f_in.(currEpoch)(1) & f < theta_f_in.(currEpoch)(2)) | (f > theta_f_in.(currEpoch)(3) & f < theta_f_in.(currEpoch)(4)) );
    end

    if numel(theta_f_out.(currEpoch)) == 4 
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2))|(f > theta_f_out.(currEpoch)(3) & f < theta_f_out.(currEpoch)(4))); 
    elseif numel(theta_f_out.(currEpoch)) == 2
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2)));
    end

    [thratio.(currEpoch), inTh.(currEpoch), thetaPeriods.(currEpoch)] = TDRatioAuto(Pxx, nStates, t, in_period_idx.(currEpoch), thfin, thfout, Par); % , thetaPeriods.(currEpoch)
    
end



% concatenating data from different epochs

theratio_all      = [];
thetaUpIdx        = [];
% thetaPeriods_all = [];
in_period_idx_all = [];

for iepoch = 1:numel(epochNames)
    currEpoch = epochNames{iepoch};

    theratio_all      = [theratio_all; thratio.(currEpoch)];
    thetaUpIdx        = [thetaUpIdx inTh.(currEpoch)];
%     thetaPeriods_all  = [thetaPeriods_all; thetaPeriods.(currEpoch)];
    in_period_idx_all = [in_period_idx_all; in_period_idx.(currEpoch)];

end

t = t(in_period_idx_all);


%% schmidt triggering  emg, slow wave


% emg
[emgPeriods, emgUpIdx] = schnmidtTrigger(emg, timePnts, t);


% slow wave
[swPeriods, swUpIdx]   = schnmidtTrigger(slowWave, timePnts, t);


%% detecting brain states using the theta, emg and slow wave

NREMIdx = swUpIdx;
REMIdx = ~NREMIdx & ~emgUpIdx & thetaUpIdx';

WAKEIdx  = ~NREMIdx & emgUpIdx & thetaUpIdx';
QWAKEIdx = ~NREMIdx & ~thetaUpIdx';


NREMperiods  = IdxtoInt(NREMIdx, t);
REMperiods   = IdxtoInt(REMIdx, t);
WAKEperiods  = IdxtoInt(WAKEIdx, t); 
QWAKEperiods = IdxtoInt(QWAKEIdx, t);



% doing some corrections

allPeriods = [NREMperiods ones(size(NREMperiods, 1), 1); ... 
    REMperiods 2*ones(size(REMperiods, 1), 1); ...
    WAKEperiods 3*ones(size(WAKEperiods, 1), 1); ...
    QWAKEperiods 4*ones(size(QWAKEperiods, 1), 1)];

[~, sortIdx] = sort(allPeriods(:, 1), 'ascend');
allPeriods = allPeriods(sortIdx, :);

durations = diff(allPeriods(:, 1:2), [], 2);


% % short wake after REM --> REM 
% idx = find(durations(2:end) < 6 & ismember(allPeriods(2:end, 3), [3 4]) & allPeriods(1:end-1, 3) == 2);
% allPeriods(idx+1, 3) = 2;
% [allPeriods, durations] = concatSameState(allPeriods);


% wake in the middle of REM --> REM
idx = find(ismember(allPeriods(2:end-1, 3), [3 4]) & (allPeriods(1:end-2, 3) == 2 & durations(1:end-2) > 10) & (allPeriods(3:end, 3) == 2 & durations(3:end) > 10));
allPeriods(idx+1, 3) = 2;
[allPeriods, durations] = concatSameState(allPeriods);


% REM in the middle of wake --> wake
idx = find(durations(2:end-1) < 6 & allPeriods(2:end-1, 3) == 2 & (ismember(allPeriods(1:end-2, 3), [3 4]) & durations(1:end-2) > 6) & (ismember(allPeriods(3:end, 3), [3 4]) & durations(3:end) > 6));
allPeriods(idx+1, 3) = 3;
[allPeriods, durations] = concatSameState(allPeriods);


% REM that is not after a NREM or another REM --> QWake
idx = find(durations(2:end) < 10 & allPeriods(2:end, 3) == 2 & (~ismember(allPeriods(1:end-1, 3), [1 2]) & durations(1:end-1) > 10)); 
allPeriods(idx+1, 3) = 4;
[allPeriods, durations] = concatSameState(allPeriods);


% short Qwake in the middle of NREM --> NREM
idx = find(durations(2:end-1) < 6 & allPeriods(2:end-1, 3) == 4 & allPeriods(1:end-2, 3) == 1 & allPeriods(3:end, 3) == 1);
allPeriods(idx+1, 3) = 1;
allPeriods = concatSameState(allPeriods);


allPeriods(allPeriods(:, 3) == 2 & allPeriods(:, 1) < epochs.pre(2), 3) = 3; % all the detected REM epochs in PRE --> wake

allPeriods(allPeriods(:, 3) == 2 & allPeriods(:, 1) > epochs.run(1) & allPeriods(:, 1) < epochs.run(2), 3) = 3; % all the detected REM epochs during MAZE --> wake
allPeriods(allPeriods(:, 3) == 1 & allPeriods(:, 1) > epochs.run(1) & allPeriods(:, 1) < epochs.run(2), 3) = 4; % all detected NREM periods during MAZE --> Qwake 

allPeriods(allPeriods(:, 3) == 2 & allPeriods(:, 1) > epochs.remaze(1) & allPeriods(:, 1) < epochs.remaze(2), 3) = 3; % all the detected REM epochs during REMAZE --> wake
allPeriods(allPeriods(:, 3) == 1 & allPeriods(:, 1) > epochs.remaze(1) & allPeriods(:, 1) < epochs.remaze(2), 3) = 4; % all detected NREM periods during REMAZE --> Qwake 

allPeriods(allPeriods(:, 3) == 2 & allPeriods(:, 1) > epochs.post(1) & allPeriods(:, 1) < 1.55e4, 3) = 3;
allPeriods(allPeriods(:, 3) == 2 & allPeriods(:, 1) > 17133 & allPeriods(:, 1) < 18329, 3) = 3;



stateIdx = allPeriods(:, 3);
NREMperiods  = allPeriods(stateIdx == 1, :);
REMperiods   = allPeriods(stateIdx == 2, :);
WAKEperiods  = allPeriods(stateIdx == 3, :); 
QWAKEperiods = allPeriods(stateIdx == 4, :);



%% plots

Pxx = Pxx(in_period_idx_all, :);

period2plot = epochs.post;
% period2plot = [epochs.post(1) epochs.post(1)+5*60*60];

tindex = find(t>period2plot(1) & t<period2plot(2));



plotheight = 450;
plotwidth  = 900;
fontsize   = 16;

figure;
clf(gcf);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



% spectrogram
ax(1) = subplot(13, 1, 1:4);

spectrogram2Plot = log10(Pxx(tindex, 1:find(f>100, 1, 'first'), 1))';

cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
imagesc(t(tindex), f(1:find(f > 100, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 12)




% raw EEG 
ax(2) = subplot(13, 1, 5);

ts = 1/lfpSampleRate;
teeg = (period2plot(1)+ts):ts:period2plot(2);

currentEeg = Eeg((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg), 'linewidth', 1, 'color', [0 0 0 0.5]);

ylabel('raw LFP', 'FontSize', fontsize)

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'XTickLabel', [], 'YTickLabel', [], 'fontsize', 12)



% theta LFP
ax(3) = subplot(13, 1, 6:7);

forder = 256;  
forder = ceil(forder/2)*2; 

fTheta = fir1(forder,theta_f_in.pre/lfpSampleRate*2); % calculate convolution func
thetaLFP = Filter0(fTheta, Eeg);


currentEeg = thetaLFP((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg), 'linewidth', 1, 'color', [0 0 0 0.5]);

ylabel(sprintf('theta (%d - %d)', theta_f_in.pre(1), theta_f_in.pre(2)), 'fontsize', fontsize)
set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'XTickLabel', [], 'fontsize', 12)



% theta ratio
ax(4) = subplot(13, 1, 8:9);

hold on

plot(t(tindex), theratio_all(tindex), 'linewidth', 1, 'color', [0 0 1], 'DisplayName', 'theta ratio')

currNREMperiods = NREMperiods(NREMperiods(:, 1) > period2plot(1) & NREMperiods(:, 1) < period2plot(2), :);  
currREMperiods = REMperiods(REMperiods(:, 1) > period2plot(1) & REMperiods(:, 1) < period2plot(2), :);  
currWAKEperiods = WAKEperiods(WAKEperiods(:, 1) > period2plot(1) & WAKEperiods(:, 1) < period2plot(2), :);
currQWAKEperiods = QWAKEperiods(QWAKEperiods(:, 1) > period2plot(1) & QWAKEperiods(:, 1) < period2plot(2), :);

for ii = 1:size(currNREMperiods, 1)
    h1 = line([currNREMperiods(ii, 1) currNREMperiods(ii, 2)], [1 1], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'NREM');
end

for ii = 1:size(currREMperiods, 1)
    h2 = line([currREMperiods(ii, 1) currREMperiods(ii, 2)], [1.5 1.5], 'linewidth', 1, 'color', [0 1 0 1], 'DisplayName', 'REM');
end

for ii = 1:size(currWAKEperiods, 1)
    h3 = line([currWAKEperiods(ii, 1) currWAKEperiods(ii, 2)], [2 2], 'linewidth', 1, 'color', [0 0 0 1], 'DisplayName', 'WAKE');
end

for ii = 1:size(currQWAKEperiods, 1)
    h4 = line([currQWAKEperiods(ii, 1) currQWAKEperiods(ii, 2)], [2.5 2.5], 'linewidth', 1, 'color', [0 0 0 0.5], 'DisplayName', 'QWAKE');
end


ylim([-1 3])
% legend([h1, h2, h3, h4])
ylabel('theta ratio', 'FontSize', fontsize, 'Color', [0 0 1])
set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'XTickLabel', [], 'fontsize', 12)




% EMG
tindex2 = find(timePnts >period2plot(1) & timePnts<period2plot(2));
ax(5) = subplot(13, 1, 10:11);

hold on
plot(timePnts(tindex2), emg(tindex2), 'linewidth', 1, 'color', [0 0 0], 'DisplayName', 'EMG')

currEMGPeridos = emgPeriods(emgPeriods(:, 1) >= period2plot(1) & emgPeriods(:, 1) <= period2plot(2), :);
for ii = 1:size(currEMGPeridos, 1)
    h1 = line([currEMGPeridos(ii, 1) currEMGPeridos(ii, 2)], [0.8 0.8], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'hi emg periods');
end

legend(h1)
ylabel('EMG', 'FontSize', fontsize)
set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'XTickLabel', [], 'fontsize', 12)




% SWS
ax(6) = subplot(13, 1, 12:13);

hold on
plot(timePnts(tindex2), slowWave(tindex2), 'linewidth', 1, 'color', [0 0 0], 'DisplayName', 'EMG')

currSWPeridos = swPeriods(swPeriods(:, 1) >= period2plot(1) & swPeriods(:, 1) <= period2plot(2), :);
for ii = 1:size(currSWPeridos, 1)
    h1 = line([currSWPeridos(ii, 1) currSWPeridos(ii, 2)], [0.8 0.8], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'hi sws periods');
end

legend(h1)
ylabel('SWS', 'FontSize', fontsize)
set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 12)


linkaxes(ax, 'x')
xlim(period2plot)



%% choose which method is most accurate

detectionParams.theta_f_in      = theta_f_in;
detectionParams.theta_f_out     = theta_f_out;
detectionParams.winLen          = Window;
detectionParams.WinOverlapRatio = overlapRatio;
detectionParams.thetaChannel    = thetaChannel;
detectionParams.date            = date;


allPeriods = [NREMperiods ones(size(NREMperiods, 1), 1); ... 
    REMperiods 2*ones(size(REMperiods, 1), 1); ...
    WAKEperiods 3*ones(size(WAKEperiods, 1), 1); ...
    QWAKEperiods 4*ones(size(QWAKEperiods, 1), 1)];

[~, sortIdx] = sort(allPeriods(:, 1), 'ascend');

clear brainState

brainState.bouts = allPeriods(sortIdx, :);
brainState.names = {'NREM'; 'REM'; 'WAKE'; 'QWAKE'};


timePnts = t;
theratio = theratio_all;

save(fullfile(outputDir, [sessionName '.brainStateDetection_HMMtheta_EMG_SWS_SchmidtTrigger.mat']), 'brainState', 'detectionParams', 'timePnts', 'theratio') % 'inTh'


end



function [thratio, InTh, thetaPeriods, thratio_st] = TDRatioAuto(Pxx, nStates, t, in_period_idx, thfin, thfout, Par)

%automatic theta periods detection just using the thetaratio

% thfin = find(f>5 & f<11);
% thfout = find((f>1 & f<4)|(f> 12& f<15));
% thfout = find((f>1 & f<4));


thratio = log(mean(squeeze(Pxx(in_period_idx,thfin,1)),2))-log(mean(squeeze(Pxx(in_period_idx,thfout,1)),2));

gw = gausswindow(6,15);
thratio = conv(thratio, gw, 'same'); % smooth the theta delta ratio

    
% fit gaussian mixture and to HMM - experimental version .. uses only thetaratio

TheState = gausshmm(thratio,{nStates,1,0});

thratio_st = zeros(nStates, 1);

for i=1:nStates
    thratio_st(i) = mean(thratio(TheState==i));
end


%     [~, sortInd] = sort(thratio_st, 'ascend');

thetaStates    = find(thratio_st  == max(thratio_st)); % == max(thratio_st)
nonThetaStates = setdiff(1:nStates, thetaStates); 

InTh = zeros(size(TheState));


% the two states with the highest theta2delta ratio are considered as
% theta states

InTh(ismember(TheState, nonThetaStates)) = 0; 
InTh(ismember(TheState, thetaStates))    = 1;


minLen = 3;
thetaPeriods = calThetaPeriods(InTh, t(in_period_idx), minLen);
thetaPeriods = round(thetaPeriods * Par.lfpSampleRate);



end


function thetaPeriods = calThetaPeriods(insideTheta, t, minLen)


temp = [0 diff(insideTheta)]';

thStartTimes = find(temp == 1);
thEndTimes   = find(temp == -1);


if thStartTimes(1) > thEndTimes(1)
    thStartTimes = [1; thStartTimes];
end


if thStartTimes(end) > thEndTimes(end)
    thEndTimes = [thEndTimes; numel(insideTheta)];
end


thetaPeriods = [t(thStartTimes) t(thEndTimes)];
% thetaDur     = thetaPeriods(:, 2)-thetaPeriods(:, 1);


% thetaPeriods(thetaDur < minLen, :) = []; % detected theta periods with duration less then 2 seconds are excluded.

end



function [periods, upIdx] = schnmidtTrigger(variable, oldt, t)



variable = interp1(oldt, variable, t);


rmvIdx = variable > prctile(variable, 97.5) | variable < prctile(variable, 2.5);
variable2 = variable;
variable2(rmvIdx) = [];



numpeaks = 1;
numbins = 10; 
while numpeaks ~=2
    [bihist.hist,bihist.bins]= hist(variable2,numbins);

    [PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

% % numbins = 30;
% [bihist.hist,bihist.bins]= hist(variable, numbins);
% [PKS,LOCS] = findpeaks([0 bihist.hist 0]);
% LOCS = sort(LOCS)-1;


betweenpeaks = bihist.bins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-bihist.hist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

thresh = betweenpeaks(diploc);


threshUP = thresh + 0.2.*(betweenpeaks(end)-thresh);
threshDOWN = thresh - 0.2.*(thresh- betweenpeaks(1));
    


overUP   = variable > threshUP;
overDOWN = variable > threshDOWN;

crossup = find(diff(overUP)==1);
crossdown = find(diff(overDOWN)==-1);



%Delete incomplete (repeat) crossings
allcrossings = [crossup ones(size(crossup)) ;...
    crossdown zeros(size(crossdown))];
[~,sortorder] = sort(allcrossings(:,1));
allcrossings = allcrossings(sortorder,:);
updownswitch = diff(allcrossings(:,2));
samestate = find(updownswitch==0)+1;
allcrossings(samestate,:)=[];

crossup = allcrossings(allcrossings(:,2)==1,1);
crossdown = allcrossings(allcrossings(:,2)==0,1);



if crossdown(1) < crossup(1) 
    crossdown(1) = [];
end
if crossup(end) > crossdown(end) 
    crossup(end) = [];
end


periods = t([crossup crossdown]);




% convert to indices for successive time bins (divided based on the desired time resolution)

upIdx = zeros(numel(t), 1);

for ip = 1:size(periods, 1)
    upIdx(t >= periods(ip, 1) & t <= periods(ip, 2)) = 1;
end

end


function  intervals = IdxtoInt(Idx, timestamps)

    statetimes = Idx==1;
    %Get indices of state on/offsets

    crossup = find([0; diff(statetimes)] == 1);
    crossdown = find([0; diff(statetimes)] == -1);

    if crossdown(1) < crossup(1) 
        crossdown(1) = [];
    end
    if crossup(end) > crossdown(end) 
        crossup(end) = [];
    end

    stateints = [crossup crossdown];
    intervals = [timestamps(stateints(:,1)) timestamps(stateints(:,2))];

end


function [allPeriods, durations] = concatSameState(allPeriods)

totalCounts = numel(find(allPeriods(1:end-1, 2) == allPeriods(2:end, 1) & allPeriods(1:end-1, 3) == allPeriods(2:end, 3)));

for ii = 1:totalCounts

    currTransition = find(allPeriods(1:end-1, 2) == allPeriods(2:end, 1) & allPeriods(1:end-1, 3) == allPeriods(2:end, 3), 1, 'first');
    allPeriods(currTransition, 2) = allPeriods(currTransition+1, 2);
    allPeriods(currTransition+1, :) = [];
end

durations = diff(allPeriods(:, 1:2), [], 2);


end





