function [thetaPeriods, thStatesMeans] = ThetaDetect_EMGadded_forREMDetection2(basePath)


% theta detection separately performed on each epoch. The reason was that
% the theta frequency band and the frequency gap between thata and its
% first harmonic is different between the sleep and maze epochs. We need to
% adjust the frequency bands separately for each epoch.

% DON'T TRUST THE THETA PERIODS CALCUALTED USING HMM
% In epochs dominated with NREM, two states might not be fully enough to capture the variablity in theta/delta ratio, so the state detection might lack enough accuracy. 
% Either exclude the NREM before HMM or increase the number of states (how
% many would be right?). 

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
epochNames = {'pre'; 'run'; 'post'};



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


% for Grosmark's sessions
% 
% load(fullfile(basePath, [sessionName '_sessInfo.mat']), 'sessInfo') % for Grosmark sessions
% 
% pre  = sessInfo.Epochs.PREEpoch;
% maze = sessInfo.Epochs.MazeEpoch;
% post = sessInfo.Epochs.POSTEpoch;
% 
% epochs.pre = pre;
% epochs.run = maze;
% epochs.post = post;



% % theta and outside-theta frequency bands

% Rat S
theta_f_in.pre = [5 10];
theta_f_out.pre = [1 4 10 12];

theta_f_in.run = [5 10];
theta_f_out.run = [1 5 11 13];

theta_f_in.post = [5 10];
theta_f_out.post = [1 4 10 12];

fileName = fullfile('~/Documents/NCMLproject/assemblyTuning_finalResults', sessionName, 'spikes', [sessionName '.brainState.mat']);
load(fileName, 'brainState')

emg      = brainState.emg;
slowWave = brainState.slowWave;

emgThresh = brainState.emgThresh;
swThresh  = brainState.swthresh;



% % loading the eeg file (I should add a procedure to calculate the best theta channel)

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



%% apply HMM on theta ratio

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




%% 



% for iepoch = 1:numel(epochNames)
%     currEpoch = epochNames{iepoch};
%     thratio.(currEpoch) = conv(thratio.(currEpoch), gw, 'same');
% end

in_period_idx_all = [in_period_idx.pre; in_period_idx.run; in_period_idx.post];
t = t(in_period_idx_all);

theratio_all      = [thratio.pre; thratio.run; thratio.post];
thetaPeriods_all = [thetaPeriods.pre; thetaPeriods.run; thetaPeriods.post];
inTh_all          = [inTh.pre inTh.run inTh.post];


% %% hard threshold
%
% thetaThresh = 1;
%
% temp = theratio_all;
% idx = temp > thetaThresh;
% 
% crossup   = find([0; diff(idx)] == 1);
% crossdown = find([0; diff(idx)] == -1);
% 
% if crossdown(1) < crossup(1)
%     crossdown(1) = [];
% end
% 
% if crossup(end) > crossdown(end)
%     crossup(end) = [];
% end
% 
% 
% thetaPeriods_all = t([crossup crossdown]); 
% thetaPeriods_all(diff(thetaPeriods_all')' < 3, :) = []; 



%% detecting brain states by schmidt triggering theta, emg, slow wave


%% theta 

% nBins = 60;
% [thetaPeriods, thetaUpIdx] = schnmidtTrigger(theratio_all, thetaThresh, t, t, nBins);
% thetaPeriods(diff(thetaPeriods')' < 3, :) = [];


%% emg
nBins = 30;
[emgPeriods, emgUpIdx] = schnmidtTrigger(emg, emgThresh, brainState.sw_emg_timePnts, t, nBins);



%% slow wave
nBins = 30;
[swPeriods, swUpIdx] = schnmidtTrigger(slowWave, swThresh, brainState.sw_emg_timePnts, t, nBins);






%% plots

t = t(in_period_idx_all);
Pxx = Pxx(in_period_idx_all, :);


% period2plot = epochs.pre;
period2plot = [epochs.post(1) epochs.post(1)+5*60*60];


tindex = find(t>period2plot(1) & t<period2plot(2));
% plot the spectrogram and theta boundaries


plotheight = 450;
plotwidth  = 900;
fontsize   = 8;

figure;
clf(gcf);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);




% spectrogram
ax(1) = subplot(13, 1, 1:4);

spectrogram2Plot = log10(Pxx(tindex, 1:find(f>100, 1, 'first'), 1))'; %

cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
imagesc(t(tindex), f(1:find(f > 100, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 

caxis([0.7 6])

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)




% raw EEG 
ax(2) = subplot(13, 1, 5);

ts = 1/lfpSampleRate;
teeg = (period2plot(1)+ts):ts:period2plot(2);

currentEeg = Eeg((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg), 'linewidth', 1, 'color', [0 0 0 0.5]);

ylabel('raw LFP', 'FontSize', fontsize)

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)



% theta LFP

ax(3) = subplot(13, 1, 6:7);

fNum = 2^floor(log(600*lfpSampleRate)/log(2)+1);

forder = 256;  
forder = ceil(forder/2)*2; % to make sure filter order is even

fTheta = fir1(forder,theta_f_in.pre/lfpSampleRate*2); % calculate convolution func
thetaLFP = Filter0(fTheta, Eeg);

% dTheta = filtfilt(fTheta, 1, Eeg);

currentEeg = thetaLFP((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg), 'linewidth', 1, 'color', [0 0 0 0.5]);

ylabel(sprintf('theta (%d - %d)', theta_f_in.pre(1), theta_f_in.pre(2)), 'fontsize', fontsize)

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)



% theta ratio

ax(4) = subplot(13, 1, 8:9);

hold on

h1 = plot(t(tindex), theratio_all(tindex), 'linewidth', 1, 'color', [0 0 1], 'DisplayName', 'theta ratio');
h2 = stairs(t(tindex), inTh_all(tindex)*2, 'linewidth', 1, 'color', [0 0 0 1] , 'DisplayName', 'HMM states');
% currThetaPeridos = thetaPeriods_all(thetaPeriods_all(:, 1) > period2plot(1) & thetaPeriods_all(:, 1) < period2plot(2), :);
currThetaPeridos = thetaPeriods(thetaPeriods(:, 1) >= period2plot(1) & thetaPeriods(:, 1) <= period2plot(2), :);

for ii = 1:size(currThetaPeridos, 1)
    h3 = line([currThetaPeridos(ii, 1) currThetaPeridos(ii, 2)], [2.5 2.5], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'theta periods(hard threshold)');
end

h4 = plot(t(tindex), thetaThresh*ones(numel(find(tindex)), 1), 'linewidth', 1, 'color', [0 0 1], 'linestyle', ':', 'DisplayName', 'theta threshold');

ylim([-1 3])

% legend([h2, h3, h4])
legend([h2, h4])

ylabel('theta ratio', 'FontSize', fontsize, 'Color', [0 0 1])

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)



% EMG

tindex2 = find(brainState.sw_emg_timePnts>period2plot(1) & brainState.sw_emg_timePnts<period2plot(2));
ax(5) = subplot(13, 1, 10:11);

hold on
plot(brainState.sw_emg_timePnts(tindex2), brainState.emg(tindex2), 'linewidth', 1, 'color', [0 0 0], 'DisplayName', 'EMG')

h1 = plot(t(tindex), emgThresh*ones(numel(find(tindex)), 1), 'linewidth', 1, 'color', [0 0 0], 'linestyle', ':', 'DisplayName', 'EMG threshold');

currEMGPeridos = emgPeriods(emgPeriods(:, 1) >= period2plot(1) & emgPeriods(:, 1) <= period2plot(2), :);

for ii = 1:size(currEMGPeridos, 1)
    h2 = line([currEMGPeridos(ii, 1) currEMGPeridos(ii, 2)], [0.8 0.8], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'hi emg periods');
end

legend([h1 h2])


ylabel('EMG', 'FontSize', fontsize)

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)


% SWS
ax(6) = subplot(13, 1, 12:13);

hold on
plot(brainState.sw_emg_timePnts(tindex2), brainState.slowWave(tindex2), 'linewidth', 1, 'color', [0 0 0], 'DisplayName', 'EMG')

h1 = plot(t(tindex), brainState.swthresh*ones(numel(find(tindex)), 1), 'linewidth', 1, 'color', [0 0 0], 'linestyle', ':', 'DisplayName', 'SW threshold');

currSWPeridos = swPeriods(swPeriods(:, 1) >= period2plot(1) & swPeriods(:, 1) <= period2plot(2), :);

for ii = 1:size(currSWPeridos, 1)
    h2 = line([currSWPeridos(ii, 1) currSWPeridos(ii, 2)], [0.8 0.8], 'linewidth', 1, 'color', [1 0 0 1], 'DisplayName', 'hi sws periods');
end

legend([h1 h2])


ylabel('SWS', 'FontSize', fontsize)

set(gca,'box', 'off', 'TickDir', 'out','TickLength',[0.01, 0.02], 'fontsize', 6)


linkaxes(ax, 'x')
xlim(period2plot)


% set(ax2, 'box', 'off', 'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
% legend(['position', 'theta/non-theta ratio', 'theta periods']); % , 'theta: 4-8Hz, non-theta:1-3H & 9-11Hz, nStates=4', 'either' 

%, 'theta: 5-10Hz, deltaplus:1-4 Hz & 10-14 Hz, nStates=2''raw LFP',


%% choose which method is most accurate

SaveFileName = fullfile(outputDir, [sessionName '.theta.1']);
                        
dlgoptions.Resize = 'off';
saveans = inputdlg({'Enter filename '},'Save Periods', 1,{SaveFileName},dlgoptions);

msave(saveans{1}, thetaPeriods_all);

MakeEvtFile(thetaPeriods_all,fullfile(outputDir, [sessionName '.the.evt']), {'beg', 'end'}, lfpSampleRate); % make evt file for neuroscope browsing


% inTh = inTh_all';
timePnts = t;
theratio = theratio_all;


detectionParams.theta_f_in      = theta_f_in;
detectionParams.theta_f_out     = theta_f_out;
detectionParams.winLen          = Window;
detectionParams.WinOverlapRatio = overlapRatio;
detectionParams.thetaChannel    = thetaChannel;
detectionParams.date            = date;


save(fullfile(outputDir, [sessionName '.thetaStates_including_remaze_stricterDetection.mat']), 'timePnts', 'theratio', 'thetaPeriods_all', 'detectionParams') % 'inTh'


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



function [periods, upIdx] = schnmidtTrigger(variable, thresh, oldt, t, numbins)



variable = interp1(oldt, variable, t);

% numpeaks = 1;
% numbins = 10; 
% while numpeaks ~=2
%     [bihist.hist,bihist.bins]= hist(variable,numbins);
% 
%     [PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
%     LOCS = sort(LOCS)-1;
%     numbins = numbins+1;
%     numpeaks = length(LOCS);
% end

% numbins = 30;
[bihist.hist,bihist.bins]= hist(variable, numbins);
[PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
LOCS = sort(LOCS)-1;


betweenpeaks = bihist.bins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-bihist.hist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

thresh = betweenpeaks(diploc);


threshUP = thresh + 0.5.*(betweenpeaks(end)-thresh);
threshDOWN = thresh + 0.5.*(betweenpeaks(1)-thresh);
    


overUP = variable>threshUP;
overDOWN = variable>threshDOWN;

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



if crossup(1) < crossdown(1)
    crossup(1) = [];
end
if crossdown(end) > crossup(end)
    crossdown(end) = [];
end
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

