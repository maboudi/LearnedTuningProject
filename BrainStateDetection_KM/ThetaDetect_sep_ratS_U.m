function [thetaPeriods, thStatesMeans] = ThetaDetect_sep_ratS_U(basePath)


% theta detection separately performed on each epoch. The reason was that
% the theta frequency band and the frequency gap between thata and its
% first harmonic is different between the sleep and maze epochs. We need to
% adjust the frequency bands separately for each epoch.

% DON'T TRUST THE THETA PERIODS CALCUALTED USING HMM
% In epochs dominated with NREM, two states might not be fully enough to capture the variablity in theta/delta ratio, so the state detection might lack enough accuracy. 
% Either exclude the NREM before HMM or increase the number of states (how
% many would be enough?). 



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
epochNames = {'pre'; 'run'; 'post'; 'run2'};


% for Bapun sessions add the follwoing steps
 
load(fullfile(basePath, [sessionName '.epochs.mat']))

pre  = double(pre);
post = double(post);

if exist('maze', 'var')
    maze = double(maze);
elseif exist('maze1', 'var')
    maze = double(maze1);
end

if exist('maze2')
    remaze = double(maze2);
elseif exist('re_maze')
    remaze = double(re_maze);
end


epochs.pre = pre;

epochs.run(1) = pre(2);
epochs.run(2) = maze(2);

epochs.post(1) = maze(2);
epochs.post(2) = post(2);

epochs.run2(1) = post(2);
epochs.run2(2) = remaze(2); 


% position data
fileName = fullfile('~/Documents/NCMLproject/assemblyTuning_finalResults', sessionName, 'spikes', [sessionName '.spikes.mat']);
load(fileName, 'fileInfo')


linearPos = fileInfo.linearPos;


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

theta_f_in.run2 = [5 10];
theta_f_out.run2 = [1 5 11 13];


% % Rat U
% theta_f_in.pre = [5 10];
% theta_f_out.pre = [1 4 10 12];
% 
% theta_f_in.run = [6 10];
% theta_f_out.run = [1 5 11 14];
% 
% theta_f_in.post = [5 10];
% theta_f_out.post = [1 4 10 12];
% 
% theta_f_in.run2 = [6 10];
% theta_f_out.run2 = [1 5 11 14];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = fullfile('~/Documents/NCMLproject/assemblyTuning_finalResults', sessionName, 'spikes', [sessionName '.brainState.mat']);
load(fileName, 'brainState')

load(fullfile('~/Documents/NCMLproject/assemblyTuning_finalResults', sessionName 'PopulationBurstEvents', [sessionName '.PBEInfo.mat']), 'sdat')


sdat_t = (1:length(sdat))*1e-3;

minREMDur    = 0.5;
min_interREM = 0;


gw = gausswindow(3,9);
brainState.theratio = conv(brainState.theratio, gw, 'same');


theta_intrp = interp1(brainState.thetaTimePnts, brainState.theratio, sdat_t);
sw_intrp    = interp1(brainState.sw_emg_timePnts, brainState.slowWave, sdat_t);
emg_intrp   = interp1(brainState.sw_emg_timePnts, brainState.emg, sdat_t);


% REM periods
REMIdx = sw_intrp' < brainState.swthresh & theta_intrp' > 1 & emg_intrp' < brainState.emgThresh; % brainState.swthresh 

% REM bouts
try 
    crossup   = find([0; diff(REMIdx)] == 1);
    crossdown = find([0; diff(REMIdx)] == -1);
catch 
    crossup   = find([0 diff(REMIdx)] == 1);
    crossdown = find([0 diff(REMIdx)] == -1);
end


if crossdown(1) < crossup(1); crossdown(1) = []; end
if crossup(end) > crossdown(end); crossup(end) = []; end

REMBouts = sdat_t([crossup crossdown]);



% concatenate REM bouts that are less than min_interREM apart
% We don't do it for NREM periods as the concatenation might lead to
% inclusion of population burst events .. 
    
tmp_rem = REMBouts(1,:);
finalREMBouts = [];

for ii = 2:size(REMBouts,1)
    if (REMBouts(ii,1)-tmp_rem(2)) < min_interREM
        tmp_rem = [tmp_rem(1) REMBouts(ii, 2)]; % merge two PBEs
    else
        finalREMBouts = [finalREMBouts; tmp_rem];
        tmp_rem = REMBouts(ii,:);
    end
end

finalREMBouts = [finalREMBouts;tmp_rem]; 

REMBouts = finalREMBouts;

REMdur  = REMBouts(:, 2) - REMBouts(:, 1);
% REMdist = REMBouts(2:end, 1) - REMBouts(1:end-1, 2); 

REMBouts(REMdur < minREMDur, :) = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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




%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 
% Window = 1; %sec - for spectrogram calculation
% FreqRange = [1 300]; % Hz - range of freq. for the spectrogram
% SpecWindow = 2^round(log2(Window*lfpSampleRate));% choose window length as power of two
% overlapRatio = 0.75;
% 
% 
% nFFT = SpecWindow*4;
% weeg = WhitenSignal(Eeg,lfpSampleRate*2000,1);
% 
% [Pxx,f,t] = mtcsglong(weeg,nFFT, lfpSampleRate, SpecWindow, overlapRatio*SpecWindow, 2,'linear',[],FreqRange);
% 
% % t starts from zero; the start time should be adjusted to the center of
% % the first window
% t = t + SpecWindow/2/lfpSampleRate;
% 
% %%
% 
% for iepoch = 1:3
%     currEpoch = epochNames{iepoch};
%     period    = epochs.(currEpoch); % period during which the theta detection is performed
%             
%     in_period_idx.(currEpoch) = find(t >= period(1) & t < period(2)); 
% end
% 
% 
% in_period_idx_all = [in_period_idx.pre; in_period_idx.run; in_period_idx.post];
% 
% t = t(in_period_idx_all);
% Pxx = Pxx(in_period_idx_all, :);
% 
% 
% 
% 
% period2plot = epochs.run; 
% 
% tindex = find(t>period2plot(1) & t<period2plot(2));
% % plot the spectrogram and theta boundaries
% 
% figure; 
% 
% ax1 = subplot(2,1,1);
% 
% spectrogram2Plot = log10(Pxx(tindex, 1:find(f>300, 1, 'first'), 1))'; %
% 
% cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
% imagesc(t(tindex), f(1:find(f > 300, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

nStates = 2; % parameters for HMM

for iepoch = 1:numel(epochNames)
    
    currEpoch = epochNames{iepoch};
    period    = epochs.(currEpoch); % period during which the theta detection is performed
            
    in_period_idx.(currEpoch) = find(t >= period(1) & t < period(2)); 
    

    if numel(theta_f_in.(currEpoch)) == 2
        thfin = find(f > theta_f_in.(currEpoch)(1) & f < theta_f_in.(currEpoch)(2));    
    elseif numel(theta_f_in.(currEpoch)) == 4
        thfin = find((f > theta_f_in.(currEpoch)(1) & f < theta_f_in.(currEpoch)(2)) | (f > theta_f_in.(currEpoch)(3) & f < theta_f_in.(currEpoch)(4)) );
    end


    if numel(theta_f_out.(currEpoch)) == 4
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2))|(f > theta_f_out.(currEpoch)(3) & f < theta_f_out.(currEpoch)(4))); 
    elseif numel(theta_f_out.(currEpoch)) == 2
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2)));
    end

    [thratio.(currEpoch), inTh.(currEpoch), thetaPeriods.(currEpoch)] = TDRatioAuto(Pxx, nStates, t, in_period_idx.(currEpoch), thfin, thfout, Par);
    
end



%%


gw = gausswindow(3,9);
thetaThresh = 1;

for iepoch = 1:numel(epochNames)
    
    currEpoch = epochNames{iepoch};
    thratio.(currEpoch) = conv(thratio.(currEpoch), gw, 'same');

end


in_period_idx_all = [in_period_idx.pre; in_period_idx.run; in_period_idx.post; in_period_idx.run2];
theratio_all      = [thratio.pre; thratio.run; thratio.post; thratio.run2];
inTh_all          = [inTh.pre inTh.run inTh.post inTh.run2];

temp = theratio_all;
idx = temp > thetaThresh;
  
crossup   = find([0; diff(idx)] == 1);
crossdown = find([0; diff(idx)] == -1);

if crossdown(1) < crossup(1)
    crossdown(1) = [];
end

if crossup(end) > crossdown(end)
    crossup(end) = [];
end

thetaPeriods_all = t([crossup crossdown]); 
thetaPeriods_all(diff(thetaPeriods_all')' < 3, :) = []; % minimum duration = 3


% thetaPeriods_all = [thetaPeriods.pre; thetaPeriods.run; thetaPeriods.post; thetaPeriods.run2];

t = t(in_period_idx_all);
Pxx = Pxx(in_period_idx_all, :);



period2plot = epochs.pre;


tindex = find(t>period2plot(1) & t<period2plot(2));
% plot the spectrogram and theta boundaries


plotheight = 450;
plotwidth  = 900;

figure;
clf(gcf);
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



ax1 = subplot(2,1,1);

spectrogram2Plot = log10(Pxx(tindex, 1:find(f>100, 1, 'first'), 1))'; %

cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
imagesc(t(tindex), f(1:find(f > 100, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 

% ylim(yl);
% xlim(xl);

caxis([0.7 6])


ax2 = subplot(2,1,2);

ts = 1/lfpSampleRate;
teeg = (period2plot(1)+ts):ts:period2plot(2);
% plot((1:length(Eeg))/lfpSampleRate, Eeg(:, 1)/range(Eeg(:, 1))*30, 'linewidth', 1, 'color', [1 0 0 0.5])

currentEeg = Eeg((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
h1 = plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg)*50-3, 'linewidth', 1, 'color', [1 0 0 0.5]);

hold on



% idx = find(linearPos(:,2) > period2plot(1) & linearPos(:,2) < period2plot(2));
% h1 = plot(linearPos(idx, 2), linearPos(idx, 1)/range(linearPos(:, 1))*10, 'linewidth', 1, 'color', [0 1 0 0.5], 'DisplayName', 'position');




hold on

h2 = plot(t(tindex), theratio_all(tindex), 'linewidth', 1, 'color', [0 0 1], 'DisplayName', 'theta ratio');
hold on
stairs(t(tindex), inTh_all(tindex)*3, 'linewidth', 1, 'color', [0 0 0 1])

currThetaPeridos = thetaPeriods_all(thetaPeriods_all(:, 1) > period2plot(1) & thetaPeriods_all(:, 1) < period2plot(2), :);

remaze_theta  = thetaPeriods_all(thetaPeriods_all(:, 1) > maze(1) & thetaPeriods_all(:, 1) < maze(2), :);


for ii = 1:size(currThetaPeridos, 1)
    h4 = line([currThetaPeridos(ii, 1) currThetaPeridos(ii, 2)], [2 2], 'linewidth', 1, 'color', [0 1 0 1], 'DisplayName', 'Thresholded theta');
end


hold on
h3 = plot(t(tindex), thetaThresh*ones(numel(find(tindex)), 1), 'linewidth', 1, 'color', [0 0 1], 'linestyle', ':', 'DisplayName', 'theta threshold');


ylim([-6 6])


linkaxes([ax1 ax2], 'x')

legend([h2, h3, h4])

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

gw = gausswindow(3,9);
thratio = conv(thratio, gw, 'same');

    
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
