function [thetaPeriods, thStatesMeans] = ThetaDetect_sep(basePath)


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
theta_f_in.pre = [5 10];
theta_f_out.pre = [1 4 10 12];

theta_f_in.run = [5 10];
theta_f_out.run = [1 5 10 14];

theta_f_in.post = [5 10];
theta_f_out.post = [1 4 10 12];



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

for iepoch = 1:3
    
    currEpoch = epochNames{iepoch};
    period    = epochs.(currEpoch); % period during which the theta detection is performed
            
    in_period_idx.(currEpoch) = find(t >= period(1) & t < period(2)); 
    
    thfin = find(f > theta_f_in.(currEpoch)(1) & f < theta_f_in.(currEpoch)(2));
    
    if numel(theta_f_out.(currEpoch)) == 4
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2))|(f > theta_f_out.(currEpoch)(3) & f < theta_f_out.(currEpoch)(4))); 
    elseif numel(theta_f_out.(currEpoch)) == 2
        thfout = find((f > theta_f_out.(currEpoch)(1) & f < theta_f_out.(currEpoch)(2)));
    end

    [thratio.(currEpoch), inTh.(currEpoch), thetaPeriods.(currEpoch)] = TDRatioAuto(Pxx, nStates, t, in_period_idx.(currEpoch), thfin, thfout, Par);
    
end


%%

in_period_idx_all = [in_period_idx.pre; in_period_idx.run; in_period_idx.post];
theratio_all      = [thratio.pre; thratio.run; thratio.post];
inTh_all          = [inTh.pre inTh.run inTh.post];

thetaPeriods_all = [thetaPeriods.pre; thetaPeriods.run; thetaPeriods.post];


t = t(in_period_idx_all);
Pxx = Pxx(in_period_idx_all, :);



period2plot = epochs.pre;

tindex = find(t>period2plot(1) & t<period2plot(2));
% plot the spectrogram and theta boundaries

figure; 

ax1 = subplot(2,1,1);

spectrogram2Plot = log10(Pxx(tindex, 1:find(f>100, 1, 'first'), 1))'; %

cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
imagesc(t(tindex), f(1:find(f > 100, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 




ax2 = subplot(2,1,2);

ts = 1/lfpSampleRate;
teeg = (period2plot(1)+ts):ts:period2plot(2);
% plot((1:length(Eeg))/lfpSampleRate, Eeg(:, 1)/range(Eeg(:, 1))*30, 'linewidth', 1, 'color', [1 0 0 0.5])

currentEeg = Eeg((period2plot(1)*lfpSampleRate+1):period2plot(2)*lfpSampleRate, 1);
plot(teeg(1:length(currentEeg)), currentEeg/range(currentEeg)*30, 'linewidth', 1, 'color', [1 0 0 0.5])

hold on

plot(t(tindex), theratio_all(tindex), 'linewidth', 3, 'color', [0 0 1])
hold on
stairs(t(tindex), inTh_all(tindex)*3, 'linewidth', 3, 'color', [0 0 0 1])

ylim([-6 6])


linkaxes([ax1 ax2], 'x')

% set(ax2, 'box', 'off', 'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
legend('raw LFP','theta/non-theta ratio', 'theta: 5-11Hz, deltaplus:1-4 Hz & 12-15 Hz, nStates=2'); % , 'theta: 4-8Hz, non-theta:1-3H & 9-11Hz, nStates=4', 'either' 




%% choose which method is most accurate

SaveFileName = fullfile(outputDir, [sessionName '.theta.1']);
                        
dlgoptions.Resize = 'off';
saveans = inputdlg({'Enter filename '},'Save Periods', 1,{SaveFileName},dlgoptions);

msave(saveans{1}, thetaPeriods_all);

MakeEvtFile(thetaPeriods_all,fullfile(outputDir, [sessionName '.the.evt']), {'beg', 'end'}, lfpSampleRate); % make evt file for neuroscope browsing


inTh = inTh_all';
timePnts = t;
thetaPeriods = thetaPeriods_all/lfpSampleRate;
theratio = theratio_all;


detectionParams.theta_f_in  = theta_f_in;
detectionParams.theta_f_out = theta_f_out;
detectionParams.winLen      = Window;
detectionParams.WinOverlapRatio = overlapRatio;
detectionParams.thetaChannel = thetaChannel;
detectionParams.date = date;


save(fullfile(outputDir, [sessionName '.thetaStates.mat']), 'inTh', 'timePnts', 'theratio', 'thetaPeriods', 'detectionParams')


end


function [thratio, InTh, thetaPeriods, thratio_st] = TDRatioAuto(Pxx, nStates, t, in_period_idx, thfin, thfout, Par)

%automatic theta periods detection just using the thetaratio

% thfin = find(f>5 & f<11);
% thfout = find((f>1 & f<4)|(f> 12& f<15));
% thfout = find((f>1 & f<4));


thratio = log(mean(squeeze(Pxx(in_period_idx,thfin,1)),2))-log(mean(squeeze(Pxx(in_period_idx,thfout,1)),2));

    
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
