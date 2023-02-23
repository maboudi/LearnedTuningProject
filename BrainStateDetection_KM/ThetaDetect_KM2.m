function [thetaPeriods, thStatesMeans] = ThetaDetect_KM2(fileBase, thetaChannels)



Par = LoadXml([fileBase '.xml']);
lfpSampleRate = Par.lfpSampleRate;
nCh = Par.nChannels;

% epochs

load(fullfile(fileBase, 'mat_files', [fileBase '.epochs.mat']), 'pre', 'post', 'maze')
pre = double(pre);
maze = double(maze);
post = double(post);


% minLen = 2; %seconds
Window = 2; %sec - for spectrogram calculation
FreqRange = [1 100]; % Hz - range of freq. for the spectrogram



% % calculate the spectrogram

Eeg = readmulti([fileBase '.eeg'], nCh, thetaChannels);
preEpoch = pre*lfpSampleRate;

Eeg = Eeg(preEpoch(1)+1: preEpoch(2));

SpecWindow = 2^round(log2(Window*lfpSampleRate));% choose window length as power of two
nFFT = SpecWindow*4;
weeg = WhitenSignal(Eeg,lfpSampleRate*2000,1);
[Pxx,f,t]=mtcsglong(weeg,nFFT,lfpSampleRate,SpecWindow,[],2,'linear',[],FreqRange);



% [smoothdata,filtwts] = eegfilt(Eeg',1250,59,61,0,0,1, 'fir1');

%% calculate theta2surround power ratio and run HMM to detect the theta periods

nStates = 2;
thfin = find(f>5 & f<10);
thfout = find((f>1 & f<4)|(f> 11 & f<13)); % |(f> 11& f<13)

[thratio, inTh, thetaPeriods] = TDRatioAuto(Pxx, nStates, t, thfin, thfout, Par);
% RUNthetaPeriods = calThetaPeriods(RUNinsideTheta, t, minLen);

thetaPeriods = thetaPeriods + preEpoch(1);


% plot the spectrogram and theta boundaries using each
% method

figure; 

ax1 = subplot(2,1,1);

spectrogram2Plot = log10(Pxx(:,1:find(f>100, 1, 'first'), 1))'; %

cl = [min(spectrogram2Plot(:)) min(spectrogram2Plot(:))+1*(range(spectrogram2Plot(:)))];
imagesc(t,f(1:find(f>100, 1, 'first')), spectrogram2Plot, cl); colormap('jet'); set(gca, 'YDir', 'normal'); 
% 

ax2 = subplot(2,1,2);


plot((1:length(Eeg))/lfpSampleRate, Eeg(:, 1)/range(Eeg(:, 1))*30, 'linewidth', 1, 'color', [1 0 0 0.5])
hold on

plot(t, thratio*4, 'linewidth', 3, 'color', [0 0 1])
hold on
stairs(t, inTh*12, 'linewidth', 3, 'color', [0 0 0 1])
% hold on
% plot(t, (SLEEPinsideTheta-1.1)*12, 'linewidth', 3, 'color', [0 0 0 0.5])
% hold on
% plot(t, (double(RUNinsideTheta | SLEEPinsideTheta) + 1.1)*12, 'linewidth', 3, 'color', 'g')
ylim([-20 20])


linkaxes([ax1 ax2], 'x')
set(ax2, 'box', 'off', 'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

legend('raw LFP','theta/non-theta ratio', 'theta: 5-11Hz, deltaplus:1-4 Hz & 12-15 Hz, nStates=2'); % , 'theta: 4-8Hz, non-theta:1-3H & 9-11Hz, nStates=4', 'either' 
% xlim(sleepPeriod)

% mkdir(fullfile(fileBase, 'thetaDetection'))
% saveas(gcf, fullfile(fileBase, 'thetaDetection', 'samplePeriod2_nStates=2.png'))



%% choose which method is most accurate


save('Maze_thetaStates.mat', 't', 'inTh', 'thetaPeriods')

SaveFileName = [fileBase '.MAZE.theta.1' ];
                        
dlgoptions.Resize = 'off';
saveans = inputdlg({'Enter filename '},'Save Periods', 1,{SaveFileName},dlgoptions);

msave(saveans{1}, thetaPeriods);


MakeEvtFile(thetaPeriods,[fileBase '.MAZE.the.evt'], {'beg', 'end'}, lfpSampleRate); % make evt file for neuroscope browsing





% Add theta-delta ratio to the DAT file

% eegInfo = dir([fileBase '.eeg']);
% noTimePnts = eegInfo.bytes/nCh/2;
% 
% samplePnts = (1: noTimePnts)* 1/lfpSampleRate;
% 
% sdtAdd2eeg = interp1((1:numel(insideTheta))*deltT, thratio, samplePnts);
% % sdtAdd2eeg = thratio;
% 
% sdtAdd2eeg = (2 * (sdtAdd2eeg - min(sdtAdd2eeg)) / (max(sdtAdd2eeg) - min(sdtAdd2eeg)) - 1)*12000;
% 
% 
% inputFileHandle = fopen([fileBase,'.eeg']);
% outputFileHandle=fopen([fileBase,'-2.eeg'],'w');
% 
% 
% bufferSize=4096;
% 
% doneFrame=0;
% while ~feof(inputFileHandle)
%     
%     data = fread(inputFileHandle,[nCh,bufferSize],'int16')';
%     for frame=1:size(data,1)
%         fwrite(outputFileHandle,[data(frame,:), sdtAdd2eeg(frame+doneFrame)]','int16');
%     end
%     doneFrame=doneFrame+frame;
%     
% end
% 
% 
% fclose(inputFileHandle);
% fclose(outputFileHandle);


end


function [thratio, InTh, thetaPeriods, thratio_st] = TDRatioAuto(Pxx, nStates, t, thfin, thfout, Par)

%automatic theta periods detection just using the thetaratio

% thfin = find(f>5 & f<11);
% thfout = find((f>1 & f<4)|(f> 12& f<15));
% thfout = find((f>1 & f<4));


thratio = log(mean(squeeze(Pxx(:,thfin,1)),2))-log(mean(squeeze(Pxx(:,thfout,1)),2));

    
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


% for i = 1:nStates
%     InTh(TheState == sortInd(i)) = i-1;
% end


thebeg= SchmittTrigger(InTh,0.9, 0.9);
theend= SchmittTrigger(-InTh,-0.9, -0.9);

if thebeg(1)>theend(1)
    theend =theend(2:end);
end

if thebeg(end)>theend(end)
    thebeg =thebeg(1:end-1);
end

thetaPeriods =round([t(thebeg) t(theend)]*Par.lfpSampleRate);

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
thetaDur     = thetaPeriods(:, 2)-thetaPeriods(:, 1);


thetaPeriods(thetaDur < minLen, :) = []; % detected theta periods with duration less then 2 seconds are excluded.

end

% intervals = [nan ;t(thStartTimes(2:end)) - t(thEndTimes(1:end-1))];
% periodsWithDistless = find(intervals < 1);
% 
% thStartTimes2 = thStartTimes;
% thStartTimes2(periodsWithDistless) = 0;
% 
% thEndTimes2 = thEndTimes;
% thEndTimes2(periodsWithDistless - 1) = 0;
% 
% thStartTimes2(thStartTimes2 == 0) = [];
% thEndTimes2(thEndTimes2 == 0) = [];
