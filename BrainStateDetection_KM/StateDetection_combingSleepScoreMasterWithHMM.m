load('SleepScoreMetrics.mat')
broadbandSlowWave   = SleepScoreMetrics.broadbandSlowWave;
t_clus   = SleepScoreMetrics.t_clus;
swthresh = SleepScoreMetrics.histsandthreshs.swthresh;

EMG = SleepScoreMetrics.EMG;
thratio = SleepScoreMetrics.thratio;


% load theta state indices calculated using checkEegState
load('ThetaStates_CheckEegStates.mat', 'inTh', 't')
isTheta = inTh;
isTheta_t = t;
isTheta_ip = round(interp1(isTheta_t, isTheta, t_clus)); % interpolated

if size(isTheta_ip, 2) > size(isTheta_ip, 1)
    isTheta_ip = isTheta_ip';
end

figure; scatter(thratio, EMG, 5, 'filled')

% figure; hist(bbSlowWave(isTheta_ip == 1), 20)

EMGthresh = 0.2;


% % the following pieces of code were copied from Buzcode state detection

%% Re-Do this code (should be same as in ClusterStates_GetParams.m) to see if theta is bimodal

%This switch turns on a "schmidt trigger", or sticky trigger,
%which means that threshold crossings have to reach the
%midpoint between the dip and the opposite peak, this
%reduces noise. Passed through via histsandthreshs from checkboxes in
%TheStateEditor or 'stickytrigger',true in SleepScoreMaster via GetMetrics
if ~exist('stickySW','var'); stickySW = false; end
if ~exist('stickyTH','var'); stickyTH = false; end
if ~exist('stickyEMG','var'); stickyEMG = false; end


[~,~,~,~,NREMtimes] = bz_BimodalThresh(broadbandSlowWave(:),'startbins',15,...
    'setthresh',swthresh,'diptest',false,'Schmidt',stickySW,'0Inf',true);

% % modifications by KM
% [~,~,~,~,hightheta] = bz_BimodalThresh(thratio(:),'startbins',15,...
%     'setthresh',THthresh,'diptest',false,'Schmidt',stickyTH,'0Inf',true);
hightheta = isTheta_ip == 1;


[~,~,~,~,highEMG] = bz_BimodalThresh(EMG(:),'startbins',15,...
    'setthresh',EMGthresh,'diptest',false,'Schmidt',stickyEMG,'0Inf',true);

REMtimes = (~NREMtimes & ~highEMG & hightheta);

%ACTIVE/QUIET WAKE:
WAKEtimes = ~NREMtimes & ~REMtimes;
QWAKEtimes =  WAKEtimes & ~hightheta; %Used later if QWake scored

%%
%Start/end offset due to FFT

%Construct IDX STRUCTURE FOR bz_IDXtoINT
IDX.statenames = {'WAKE','','NREM','','REM'};
IDX.timestamps = t_clus; %Timestamps pulled from clustering (in ClusterStates_GetMetrics)
IDX.states = zeros(size(IDX.timestamps));
IDX.states(NREMtimes) = 3;
IDX.states(REMtimes) = 5;
IDX.states(WAKEtimes) = 1;

if ~exist('scoreQW','var'); scoreQW = false; end

scoreQW = true;

if scoreQW
    IDX.states(QWAKEtimes) = 2;
    IDX.statenames{2} = 'QWAKE';
end

%% Fill in gaps with 0, by converting to intervals and then back to indices/timestamps
INT = bz_IDXtoINT(IDX);
IDX = bz_INTtoIDX(INT,'statenames',{'WAKE','QWAKE','NREM','','REM'});
%% Minimum Interuptions

minSWSsecs = 6;
minWnexttoREMsecs = 6;
minWinREMsecs = 6;       
minREMinWsecs = 6;
minREMsecs = 6;
minWAKEsecs = 6;

%Make the following repeated chunks of code into a single function.

%Short NREM -> WAKE
Sdur = diff(INT.NREMstate,[],2);
shortSints = Sdur<=minSWSsecs;
shortSidx = InIntervals(IDX.timestamps,INT.NREMstate(shortSints,:));
IDX.states(shortSidx) = 1;   
INT = bz_IDXtoINT(IDX);

%Short WAKE (next to REM) -> REM
Wdur = diff(INT.WAKEstate,[],2);
shortWints = Wdur<=minWnexttoREMsecs;
[~,~,WRints,~] = FindIntsNextToInts(INT.WAKEstate,INT.REMstate,1);
[~,~,~,RWints] = FindIntsNextToInts(INT.REMstate,INT.WAKEstate,1);
shortWRidx = InIntervals(IDX.timestamps,INT.WAKEstate(shortWints & (WRints | RWints),:));
IDX.states(shortWRidx) = 5;   
INT = bz_IDXtoINT(IDX);

%Short REM (in WAKE) -> WAKE
Rdur = diff(INT.REMstate,[],2);
shortRints =Rdur<=minREMinWsecs;
[~,~,~,WRints] = FindIntsNextToInts(INT.WAKEstate,INT.REMstate,1);
[~,~,RWints,~] = FindIntsNextToInts(INT.REMstate,INT.WAKEstate,1);
shortRWidx = InIntervals(IDX.timestamps,INT.REMstate(shortRints & (WRints & RWints),:));
IDX.states(shortRWidx) = 1;   
INT = bz_IDXtoINT(IDX);

%Remaining Short REM (in NREM) -> WAKE
Rdur = diff(INT.REMstate,[],2);
shortRints = Rdur<=minREMsecs;
shortRidx = InIntervals(IDX.timestamps,INT.REMstate(shortRints,:));
IDX.states(shortRidx) = 1;   
INT = bz_IDXtoINT(IDX);

%WAKE   (to SWS)     essentiall a minimum MA time
Wdur = diff(INT.WAKEstate,[],2);
shortWints = Wdur<=minWAKEsecs;
shortWidx = InIntervals(IDX.timestamps,INT.WAKEstate(shortWints,:));
IDX.states(shortWidx) = 3;   
INT = bz_IDXtoINT(IDX);

%SWS  (to NonMOV)
Sdur = diff(INT.NREMstate,[],2);
shortSints = Sdur<=minSWSsecs;
shortSidx = InIntervals(IDX.timestamps,INT.NREMstate(shortSints,:));
IDX.states(shortSidx) = 1;   
INT = bz_IDXtoINT(IDX);

ints = INT;
idx = IDX;